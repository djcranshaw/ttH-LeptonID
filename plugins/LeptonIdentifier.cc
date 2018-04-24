// -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 3; -*-
// vim: set expandtab shiftwidth=3:
//
// Package:    ttH/LeptonIdentifier
// Class:      LeptonIdentifier
//
/**\class LeptonIdentifier LeptonIdentifier.cc ttH/LeptonID/plugins/LeptonIdentifier.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Wolf
//         Created:  Mon, 27 Jul 2015 12:22:46 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TLorentzVector.h"
#include "TMVA/Reader.h"

//
// class declaration
//

enum ID { nonIsolated, preselection, fakeable, mvabased, selection };

class LeptonIdentifier : public edm::one::EDProducer<>
{
public:
   explicit LeptonIdentifier(const edm::ParameterSet &);
   ~LeptonIdentifier();

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
   virtual void beginJob() override;
   virtual void produce(edm::Event &, const edm::EventSetup &) override;
   virtual void endJob() override;

   bool passes(const pat::Electron &e, ID id);
   bool passes(const pat::Muon &mu, ID id);
   bool passes(const pat::Tau &tau, ID id);

   float mva(const pat::Muon &mu);
   float mva(const pat::Electron &ele);

   float effectiveArea(const pat::Muon&);
   float effectiveArea(const pat::Electron&);
   float effectiveAreaDR04(const pat::Muon&);
   float effectiveAreaDR04(const pat::Electron&);

   float computeRelIso04(const pat::Muon&);
   float computeRelIso04(const pat::Electron&);

   template<typename T> void addCommonUserFloats(T& lepton);
   template<typename T> void addIsolationFloats(T& lepton);
   template<class C1, class C2>
   bool matchByCommonSourceCandidatePtr(const C1 & c1, const C2 & c2);

   // ----------member data ---------------------------

   TMVA::Reader *mu_reader_;
   TMVA::Reader *ele_reader_;

   Float_t varpt;
   Float_t vareta;
   Float_t varneuRelIso;
   Float_t varchRelIso;
   Float_t varjetPtRel_in;
   Float_t varjetPtRatio_in;
   Float_t varjetBTagCSV_in;
   Float_t varjetNDauCharged_in;
   Float_t varsip3d;
   Float_t varmvaId;
   Float_t vardxy;
   Float_t vardz;
   Float_t varSegCompat;

   edm::EDGetTokenT<double> rho_token_;
   edm::EDGetTokenT<pat::PackedCandidateCollection> packedCand_token_;
   edm::EDGetTokenT<pat::ElectronCollection> ele_token_;
   edm::EDGetTokenT<pat::JetCollection> jet_token_;
   edm::EDGetTokenT<pat::MuonCollection> mu_token_;
   edm::EDGetTokenT<pat::TauCollection> tau_token_;
   edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
   edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesMapToken_;
   edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesMapToken_;

   reco::Vertex vertex_;
   pat::JetCollection jets_;

   double mu_minpt_;
   double ele_minpt_;
   double tau_minpt_;
   double loose_csv_wp; //= .46;
   double medium_csv_wp; //= .80;

   std::string jectag_;

   bool hip_safe_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
LeptonIdentifier::LeptonIdentifier(const edm::ParameterSet &config)
      : mu_minpt_(config.getParameter<double>("muonMinPt")),
        ele_minpt_(config.getParameter<double>("electronMinPt")),
        tau_minpt_(config.getParameter<double>("tauMinPt")),
        loose_csv_wp(config.getParameter<double>("LooseCSVWP")),
        medium_csv_wp(config.getParameter<double>("MediumCSVWP")),
        jectag_(config.getParameter<std::string>("JECTag")),
        hip_safe_(config.getParameter<bool>("IsHIPSafe"))
{
   produces<pat::ElectronCollection>();
   produces<pat::MuonCollection>();
   produces<pat::TauCollection>();

   rho_token_ = consumes<double>(config.getParameter<edm::InputTag>("rhoParam"));
   packedCand_token_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
   ele_token_ = consumes<pat::ElectronCollection>(config.getParameter<edm::InputTag>("electrons"));
   jet_token_ = consumes<pat::JetCollection>(config.getParameter<edm::InputTag>("jets"));
   mu_token_ = consumes<pat::MuonCollection>(config.getParameter<edm::InputTag>("muons"));
   tau_token_ = consumes<pat::TauCollection>(config.getParameter<edm::InputTag>("taus"));
   vtx_token_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));

   mvaValuesMapToken_ =
      consumes<edm::ValueMap<float>>(config.getParameter<edm::InputTag>("mvaValuesMap"));
   mvaCategoriesMapToken_ =
      consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("mvaCategoriesMap"));

   mu_reader_ = new TMVA::Reader("!Color:!Silent");
   ele_reader_ = new TMVA::Reader("!Color:!Silent");

   std::vector<TMVA::Reader *> mvas = { ele_reader_, mu_reader_ };

   for (auto &m : mvas) {
      m->AddVariable("LepGood_pt", &varpt);
      m->AddVariable("LepGood_eta", &vareta);
      m->AddVariable("LepGood_jetNDauChargedMVASel", &varjetNDauCharged_in);
      m->AddVariable("LepGood_miniRelIsoCharged", &varchRelIso);
      m->AddVariable("LepGood_miniRelIsoNeutral", &varneuRelIso);
      m->AddVariable("LepGood_jetPtRelv2", &varjetPtRel_in);
      m->AddVariable("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", &varjetBTagCSV_in);
      //m->AddVariable("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", &varjetPtRatio_in);
      m->AddVariable("LepGood_jetPtRatiov2 := (LepGood_jetBTagCSV>-5)*min(LepGood_jetPtRatiov2,1.5)+(LepGood_jetBTagCSV<-5)/(1+LepGood_relIso04)", &varjetPtRatio_in);
      m->AddVariable("LepGood_sip3d", &varsip3d);
      m->AddVariable("LepGood_dxy := log(abs(LepGood_dxy))", &vardxy);
      m->AddVariable("LepGood_dz := log(abs(LepGood_dz))", &vardz);
   }

   ele_reader_->AddVariable("LepGood_mvaIdFall17noIso", &varmvaId);
   mu_reader_->AddVariable("LepGood_segmentCompatibility", &varSegCompat);

   const std::string base = std::string(getenv("CMSSW_BASE")) + "/src/ttH/LeptonID/data";

   mu_reader_->BookMVA("BDTG method", base + "/mu_BDTG.weights.xml");
   ele_reader_->BookMVA("BDTG method", base + "/el_BDTG.weights.xml");
}

LeptonIdentifier::~LeptonIdentifier()
{
   delete ele_reader_;
   delete mu_reader_;
}

//
// member functions
//

float
LeptonIdentifier::mva(const pat::Muon &mu)
{
   varpt = mu.pt();
   vareta = mu.eta();
   varchRelIso = mu.userFloat("miniAbsIsoCharged") / mu.pt();
   varneuRelIso = mu.userFloat("miniAbsIsoNeutralcorr") / mu.pt();
   varjetPtRel_in = mu.userFloat("nearestJetPtRel");
   varjetPtRatio_in = (mu.userFloat("nearestJetCsv") > -99.) ?
      std::min(mu.userFloat("nearestJetPtRatio"), 1.5f) :  // found matched jet
      (1./(1.+mu.userFloat("relIso04")));  // no matched jet
   varjetBTagCSV_in = std::max(mu.userFloat("nearestJetCsv"), 0.f);
   varjetNDauCharged_in = mu.userFloat("nearestJetNDauCharged");
   varsip3d = mu.userFloat("sip3D");
   vardxy = log(fabs(mu.userFloat("dxy")));
   vardz = log(fabs(mu.userFloat("dz")));
   varSegCompat = mu.segmentCompatibility();

   return mu_reader_->EvaluateMVA("BDTG method");
}

float
LeptonIdentifier::mva(const pat::Electron &ele)
{
   varpt = ele.pt();
   vareta = ele.eta();
   varchRelIso = ele.userFloat("miniAbsIsoCharged") / ele.pt();
   varneuRelIso = ele.userFloat("miniAbsIsoNeutralcorr") / ele.pt();
   varjetPtRel_in = ele.userFloat("nearestJetPtRel");
   varjetPtRatio_in = (ele.userFloat("nearestJetCsv") > -99.) ?
      std::min(ele.userFloat("nearestJetPtRatio"), 1.5f) :  // found matched jet
      (1./(1.+ele.userFloat("relIso04")));  // no matched jets
   varjetBTagCSV_in = std::max(ele.userFloat("nearestJetCsv"), 0.f);
   varjetNDauCharged_in = ele.userFloat("nearestJetNDauCharged");
   varsip3d = ele.userFloat("sip3D");
   vardxy = log(fabs(ele.userFloat("dxy")));
   vardz = log(fabs(ele.userFloat("dz")));
   varmvaId = ele.userFloat("eleMvaId");

   return ele_reader_->EvaluateMVA("BDTG method");
}

// ------------ id functions ------------

// intermediate function for CHEP
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Short_Term_Medium_Muon_Definitio
bool isMediumMuon(const reco::Muon & recoMu, bool isHIPSafe)
{
   std::cout << "Method deprecated!! Use isMediumMuon() directly from miniAOD instead." << std::endl;
   assert(0);
   
   bool goodGlob = recoMu.isGlobalMuon() &&
      recoMu.globalTrack()->normalizedChi2() < 3 &&
      recoMu.combinedQuality().chi2LocalPosition < 12 &&
      recoMu.combinedQuality().trkKink < 20;
   bool isMedium = muon::isLooseMuon(recoMu) &&
      recoMu.innerTrack()->validFraction() > (isHIPSafe ? 0.49 : 0.80) &&
      muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451);
   return isMedium;
}

bool
LeptonIdentifier::passes(const pat::Muon &mu, ID id)
{
   //double minMuonPt = id == preselection ? 5.0 : 15.0; // iMinPt;
   double minMuonPt = 5.0;
   double minMuonConePt = 10.;
   double maxMuonEta = 2.4;

   bool passesMuonBestTrackID = false;
   if (mu.innerTrack().isAvailable()) { // innerTrack() // muonBestTrack // isAvailable
      passesMuonBestTrackID = (fabs(mu.userFloat("dxy")) < 0.05 && fabs(mu.userFloat("dz")) < 0.1 && (mu.userFloat("sip3D") < 8));
   }

   bool passesKinematics = (mu.pt() > minMuonPt) and (fabs(mu.eta()) < maxMuonEta);
   bool passesIso = (mu.userFloat("miniRelIso") < 0.4);
   bool passesPreselection = mu.isLooseMuon() && passesMuonBestTrackID;

   bool passesID = false;

   switch (id) {
      case preselection:
         passesID = passesPreselection;
         break;
      case fakeable:
         if (mu.userFloat("leptonMVA") > 0.90)
            passesID = passesPreselection and mu.userFloat("nearestJetCsv") < medium_csv_wp;
         else
            passesID = passesPreselection and mu.userFloat("nearestJetCsv") < 0.3 and mu.userFloat("nearestJetPtRatio") > 0.5 and mu.segmentCompatibility() > 0.3;
         //passesIso = true;
         passesKinematics =
            passesKinematics and (mu.userFloat("correctedPt") > minMuonConePt);
         break;
      case mvabased:
         //passesIso = true;
         passesID = passesPreselection and
            (mu.userFloat("leptonMVA") > 0.90) and
            (mu.userFloat("nearestJetCsv") < medium_csv_wp) and
            mu.userFloat("isMediumMuon");
            //isMediumMuon(mu, hip_safe_);
         passesKinematics =
            passesKinematics and (mu.userFloat("correctedPt") > minMuonConePt);
         break;
      case nonIsolated:
         edm::LogError("LeptonID") << "Invalid ID 'nonIsolated' for muons!";
         return false;
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID);
}

bool
LeptonIdentifier::passes(const pat::Electron &ele, ID id)
{
   double minElectronPt = 7.0;
   double minElectronConePt = 10.;
   
   bool passesKinematics = (ele.pt() > minElectronPt) and (fabs(ele.eta()) < 2.5);
   bool passesIso = ele.userFloat("miniRelIso") < 0.4;
   
   bool passGsfTrackID = false;
   if (ele.gsfTrack().isAvailable())
      passGsfTrackID = fabs(ele.userFloat("dxy")) < 0.05 and fabs(ele.userFloat("dz")) < 0.1 and ele.userFloat("numMissingHits") <= 1;

   float scEta = fabs(ele.userFloat("superClusterEta"));
   double eleMva = ele.userFloat("eleMvaId");
   // Loose Fall17noIso
   bool passesMvaWP = false;
   if (ele.pt() <= 10.) {
      if (scEta < 0.8)
         passesMvaWP = eleMva > -0.13285867293779202;
      else if (scEta < 1.479)
         passesMvaWP = eleMva > -0.31765300958836074;
      else
         passesMvaWP = eleMva > -0.0799205914718861;
   }
   else {
      if (scEta < 0.8)
         passesMvaWP = eleMva > -0.856871961305474;
      else if (scEta < 1.479)
         passesMvaWP = eleMva > -0.8107642141584835;
      else
         passesMvaWP = eleMva > -0.7179265933023059;
   }

   bool passesPreselection = passGsfTrackID and ele.userFloat("sip3D") < 8 and passesMvaWP;

   bool passesCuts = false;
   if (fabs(ele.eta()) < 0.8) {
      passesCuts = ele.full5x5_sigmaIetaIeta() < 0.011 &&
         ele.hcalOverEcal() < 0.10 &&
         fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.01 &&
         fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.04 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() > -0.05 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() < 0.010;
   }
   else if (fabs(ele.eta()) < 1.479) {
      passesCuts = ele.full5x5_sigmaIetaIeta() < 0.011 &&
         ele.hcalOverEcal() < 0.10 &&
         fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.01 &&
         fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.04 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() > -0.05 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() < 0.010;
   }
   else if (fabs(ele.eta()) < 2.5) {
      passesCuts = ele.full5x5_sigmaIetaIeta() < 0.030 &&
         ele.hcalOverEcal() < 0.07 &&
         fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.008 &&
         fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.07 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() > -0.05 &&
         1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() < 0.005;
   }


   bool passesID = false;
   bool passesJetCSV = false;

   switch (id) {
      case preselection:
         passesID = passesPreselection;
         break;
      case fakeable:
         if (ele.userFloat("leptonMVA") > 0.90)
            passesJetCSV = ele.userFloat("nearestJetCsv") < medium_csv_wp;
         else
            passesJetCSV = ele.userFloat("nearestJetCsv") < 0.3 && ele.userFloat("nearestJetPtRatio") > 0.5;
         
         passesID = passesPreselection and
                    passesCuts and
                    passesJetCSV and
                    ele.userFloat("numMissingHits") == 0;
         passesKinematics = passesKinematics and
                            (ele.userFloat("correctedPt") > minElectronConePt);
         break;
      case mvabased:
         passesID = passesPreselection and
                    passesCuts and
                    ele.userFloat("leptonMVA") > 0.90 and
                    ele.userFloat("nearestJetCsv") < medium_csv_wp and
                    ele.userFloat("numMissingHits") == 0 and
                    ele.passConversionVeto();
         passesKinematics = passesKinematics and
                            (ele.userFloat("correctedPt") > minElectronConePt);
         break;
      case nonIsolated:
         edm::LogError("LeptonID") << "Invalid ID 'nonIsolated' for electrons!";
         return false;
         break;
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID);
}

bool
LeptonIdentifier::passes(const pat::Tau &tau, ID id)
{
   double minTauPt = 5.; // iMinPt;

   bool passesIso = false;
   bool passesID = false;

   bool passesKinematics = ((tau.pt() > minTauPt) && (fabs(tau.eta()) < 2.3));
   bool passesPVassoc = ((fabs(tau.userFloat("dxy")) < 1000) && (fabs(tau.userFloat("dz")) < 0.2));


   switch (id) {
      case nonIsolated:
         passesIso = true;
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      case preselection:
         //passesIso = (tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT") > 0.5);
         //passesIso = (tau.tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017") > 0.5);
         passesIso = (tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017") > 0.5);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      case selection:
         //passesIso = (tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT") > 0.5);
         //passesIso = (tau.tauID("byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017") > 0.5);
         passesIso = (tau.tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017") > 0.5);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID);
}

template<class C1, class C2> bool
LeptonIdentifier::matchByCommonSourceCandidatePtr(const C1 & c1, const C2 & c2) {
    for(unsigned int i1 = 0 ; i1 < c1.numberOfSourceCandidatePtrs();i1++){
        auto  c1s=c1.sourceCandidatePtr(i1);
            for(unsigned int i2 = 0 ; i2 < c2.numberOfSourceCandidatePtrs();i2++) {
                if(c2.sourceCandidatePtr(i2)==c1s) return true;
            }
    }
    return false;
}

template<typename T> void
LeptonIdentifier::addCommonUserFloats(T& lepton)
{
   double corr_factor = 1.;
   double L1_SF = 1.;
   pat::Jet matchedJet;
   pat::Jet matchedJetL1;
   pat::Jet matchedJetRaw;

   // match lep and jet by jet constituents
   bool foundmatch = false;
   for (const auto &j : jets_) {
      if (matchByCommonSourceCandidatePtr(lepton, j)) {
         matchedJet = j;
         matchedJetL1 = j;
         matchedJetRaw = j;

         matchedJetL1.setP4(j.correctedJet("L1FastJet", "none", jectag_).p4());
         matchedJetRaw.setP4(j.correctedJet(0).p4());
         corr_factor = j.p4().E() / j.correctedJet(0).p4().E();
         L1_SF = matchedJetL1.p4().E() / matchedJetRaw.p4().E();

         foundmatch = true;
        
         //assert(fabs(L1_SF-j.jecFactor("L1FastJet")/j.jecFactor("Uncorrected"))<1e-4);
         break;
      }
   }

   float njet_csv = -100.;
   float njet_pt_ratio = 1.;
   float njet_pt_rel = 0.;
   float njet_ndau_charged = 0.;

   if (foundmatch) {
      
      njet_csv = matchedJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

      if ((matchedJet.correctedJet(0).p4() - lepton.p4()).Rho() >= 1e-4) {
         for (unsigned int i = 0, n = matchedJet.numberOfSourceCandidatePtrs(); i < n; ++i) {
            
            const pat::PackedCandidate &dau_jet = dynamic_cast<const pat::PackedCandidate &>(*(matchedJet.sourceCandidatePtr(i)));
            //float dR = reco::deltaR(matchedJet.eta(),matchedJet.phi(),
            //                        dau_jet.p4().eta(), dau_jet.p4().phi());
            // To match nanoAOD
            float dR = reco::deltaR(lepton.eta(),lepton.phi(),
                                    dau_jet.p4().eta(), dau_jet.p4().phi());
            
            bool isgoodtrk = false;
            try {
               const reco::Track trk = dau_jet.pseudoTrack();
               //const math::XYZPoint vtx_position = lepton.vertex();
               const auto vtx_position = vertex_.position();
               
               if(trk.pt()>1 &&
                  trk.hitPattern().numberOfValidHits()>=8 &&
                  trk.hitPattern().numberOfValidPixelHits()>=2 &&
                  trk.normalizedChi2()<5 &&
                  std::fabs(trk.dxy(vtx_position))<0.2 &&
                  std::fabs(trk.dz(vtx_position))<17
                  ) isgoodtrk = true;
               
               if( dR<=0.4 && dau_jet.charge()!=0 && dau_jet.fromPV()>1 && isgoodtrk) njet_ndau_charged++;
            } catch(...){}
            
         }
      
         auto lepAwareJetp4 = (matchedJet.p4() * (1. / corr_factor) - lepton.p4() * (1. / L1_SF)) * corr_factor + lepton.p4(); // "lep-aware" JEC
         if ((matchedJet.p4() * (1. / corr_factor) - lepton.p4()).Rho() < 1e-4)
            lepAwareJetp4 = lepton.p4();

         njet_pt_ratio = std::min(lepton.pt() / lepAwareJetp4.pt(), 1.5);

         TLorentzVector l4 = TLorentzVector(lepton.p4().Px(), lepton.p4().Py(), lepton.p4().Pz(), lepton.p4().E());
         TLorentzVector j4 = TLorentzVector(lepAwareJetp4.Px(), lepAwareJetp4.Py(), lepAwareJetp4.Pz(), lepAwareJetp4.E());

         if ((j4 - l4).Rho() < 1e-4) njet_pt_rel = 0.;
         else njet_pt_rel = l4.Perp((j4 - l4).Vect());
      }
   }
   
   lepton.addUserFloat("nearestJetCsv", njet_csv);
   lepton.addUserFloat("nearestJetPtRatio", njet_pt_ratio);
   lepton.addUserFloat("nearestJetPtRel", njet_pt_rel);
   lepton.addUserFloat("nearestJetNDauCharged", njet_ndau_charged);
   auto mva_value = mva(lepton);
   lepton.addUserFloat("leptonMVA", mva_value);

   if (  (abs(lepton.pdgId()) != 13 || lepton.userFloat("isMediumMuon")) && lepton.userFloat("leptonMVA") > 0.90 )
      {
         lepton.addUserFloat("correctedPt", lepton.pt());
      }
   else
      {
         lepton.addUserFloat("correctedPt", .90 * lepton.pt() / njet_pt_ratio);
      }
      
   lepton.addUserFloat("idPreselection", passes(lepton, preselection));

   if (lepton.userFloat("idPreselection") > .5) {
      lepton.addUserFloat("idFakeable", passes(lepton, fakeable));
      lepton.addUserFloat("idMVABased", passes(lepton, mvabased));
   } else {
      lepton.addUserFloat("idFakeable", 0.);
      lepton.addUserFloat("idMVABased", 0.);
   }
}

template<typename T> void
LeptonIdentifier::addIsolationFloats(T& lepton)
{
   // dR = 0.3
   lepton.addUserFloat("miniAbsIsoCharged",
                       lepton.miniPFIsolation().chargedHadronIso());  
   lepton.addUserFloat("miniAbsIsoNeutral",
                       lepton.miniPFIsolation().neutralHadronIso() +
                       lepton.miniPFIsolation().photonIso());
   
   // PU correction
   lepton.addUserFloat("miniIsoR", std::max(std::min(10.0/lepton.pt(), 0.2), 0.05) ); 
   lepton.addUserFloat("effArea", effectiveArea(lepton));
   assert(lepton.hasUserFloat("rho"));
   float PUCorrection =
      lepton.userFloat("rho") * lepton.userFloat("effArea") * 
      (lepton.userFloat("miniIsoR")/0.3) * (lepton.userFloat("miniIsoR")/0.3);
   
   lepton.addUserFloat("miniAbsIsoNeutralcorr",
                       std::max(lepton.userFloat("miniAbsIsoNeutral") - 
                                PUCorrection, 0.f));
   lepton.addUserFloat("miniRelIso", (lepton.userFloat("miniAbsIsoCharged")+
                                      lepton.userFloat("miniAbsIsoNeutralcorr")
                                      ) / lepton.pt());


   lepton.addUserFloat("relIso04", computeRelIso04(lepton));
}

float LeptonIdentifier::computeRelIso04(const pat::Muon& mu)
{
   float absIso04Charged = mu.pfIsolationR04().sumChargedHadronPt;
   float absIso04Neutral =
      mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt;

   // deltaBeta PU correction for muon
   float absIso04NeutralCorr =
      std::max(absIso04Neutral - mu.pfIsolationR04().sumPUPt/2, 0.f);

   return (absIso04Charged + absIso04NeutralCorr) / mu.pt();
}

float LeptonIdentifier::computeRelIso04(const pat::Electron& ele)
{
   // These are for dR = 0.3
   //float absIso04Charged = ele.pfIsolationVariables().sumChargedHadronPt;
   //float absIso04Neutral = ele.pfIsolationVariables().sumNeutralHadronEt +
   //   ele.pfIsolationVariables().sumPhotonEt;

   // dR = 0.4
   float absIso04Charged = ele.chargedHadronIso();
   float absIso04Neutral = ele.neutralHadronIso() + ele.photonIso();

   // effective area PU correction for electron
   assert(ele.hasUserFloat("rho"));  
   float absIso04NeutralCorr =
      std::max(absIso04Neutral - ele.userFloat("rho")*effectiveAreaDR04(ele), 0.f);

   return (absIso04Charged + absIso04NeutralCorr) / ele.pt();
}

float LeptonIdentifier::effectiveArea(const pat::Muon& mu)
{
   float eta = mu.eta();
   
   // Fall17 EA dR=0.3
   // have to be used with rho = fixedGridRhoFastjetAll
   //assert(fabs(eta) <= 2.5);
   
   if (fabs(eta) < 0.8)
      return 0.0566;
   else if (fabs(eta) < 1.3)
      return 0.0562;
   else if (fabs(eta) < 2.0)
      return 0.0363;
   else if (fabs(eta) < 2.2)
      return 0.0119;
   else if (fabs(eta) <= 2.5)
      return 0.0064;
   else
      return 0.;
}

float LeptonIdentifier::effectiveAreaDR04(const pat::Muon& mu)
{
   return effectiveArea(mu) * 16./9.;
}

float LeptonIdentifier::effectiveArea(const pat::Electron& ele)
{
   float eta = ele.eta();
   
   // Fall17 EA dR=0.3
   // have to be used with rho = fixedGridRhoFastjetAll
   //assert(fabs(eta) <= 2.5);

   if (fabs(eta) < 1.0)
      return 0.1566;
   else if (fabs(eta) < 1.479)
      return 0.1626;
   else if (fabs(eta) < 2.0)
      return 0.1073;
   else if (fabs(eta) < 2.2)
      return 0.0854;
   else if (fabs(eta) < 2.3)
      return 0.1051;
   else if (fabs(eta) < 2.4)
      return 0.1204;
   else if (fabs(eta) <= 2.5)
      return 0.1524;
   else
      return 0.;
}

float LeptonIdentifier::effectiveAreaDR04(const pat::Electron& ele)
{
   return effectiveArea(ele) * 16./9.;
}

// ------------ method called to produce the data  ------------
void
LeptonIdentifier::produce(edm::Event &event, const edm::EventSetup &setup)
{
   std::unique_ptr<pat::ElectronCollection> eles(new pat::ElectronCollection());
   std::unique_ptr<pat::MuonCollection> mus(new pat::MuonCollection());
   std::unique_ptr<pat::TauCollection> taus(new pat::TauCollection());

   edm::Handle<double> rho;
   edm::Handle<pat::PackedCandidateCollection> packedCands;
   edm::Handle<pat::ElectronCollection> input_ele;
   edm::Handle<pat::JetCollection> input_jet;
   edm::Handle<pat::MuonCollection> input_mu;
   edm::Handle<pat::TauCollection> input_tau;
   edm::Handle<reco::VertexCollection> input_vtx;

   edm::Handle<edm::ValueMap<float>> mvaValues;
   edm::Handle<edm::ValueMap<int>> mvaCategories;

   event.getByToken(rho_token_, rho);
   event.getByToken(packedCand_token_, packedCands);
   //event.getByToken(ele_token_, input_ele_raw);
   event.getByToken(ele_token_, input_ele);
   event.getByToken(jet_token_, input_jet);
   event.getByToken(mu_token_, input_mu);
   event.getByToken(tau_token_, input_tau);
   event.getByToken(vtx_token_, input_vtx);
   event.getByToken(mvaValuesMapToken_, mvaValues);
   event.getByToken(mvaCategoriesMapToken_, mvaCategories);

   const edm::ValueMap<float> ele_mvaValues = (*mvaValues.product());

   // determine primary vertex
   bool valid = false;
   for (const auto &v : *input_vtx) {
      if (!v.isFake() && v.ndof() >= 4 && abs(v.z()) <= 24. && abs(v.position().Rho()) <= 2.) {
         vertex_ = v;
         valid = true;
         break;
      }
   }

   if (not valid) {
      event.put(std::move(eles));
      event.put(std::move(mus));
      event.put(std::move(taus));
      return;
   }

   jets_ = *input_jet;

   for (auto mu : *input_mu) {
      if (mu.pt() < mu_minpt_)
         continue;

      // add members
      if (mu.innerTrack().isAvailable()) // muonBestTrack
      {
         mu.addUserFloat("dxy", mu.innerTrack()->dxy(vertex_.position()));
         mu.addUserFloat("dz", mu.innerTrack()->dz(vertex_.position()));
         mu.addUserFloat("numValidPixelHits", mu.innerTrack()->hitPattern().numberOfValidPixelHits());
         mu.addUserFloat("trackerLayersWithMeasurement", mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
         mu.addUserFloat("chargeFlip", mu.innerTrack()->ptError() / mu.innerTrack()->pt());
         mu.addUserFloat("validFraction", mu.innerTrack()->validFraction());
      } else {
         mu.addUserFloat("dxy", -666.);
         mu.addUserFloat("dz", -666.);
         mu.addUserFloat("numValidPixelHits", -666.);
         mu.addUserFloat("trackerLayersWithMeasurement", -666.);
         mu.addUserFloat("chargeFlip", -666.);
         mu.addUserFloat("validFraction", -666.);
      }

      if (mu.globalTrack().isAvailable()) {
         mu.addUserFloat("normalizedChiSq", mu.globalTrack()->normalizedChi2());
         mu.addUserFloat("numValidMuonHits", mu.globalTrack()->hitPattern().numberOfValidMuonHits());
      } else {
         mu.addUserFloat("normalizedChiSq", -666.);
         mu.addUserFloat("numValidMuonHits", -666.);
      }

      mu.addUserFloat("isMediumMuon", mu.isMediumMuon());//isMediumMuon(mu, hip_safe_));
      mu.addUserFloat("rho", *rho);

      addIsolationFloats(mu);

      mu.addUserFloat("localChiSq", mu.combinedQuality().chi2LocalPosition);
      mu.addUserFloat("trackKink", mu.combinedQuality().trkKink);     
      mu.addUserFloat("sip3D", fabs(mu.dB(pat::Muon::PV3D) / mu.edB(pat::Muon::PV3D)));

      addCommonUserFloats(mu);

      mus->push_back(mu);
   }

   int ele_index_for_mva = 0;
   for (auto ele : *input_ele) {
      ele_index_for_mva++;
      
      if (ele.pt() < ele_minpt_)
         continue;

      ele.addUserFloat("rho", *rho);

      addIsolationFloats(ele);

      ele.addUserFloat("superClusterEta", abs(ele.superCluster()->position().eta()));
      ele.addUserFloat("dxy", ele.gsfTrack()->dxy(vertex_.position()));
      ele.addUserFloat("dz", ele.gsfTrack()->dz(vertex_.position()));

      //ele.addUserFloat("numMissingHits", ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
      ele.addUserFloat("numMissingHits", ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      ele.addUserFloat("sip3D", fabs(ele.dB(pat::Electron::PV3D) / ele.edB(pat::Electron::PV3D)));
      ele.addUserFloat("eleMvaId", ele_mvaValues.get(ele_index_for_mva - 1));
      ele.addUserFloat("isMediumMuon", 0.);

      addCommonUserFloats(ele);

      eles->push_back(ele);
   }


   for (auto tau : *input_tau) {
      if (tau.pt() < tau_minpt_)
         continue;

      //       pat::Jet matchedJet;
      //       double dR = 666.;
      //       for (const auto& j: jets_) {
      // 	double newDR = helper_.DeltaR(&j, &tau);
      // 	if (newDR < dR) {
      // 	  dR = newDR;
      // 	  matchedJet = j;
      // 	}
      //       }

      float dxy_old = -666.;
      float dz_old = -666.;
      float dxy = -666.;
      float dz = -666.;
      float id_non_isolated = -666.;
      float id_preselection = -666.;
      float id_selection = -666.;

      if (tau.leadChargedHadrCand().isAvailable()) {
         auto track = tau.leadChargedHadrCand()->bestTrack();

         if (track) {
            dxy_old = track->dxy(vertex_.position());
            dz_old = track->dz(vertex_.position());
         }
      }

      // As in the SM Htautau analysis
      // still need to understand why dz != dz_old
      auto packedLeadTauCand =
         dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      dz = packedLeadTauCand->dz();
      dxy = packedLeadTauCand->dxy();

      tau.addUserFloat("dxy_old", dxy_old);
      tau.addUserFloat("dz_old", dz_old);
      tau.addUserFloat("dxy", dxy);
      tau.addUserFloat("dz", dz);
      
      id_non_isolated = passes(tau, nonIsolated);
      id_preselection = passes(tau, preselection);
      id_selection = passes(tau, selection);

      tau.addUserFloat("idNonIsolated", id_non_isolated);
      tau.addUserFloat("idPreselection", id_preselection);
      tau.addUserFloat("idSelection", id_selection);

      taus->push_back(tau);
   }
   event.put(std::move(eles));
   event.put(std::move(mus));
   event.put(std::move(taus));
}

// ------------ method called once each job just before starting event loop  ------------
void
LeptonIdentifier::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
LeptonIdentifier::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
LeptonIdentifier::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
LeptonIdentifier::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeptonIdentifier::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeptonIdentifier::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonIdentifier::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   // The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(LeptonIdentifier);
