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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

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

   template<typename T> void addCommonUserFloats(T& lepton);

   // ----------member data ---------------------------
   MiniAODHelper helper_;

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
   edm::EDGetTokenT<edm::View<pat::Electron>> ele_token_;
   edm::EDGetTokenT<pat::JetCollection> jet_token_;
   edm::EDGetTokenT<pat::MuonCollection> mu_token_;
   edm::EDGetTokenT<pat::TauCollection> tau_token_;
   edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
   edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesHZZMapToken_;
   edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesHZZMapToken_;
   // edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesGPMapToken_;
   // edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesGPMapToken_;

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
   ele_token_ = consumes<edm::View<pat::Electron>>(config.getParameter<edm::InputTag>("electrons"));
   jet_token_ = consumes<pat::JetCollection>(config.getParameter<edm::InputTag>("jets"));
   mu_token_ = consumes<pat::MuonCollection>(config.getParameter<edm::InputTag>("muons"));
   tau_token_ = consumes<pat::TauCollection>(config.getParameter<edm::InputTag>("taus"));
   vtx_token_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));

   mvaValuesHZZMapToken_ =
      consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"));
   mvaCategoriesHZZMapToken_ =
      consumes<edm::ValueMap<int>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"));

   // mvaValuesGPMapToken_ =
   //    consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"));
   // mvaCategoriesGPMapToken_ =
   //    consumes<edm::ValueMap<int>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"));

   // Who gives a FUCK about these parameters?  They are not used in the
   // methods we access, which could be spun off, anyways.
   helper_.SetUp("2015_74x", 666, analysisType::LJ, false);

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
      m->AddVariable("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", &varjetPtRatio_in);
      m->AddVariable("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", &varjetBTagCSV_in);
      m->AddVariable("LepGood_sip3d", &varsip3d);
      m->AddVariable("LepGood_dxy := log(abs(LepGood_dxy))", &vardxy);
      m->AddVariable("LepGood_dz := log(abs(LepGood_dz))", &vardz);
   }

   ele_reader_->AddVariable("LepGood_mvaIdSpring16HZZ", &varmvaId);
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
   varjetPtRatio_in = std::min(mu.userFloat("nearestJetPtRatio"), 1.5f);
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
   varjetPtRatio_in = std::min(ele.userFloat("nearestJetPtRatio"), 1.5f);
   varjetBTagCSV_in = std::max(ele.userFloat("nearestJetCsv"), 0.f);
   varjetNDauCharged_in = ele.userFloat("nearestJetNDauCharged");
   varsip3d = ele.userFloat("sip3D");
   vardxy = log(fabs(ele.userFloat("dxy")));
   vardz = log(fabs(ele.userFloat("dz")));
   varmvaId = ele.userFloat("eleMvaIdHZZ");

   return ele_reader_->EvaluateMVA("BDTG method");
}

// ------------ id functions ------------

// intermediate function for CHEP
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Short_Term_Medium_Muon_Definitio
bool isMediumMuon(const reco::Muon & recoMu, bool isHIPSafe)
{
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
   double minMuonPt = id == preselection ? 5.0 : 15.0; // iMinPt;
   double maxMuonEta = 2.4;

   bool passesMuonBestTrackID = false;
   if (mu.innerTrack().isAvailable()) { // innerTrack() // muonBestTrack // isAvailable
      passesMuonBestTrackID = (fabs(mu.userFloat("dxy")) < 0.05 && fabs(mu.userFloat("dz")) < 0.1 && (mu.userFloat("sip3D") < 8));
   }

   bool passesKinematics = (mu.pt() > minMuonPt) and (fabs(mu.eta()) < maxMuonEta);
   bool passesIso = (mu.userFloat("miniIso") < 0.4);
   bool passesPreselection = mu.isLooseMuon() && passesMuonBestTrackID;

   bool passesID = false;

   float corrected_pt = mu.pt();
   if (mu.userFloat("leptonMVA") < 0.90) {
      corrected_pt = 0.90 * corrected_pt / mu.userFloat("nearestJetPtRatio");
      //passesKinematics = (corrected_pt > minMuonPt) and (fabs(mu.eta()) < 2.5);
   }

   switch (id) {
      case preselection:
         passesID = passesPreselection;
         break;
      case fakeable:
         if (mu.userFloat("leptonMVA") > 0.90)
            passesID = passesPreselection and mu.userFloat("nearestJetCsv") < medium_csv_wp;
         else
            passesID = passesPreselection and mu.userFloat("nearestJetCsv") < 0.3 and mu.userFloat("nearestJetPtRatio") > 0.5 and mu.segmentCompatibility() > 0.3;
         passesIso = true;
         break;
      case mvabased:
         passesIso = true;
         passesID = passesPreselection and
            mu.userFloat("leptonMVA") > 0.90 and
            mu.userFloat("nearestJetCsv") < medium_csv_wp and
            mu.userFloat("isMediumMuon");
            //isMediumMuon(mu, hip_safe_);
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
   double minElectronPt = id == preselection ? 7.0 : 15.0; // iMinPt;
   bool passesKinematics = (ele.pt() > minElectronPt) and (fabs(ele.eta()) < 2.5);
   bool passesIso = ele.userFloat("miniIso") < 0.4;

   bool passesMVA = false;
   //double eleMvaGP = ele.userFloat("eleMvaId");
   double eleMvaHZZ = ele.userFloat("eleMvaIdHZZ");
   float scEta = ele.userFloat("superClusterEta");

   // VLooseIdEmu WP
   // c.f. https://github.com/CERN-PH-CMG/cmg-cmssw/blob/a76bc9fb439b6238af56649961b3952dbef95626/PhysicsTools/Heppy/python/physicsobjects/Electron.py#L361
   // double _vlow[3] = {-0.30,-0.46,-0.63};
   // double _low[3] = {-0.86,-0.85,-0.81};
   // double _high[3] = {0.,0.,0.7};

   // if (ele.pt() <= 10) {
   //    // use HZZ MVA if electron pt less than 10 GeV
   //    passesMVA = eleMvaHZZ > _high[(scEta>=0.8)+(scEta>=1.479)];
      
   // }
   // else {
   //    // _low below 15 GeV, _high above 25 GeV, interpolation in between
   //    passesMVA = eleMvaGP > _high[(scEta>=0.8)+(scEta>=1.479)];
   // }
   
   if (scEta<1.479) passesMVA = eleMvaHZZ > 0.0;
   else passesMVA = eleMvaHZZ > 0.7;



   bool passGsfTrackID = false;
   if (ele.gsfTrack().isAvailable())
      passGsfTrackID = fabs(ele.userFloat("dxy")) < 0.05 and fabs(ele.userFloat("dz")) < 0.1 and ele.userFloat("numMissingHits") <= 1;

   float corrected_pt = ele.pt();
   if (ele.userFloat("leptonMVA") < 0.90) {
      corrected_pt = 0.90 * corrected_pt / ele.userFloat("nearestJetPtRatio");
      //passesKinematics = (corrected_pt > minElectronPt) and (fabs(ele.eta()) < 2.5);
   }

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

   bool passesPreselection = passesKinematics and passesIso and passesMVA and passGsfTrackID and ele.userFloat("sip3D") < 8;

   bool passesID = false;
   bool passesJetCSV = false;

   bool passesObjectSelection = true;

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
                    passesJetCSV;
         break;
      case mvabased:
         passesID = passesPreselection and
                    passesCuts and
                    ele.userFloat("leptonMVA") > 0.90 and
                    ele.userFloat("nearestJetCsv") < medium_csv_wp;
         break;
      case nonIsolated:
         edm::LogError("LeptonID") << "Invalid ID 'nonIsolated' for electrons!";
         return false;
         break;
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID && passesObjectSelection);
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
         passesIso = (tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT") > 0.5);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      case selection:
         passesIso = (tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT") > 0.5);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID);
}

template<typename T> void
LeptonIdentifier::addCommonUserFloats(T& lepton)
{
   double corr_factor = 1.;
   double L1_SF = 1.;
   pat::Jet matchedJet;
   pat::Jet matchedJetL1;
   pat::Jet matchedJetRaw;
   double dR = 666.;

   for (const auto &j : jets_) {
      double newDR = helper_.DeltaR(&j, &lepton);
      if (newDR < dR and newDR < .4) {
         dR = newDR;
         matchedJet = j;
         matchedJetL1 = j;
         matchedJetRaw = j;

         matchedJetL1.setP4(j.correctedJet("L1FastJet", "none", jectag_).p4());
         matchedJetRaw.setP4(j.correctedJet(0).p4());
         corr_factor = j.p4().E() / j.correctedJet(0).p4().E();
         L1_SF = matchedJetL1.p4().E() / matchedJetRaw.p4().E();
      }
   }

   lepton.addUserFloat("nearestJetDr", min(dR, 0.5)); // no longer used in MVA

   float njet_csv = 0;
   float njet_pt_ratio = 1.;
   float njet_pt_rel = 0.;
   float njet_ndau_charged = 0.;

   if (jets_.size() > 0 and dR < .4) {
      njet_csv = matchedJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      if (njet_csv < 0)
         njet_csv = -10.;

      for (unsigned int i = 0, n = matchedJet.numberOfSourceCandidatePtrs(); i < n; ++i) {

         const pat::PackedCandidate &dau_jet = dynamic_cast<const pat::PackedCandidate &>(*(matchedJet.sourceCandidatePtr(i)));
         float dR = helper_.DeltaR(&matchedJet, &(dau_jet.p4()));
         
         bool isgoodtrk = false;
         try {
            const reco::Track trk = dau_jet.pseudoTrack();
            const math::XYZPoint vtx_position = lepton.vertex();
         
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


      if (dR <= 0.4 and (matchedJet.correctedJet(0).p4() - lepton.p4()).Rho() >= 1e-4) {
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
   lepton.addUserFloat("idPreselection", passes(lepton, preselection));


   if (  (abs(lepton.pdgId()) != 13 || lepton.userFloat("isMediumMuon")) && lepton.userFloat("leptonMVA") > 0.90 )
      {
         lepton.addUserFloat("correctedPt", lepton.pt());
      }
   else
      {
         lepton.addUserFloat("correctedPt", .90 * lepton.pt() / njet_pt_ratio);
      }
   
 




   if (lepton.userFloat("idPreselection") > .5) {
      lepton.addUserFloat("idFakeable", passes(lepton, fakeable));
      lepton.addUserFloat("idMVABased", passes(lepton, mvabased));
   } else {
      lepton.addUserFloat("idFakeable", -666.);
      lepton.addUserFloat("idMVABased", -666.);
   }
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
   edm::Handle<edm::View<pat::Electron>> input_ele_raw;
   edm::Handle<pat::JetCollection> input_jet;
   edm::Handle<pat::MuonCollection> input_mu;
   edm::Handle<pat::TauCollection> input_tau;
   edm::Handle<reco::VertexCollection> input_vtx;

   edm::Handle<edm::ValueMap<float>> mvaValuesHZZ;
   edm::Handle<edm::ValueMap<int>> mvaCategoriesHZZ;

   //edm::Handle<edm::ValueMap<float>> mvaValuesGP;
   //edm::Handle<edm::ValueMap<int>> mvaCategoriesGP;

   event.getByToken(rho_token_, rho);
   event.getByToken(packedCand_token_, packedCands);
   event.getByToken(ele_token_, input_ele_raw);
   event.getByToken(jet_token_, input_jet);
   event.getByToken(mu_token_, input_mu);
   event.getByToken(tau_token_, input_tau);
   event.getByToken(vtx_token_, input_vtx);
   event.getByToken(mvaValuesHZZMapToken_, mvaValuesHZZ);
   event.getByToken(mvaCategoriesHZZMapToken_, mvaCategoriesHZZ);
   //event.getByToken(mvaValuesGPMapToken_, mvaValuesGP);
   //event.getByToken(mvaCategoriesGPMapToken_, mvaCategoriesGP);

   //   const edm::ValueMap<float> ele_mvaValuesGP = (*mvaValuesGP.product());
   const edm::ValueMap<float> ele_mvaValuesHZZ = (*mvaValuesHZZ.product());

   helper_.SetRho(*rho);
   helper_.SetPackedCandidates(*packedCands);

   // determine primary vertex
   bool valid = false;
   for (const auto &v : *input_vtx) {
      if (!v.isFake() && v.ndof() >= 4 && abs(v.z()) <= 24. && abs(v.position().Rho()) <= 2.) {
         helper_.SetVertex(v);
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

   //auto input_ele = helper_.GetElectronsWithMVAid(input_ele_raw, mvaValuesGP, mvaCategoriesGP);
   auto input_ele = helper_.GetElectronsWithMVAid(input_ele_raw, mvaValuesHZZ, mvaCategoriesHZZ);

   jets_ = helper_.GetSelectedJets(*input_jet, -666., 666., jetID::none, '-'); // already corrected (?)

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

      mu.addUserFloat("isMediumMuon", isMediumMuon(mu, hip_safe_));
    
      std::map<std::string, double> miniIso_calculation_params;

      mu.addUserFloat("relIso", helper_.GetMuonRelIso(mu, coneSize::R03, corrType::rhoEA));
      mu.addUserFloat("miniIso", helper_.GetMuonRelIso(mu, coneSize::miniIso, corrType::rhoEA, &miniIso_calculation_params));
      mu.addUserFloat("miniAbsIsoCharged", miniIso_calculation_params["miniAbsIsoCharged"]);
      mu.addUserFloat("miniAbsIsoNeutral", miniIso_calculation_params["miniAbsIsoNeutral"]);
      mu.addUserFloat("rho", miniIso_calculation_params["rho"]);
      mu.addUserFloat("effArea", miniIso_calculation_params["effArea"]);
      mu.addUserFloat("miniIsoR", miniIso_calculation_params["miniIsoR"]);
      mu.addUserFloat("miniAbsIsoNeutralcorr", miniIso_calculation_params["miniAbsIsoNeutralcorr"]);
      mu.addUserFloat("localChiSq", mu.combinedQuality().chi2LocalPosition);
      mu.addUserFloat("trackKink", mu.combinedQuality().trkKink);
      // lepMVA input vars
      mu.addUserFloat("chargedRelIso", mu.pfIsolationR03().sumChargedHadronPt / mu.pt());
      mu.addUserFloat("neutralRelIso", mu.userFloat("relIso") - mu.userFloat("chargedRelIso"));

      mu.addUserFloat("sip3D", fabs(mu.dB(pat::Muon::PV3D) / mu.edB(pat::Muon::PV3D)));

      mu.addUserFloat("idLooseLJ", helper_.isGoodMuon(mu, 15., 2.4, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta) ? 1. : -666.);
      mu.addUserFloat("idTightLJ", helper_.isGoodMuon(mu, 25., 2.1, muonID::muonTight, coneSize::R04, corrType::deltaBeta) ? 1. : -666.);

      addCommonUserFloats(mu);

      mus->push_back(mu);
   }

   int ele_index_for_mva = 0;
   for (auto ele : input_ele) {
      ele_index_for_mva++;
      if (ele.pt() < ele_minpt_)
         continue;

      std::map<std::string, double> miniIso_calculation_params;

      ele.addUserFloat("superClusterEta", abs(ele.superCluster()->position().eta()));
      ele.addUserFloat("relIso", helper_.GetElectronRelIso(ele, coneSize::R03, corrType::rhoEA));
      ele.addUserFloat("miniIso",
                       helper_.GetElectronRelIso(ele, coneSize::miniIso, corrType::rhoEA, effAreaType::spring15, &miniIso_calculation_params));
      ele.addUserFloat("miniAbsIsoCharged", miniIso_calculation_params["miniAbsIsoCharged"]);
      ele.addUserFloat("miniAbsIsoNeutral", miniIso_calculation_params["miniAbsIsoNeutral"]);
      ele.addUserFloat("rho", miniIso_calculation_params["rho"]);
      ele.addUserFloat("effArea", miniIso_calculation_params["effArea"]);
      ele.addUserFloat("miniIsoR", miniIso_calculation_params["miniIsoR"]);
      ele.addUserFloat("miniAbsIsoNeutralcorr", miniIso_calculation_params["miniAbsIsoNeutralcorr"]);

      ele.addUserFloat("dxy", ele.gsfTrack()->dxy(vertex_.position()));
      ele.addUserFloat("dz", ele.gsfTrack()->dz(vertex_.position()));

      ele.addUserFloat("numMissingHits", ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
      // leptonMVA vars
      ele.addUserFloat("chargedRelIso", ele.pfIsolationVariables().sumChargedHadronPt / ele.pt());
      ele.addUserFloat("neutralRelIso", ele.userFloat("relIso") - ele.userFloat("chargedRelIso"));

      ele.addUserFloat("sip3D", fabs(ele.dB(pat::Electron::PV3D) / ele.edB(pat::Electron::PV3D)));
      ele.addUserFloat("eleMvaId", ele_mvaValuesHZZ.get(ele_index_for_mva - 1));
      ele.addUserFloat("eleMvaIdHZZ", ele_mvaValuesHZZ.get(ele_index_for_mva - 1));

      ele.addUserFloat("idLooseLJ", helper_.isGoodElectron(ele, 15., 2.4, electronID::electronEndOf15MVA80iso0p15) ? 1. : -666.);
      ele.addUserFloat("idTightLJ", helper_.isGoodElectron(ele, 30., 2.1, electronID::electronEndOf15MVA80iso0p15) ? 1. : -666.);
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
