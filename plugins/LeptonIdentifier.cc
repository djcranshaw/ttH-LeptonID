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
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

#include "TMVA/Reader.h"

//
// class declaration
//

enum ID { nonIsolated, preselection, fakeable, cutbased, mvabased, looseCut, looseMVA, tightCut, tightMVA };

class LeptonIdentifier : public edm::EDProducer
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

   bool qualityTrack(const reco::Track& t, const reco::Vertex& v) const;
   template<typename T> void addCommonUserFloats(T& lepton, bool useMINIAODjecs);

   // ----------member data ---------------------------
   MiniAODHelper helper_;

   TMVA::Reader *mu_reader_;
   TMVA::Reader *ele_reader_;

   Float_t varpt;
   Float_t vareta;
   Float_t varneuRelIso;
   Float_t varchRelIso;
   Float_t varjetDR_in;
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
   edm::EDGetTokenT<edm::ValueMap<float>> mva_trig_val_token_;
   edm::EDGetTokenT<edm::ValueMap<int>> mva_trig_cat_token_;
   edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesMapToken_;
   edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesMapToken_;

   reco::Vertex vertex_;
   pat::JetCollection jets_;

   double mu_minpt_;
   double ele_minpt_;
   double tau_minpt_;
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
        tau_minpt_(config.getParameter<double>("tauMinPt"))
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

   mva_trig_val_token_ =
      consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"));
   mva_trig_cat_token_ =
      consumes<edm::ValueMap<int>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"));

   mvaValuesMapToken_ =
      consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"));
   mvaCategoriesMapToken_ =
      consumes<edm::ValueMap<int>>(edm::InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"));

   // Who gives a FUCK about these parameters?  They are not used in the
   // methods we access, which could be spun off, anyways.
   helper_.SetUp("2015_74x", 666, analysisType::LJ, false);

   mu_reader_ = new TMVA::Reader("!Color:!Silent");
   ele_reader_ = new TMVA::Reader("!Color:!Silent");

   std::vector<TMVA::Reader *> ele_mvas = { ele_reader_ };

   for (auto &m : ele_mvas) {
      m->AddVariable("LepGood_pt", &varpt);
      m->AddVariable("LepGood_eta", &vareta);
      m->AddVariable("LepGood_jetNDauChargedMVASel", &varjetNDauCharged_in);
      m->AddVariable("LepGood_miniRelIsoCharged", &varchRelIso);
      m->AddVariable("LepGood_miniRelIsoNeutral", &varneuRelIso);
      m->AddVariable("LepGood_jetPtRelv2", &varjetDR_in);
      m->AddVariable("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", &varjetPtRatio_in);
      m->AddVariable("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", &varjetBTagCSV_in);
      m->AddVariable("LepGood_sip3d", &varsip3d);
      m->AddVariable("LepGood_dxy := log(abs(LepGood_dxy))", &vardxy);
      m->AddVariable("LepGood_dz := log(abs(LepGood_dz))", &vardz);
      m->AddVariable("LepGood_mvaIdSpring15", &varmvaId);
   }

   std::vector<TMVA::Reader *> mu_mvas = { mu_reader_ };

   for (auto &m : mu_mvas) {
      m->AddVariable("LepGood_pt", &varpt);
      m->AddVariable("LepGood_eta", &vareta);
      m->AddVariable("LepGood_jetNDauChargedMVASel", &varjetNDauCharged_in);
      m->AddVariable("LepGood_miniRelIsoCharged", &varchRelIso);
      m->AddVariable("LepGood_miniRelIsoNeutral", &varneuRelIso);
      m->AddVariable("LepGood_jetPtRelv2", &varjetDR_in);
      m->AddVariable("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", &varjetPtRatio_in);
      m->AddVariable("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", &varjetBTagCSV_in);
      m->AddVariable("LepGood_sip3d", &varsip3d);
      m->AddVariable("LepGood_dxy := log(abs(LepGood_dxy))", &vardxy);
      m->AddVariable("LepGood_dz := log(abs(LepGood_dz))", &vardz);
      m->AddVariable("LepGood_segmentCompatibility", &varSegCompat);
   }

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

bool
LeptonIdentifier::qualityTrack(const reco::Track& t, const reco::Vertex& v) const
{
   return
      (t.pt()>1 &&
       t.hitPattern().numberOfValidHits() >= 8 &&
       t.hitPattern().numberOfValidPixelHits() >= 2 &&
       t.normalizedChi2() < 5 &&
       std::fabs(t.dxy(v.position())) < 0.2 &&
       std::fabs(t.dz(v.position())) < 17);
}

float
LeptonIdentifier::mva(const pat::Muon &mu)
{
   varpt = mu.pt();
   vareta = mu.eta();
   varchRelIso = mu.userFloat("miniAbsIsoCharged");
   varneuRelIso = mu.userFloat("miniAbsIsoNeutral");
   varjetDR_in = mu.userFloat("nearestJetPtRatio");
   varjetPtRatio_in = mu.userFloat("nearestJetPtRatio");
   varjetBTagCSV_in = mu.userFloat("nearestJetCsv");
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
   varchRelIso = ele.userFloat("miniAbsIsoCharged");
   varneuRelIso = ele.userFloat("miniAbsIsoNeutral");
   varjetDR_in = ele.userFloat("nearestJetPtRatio");
   varjetPtRatio_in = ele.userFloat("nearestJetPtRatio");
   varjetBTagCSV_in = ele.userFloat("nearestJetCsv");
   varjetNDauCharged_in = ele.userFloat("nearestJetNDauCharged");
   varsip3d = ele.userFloat("sip3D");
   vardxy = log(fabs(ele.userFloat("dxy")));
   vardz = log(fabs(ele.userFloat("dz")));
   varmvaId = ele.userFloat("eleMvaId");

   return ele_reader_->EvaluateMVA("BDTG method");
}

// ------------ id functions ------------
bool
LeptonIdentifier::passes(const pat::Muon &mu, ID id)
{
   double minMuonPt = 5.0; // iMinPt;

   bool passesKinematics = false;
   bool passesIso = false;
   bool passesID = false;

   bool passesMuonBestTrackID = false;
   bool mediumID = false;
   bool goodGlb = false;
   bool passesCuts = false;

   switch (id) {
      case looseMVA:
         passesKinematics = true;
         passesIso = true;
         goodGlb = (mu.isGlobalMuon() && mu.userFloat("normalizedChiSq") < 3 && mu.userFloat("localChiSq") < 12 && mu.userFloat("trackKink") < 20);
         mediumID = (mu.userFloat("validFraction") >= 0.8 && mu.segmentCompatibility() >= (goodGlb ? 0.303 : 0.451));
         passesID = (mu.userFloat("leptonMVA") > 0.5 && mediumID);
         break;
      case tightMVA:
         passesKinematics = true;
         passesIso = true;
         goodGlb = (mu.isGlobalMuon() && mu.userFloat("normalizedChiSq") < 3 && mu.userFloat("localChiSq") < 12 && mu.userFloat("trackKink") < 20);
         mediumID = (mu.userFloat("validFraction") >= 0.8 && mu.segmentCompatibility() >= (goodGlb ? 0.303 : 0.451));
         passesID = (mu.userFloat("leptonMVA") > 0.75 && mediumID); /////// <<--- the MVA cut !!!
         break;
      case looseCut:
         passesKinematics = ((mu.pt() >= minMuonPt) && (fabs(mu.eta()) < 2.4));
         passesIso = (mu.userFloat("relIso") < 0.5);
         if (mu.innerTrack().isAvailable()) {
            passesMuonBestTrackID = (fabs(mu.userFloat("dxy")) < 0.05 && fabs(mu.userFloat("dz")) < 0.1);
         }
         passesID = (passesMuonBestTrackID && (mu.isGlobalMuon() || mu.isTrackerMuon()) && mu.isPFMuon());
         break;
      case tightCut:
         passesKinematics = ((mu.pt() >= minMuonPt) && (fabs(mu.eta()) < 2.4));
         passesIso = (mu.userFloat("relIso") < 0.1);

         if (mu.innerTrack().isAvailable() && mu.globalTrack().isAvailable()) {
            passesMuonBestTrackID = (fabs(mu.userFloat("dxy")) < 0.05 && fabs(mu.userFloat("dz")) < 0.1);
            passesCuts = (mu.isGlobalMuon() && mu.isPFMuon() && mu.userFloat("normalizedChiSq") < 10. && mu.userFloat("numValidMuonHits") > 0 &&
                          mu.numberOfMatchedStations() > 1 && mu.userFloat("numValidPixelHits") > 0 &&
                          mu.userFloat("trackerLayersWithMeasurement") > 5 && mu.userFloat("sip3D") < 4.);
         }
         passesID = (passesMuonBestTrackID && passesCuts);
         break;
      case preselection:
         passesKinematics = ((mu.pt() > minMuonPt) && (fabs(mu.eta()) < 2.4));
         passesIso = (mu.userFloat("miniIso") < 0.4);
         if (mu.innerTrack().isAvailable()) { // innerTrack() // muonBestTrack // isAvailable
            passesMuonBestTrackID = (fabs(mu.userFloat("dxy")) < 0.05 && fabs(mu.userFloat("dz")) < 0.1 && (mu.userFloat("sip3D") < 8));
         }
         // passesID = (( mu.isGlobalMuon() || mu.isTrackerMuon() ) && mu.isPFMuon();
         passesID = mu.isLooseMuon() && passesMuonBestTrackID;
         break;
      case fakeable:
         passesKinematics = mu.pt() > 10;
         if (mva(mu) > 0.75)
            passesID = mu.userFloat("nearestJetCsv") < 0.89;
         else
            passesID = mu.userFloat("nearestJetCsv") < 0.605 && mu.userFloat("nearestJetPtRatio") > 0.3;
         passesIso = true;
         break;
      case cutbased:
         passesKinematics = mu.pt() > 10;
         passesIso = mu.userFloat("miniIso") < 0.2;
         if (mu.muonBestTrack().isAvailable())  // muonBestTrack? innerTrack?
            passesMuonBestTrackID = mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt() < 0.2;
         passesID = passesMuonBestTrackID && mu.userFloat("sip3D") < 4 && mu.isMediumMuon();
         break;
      case mvabased:
         passesKinematics = mu.pt() > 10;
         passesIso = true;
         if (mu.innerTrack().isAvailable())  // muonBestTrack?
            passesMuonBestTrackID = mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt() < 0.2;
         passesID = mva(mu) > 0.75 &&
                    mu.userFloat("nearestJetCsv") < 0.89 &&
                    passesMuonBestTrackID &&
                    mu.isMediumMuon();
         break;
      case nonIsolated:
         edm::LogError("LeptonID") << "Invalid ID 'nonIsolated' for muons!";
         return false;
   }

   return (passesKinematics && passesIso && passesID);
}

bool
LeptonIdentifier::passes(const pat::Electron &ele, ID id)
{
   double minElectronPt = 7.; // iMinPt;

   // Be skeptical about this electron making it through
   bool passesKinematics = false;
   bool passesIso = false;
   bool passesID = false;

   bool passGsfTrackID = false;
   bool passesCuts = false;

   bool passesMVA = false;

   double eleMvaNonTrig = ele.userFloat("eleMvaId");
   float scEta = ele.userFloat("superClusterEta");

   switch (id) {
      case looseMVA:
         passesKinematics = true;
         passesIso = true;
         passesID = (ele.userFloat("leptonMVA") > 0.5 && ele.userFloat("numMissingHits") == 0 && ele.passConversionVeto());
         break;
      case tightMVA:
         passesKinematics = true;
         passesIso = true;
         passesID = (ele.userFloat("leptonMVA") > 0.75 && ele.userFloat("numMissingHits") == 0 && ele.passConversionVeto()); /////// <<--- the MVA cut !!!
         break;
      case looseCut:
         passesKinematics = ((ele.pt() >= minElectronPt) && (fabs(ele.eta()) < 2.5));
         passGsfTrackID = (fabs(ele.userFloat("dxy")) < 0.05 && fabs(ele.userFloat("dz")) < 0.1 && ele.userFloat("numMissingHits") <= 1);
         passesIso = (ele.userFloat("relIso") < 0.5);
         if (scEta <= 1.479) {
            passesCuts = (fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.007 && fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.8 &&
                          ele.full5x5_sigmaIetaIeta() < 0.01 && ele.hadronicOverEm() < 0.15);
         } else if (scEta > 1.479 && scEta < 2.5) {
            passesCuts = (fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.01 && fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.7 &&
                          ele.full5x5_sigmaIetaIeta() < 0.03);
         }
         passesID = (passesCuts && passGsfTrackID);
         break;
      case tightCut:
         passesKinematics = ((ele.pt() >= minElectronPt) && (fabs(ele.eta()) < 2.5));
         passGsfTrackID = (fabs(ele.userFloat("dxy")) < 0.05 && fabs(ele.userFloat("dz")) < 0.1 && ele.isGsfCtfScPixChargeConsistent() &&
                           ele.userFloat("numMissingHits") == 0 && ele.userFloat("sip3D") < 4);
         passesIso = (ele.userFloat("relIso") < 0.1);

         if (scEta <= 1.479) {
            passesCuts = (fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.004 && fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
                          ele.full5x5_sigmaIetaIeta() < 0.01 && ele.hadronicOverEm() < 0.12);
         } else if (scEta > 1.479 && scEta < 2.5) {
            passesCuts = (fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.007 && fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
                          ele.full5x5_sigmaIetaIeta() < 0.03 && ele.hadronicOverEm() < 0.10);
         }
         passesID = (passesCuts && passGsfTrackID);
         break;
      case preselection:

         // very loose WP
         if (scEta < 0.8)
            passesMVA = (eleMvaNonTrig > -0.7);
         else if (scEta < 1.479)
            passesMVA = (eleMvaNonTrig > -0.83);
         else
            passesMVA = (eleMvaNonTrig > -0.92);

         if (ele.gsfTrack().isAvailable()) {
            passGsfTrackID = (fabs(ele.userFloat("dxy")) < 0.05 && fabs(ele.userFloat("dz")) < 0.1 && ele.userFloat("numMissingHits") <= 1);
         }
         passesKinematics = ((ele.pt() > minElectronPt) && (fabs(ele.eta()) < 2.5));

         passesIso = ele.userFloat("miniIso") < 0.4;
         passesID = (passGsfTrackID && passesMVA) && (ele.userFloat("sip3D") < 8) && ele.passConversionVeto();
         break;
      case fakeable:
         passesKinematics = ele.pt() > 10;
         passesIso = true;
         if (ele.pt() > 30) {
            if (fabs(ele.eta()) < 0.8) {
               passesCuts = ele.sigmaIetaIeta() < 0.011 &&
                            ele.hcalOverEcal() < 0.10 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.04 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.010;
            }
            else if (fabs(ele.eta()) < 1.479) {
               passesCuts = ele.sigmaIetaIeta() < 0.011 &&
                            ele.hcalOverEcal() < 0.10 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.04 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.010;
            }
            else if (fabs(ele.eta()) < 2.5) {
               passesCuts = ele.sigmaIetaIeta() < 0.030 &&
                            ele.hcalOverEcal() < 0.07 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.008 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.07 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.005;
            }
            else
               passesKinematics = false;
         }
         else
            passesCuts = true;
         passGsfTrackID = ele.userFloat("numMissingHits") == 0;
         if (mva(ele) > 0.75)
            passesID = passesCuts && passGsfTrackID && ele.userFloat("nearestJetCsv") < 0.89;
         else
            passesID = passesCuts && passGsfTrackID && ele.userFloat("nearestJetCsv") < 0.605 && ele.userFloat("nearestJetPtRatio") > 0.3;
         break;
      case cutbased:
         // Tight WP
         if (scEta < 0.8)
            passesMVA = eleMvaNonTrig > 0.87;
         else if (scEta < 1.479)
            passesMVA = eleMvaNonTrig > 0.60;
         else
            passesMVA = eleMvaNonTrig > 0.17;
         passesKinematics = ele.pt() > 15;
         passesIso = ele.userFloat("miniIso") < 0.1;
         passGsfTrackID = ele.userFloat("numMissingHits") == 0 &&
                          ele.isGsfCtfScPixChargeConsistent()
                          && ele.userFloat("sip3D") < 4;
         passesID = passGsfTrackID && passesMVA && ele.passConversionVeto();
         break;
      case mvabased:
         passesKinematics = ele.pt() > 15;
         passesIso = true;
         if (ele.pt() > 30) {
            if (fabs(ele.eta()) < 0.8) {
               passesCuts = ele.sigmaIetaIeta() < 0.011 &&
                            ele.hcalOverEcal() < 0.10 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.04 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.010;
            }
            else if (fabs(ele.eta()) < 1.479) {
               passesCuts = ele.sigmaIetaIeta() < 0.011 &&
                            ele.hcalOverEcal() < 0.10 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.04 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.010;
            }
            else if (fabs(ele.eta()) < 2.5) {
               passesCuts = ele.sigmaIetaIeta() < 0.030 &&
                            ele.hcalOverEcal() < 0.07 &&
                            ele.deltaEtaSuperClusterTrackAtVtx() < 0.008 &&
                            ele.deltaPhiSuperClusterTrackAtVtx() < 0.07 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() > -0.5 &&
                            1.0/ele.ecalEnergy() - 1.0/ele.p() < 0.005;
            }
            else
               passesKinematics = false;
         }
         else
            passesCuts = true;
         passGsfTrackID = ele.userFloat("numMissingHits") == 0 && ele.isGsfCtfScPixChargeConsistent();
         passesID = passesCuts &&
                    passGsfTrackID &&
                    ele.passConversionVeto() &&
                    mva(ele) > 0.75 &&
                    ele.userFloat("nearestJetCsv") < 0.89;
         break;
      case nonIsolated:
         edm::LogError("LeptonID") << "Invalid ID 'nonIsolated' for electrons!";
         return false;
   }

   return (passesKinematics && passesIso && passesID);
}

bool
LeptonIdentifier::passes(const pat::Tau &tau, ID id)
{
   double minTauPt = 5.; // iMinPt;

   bool passesKinematics = false;
   bool passesIso = false;
   bool passesID = false;

   bool passesPVassoc = false;

   switch (id) {
      case nonIsolated:
         passesKinematics = ((tau.pt() > minTauPt) && (fabs(tau.eta()) < 2.3));
         passesIso = true;
         passesPVassoc = (fabs(tau.userFloat("dxy")) < 1000) && (fabs(tau.userFloat("dz")) < 0.2);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      case preselection:
         passesKinematics = ((tau.pt() > minTauPt) && (fabs(tau.eta()) < 2.3));
         passesIso = (tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5);
         passesPVassoc = (fabs(tau.userFloat("dxy")) < 1000) && (fabs(tau.userFloat("dz")) < 0.2);
         passesID = (tau.tauID("decayModeFinding") > 0.5) && passesPVassoc;
         break;
      case looseCut:
      case looseMVA:
      case tightCut:
      case tightMVA:
      default:
         break;
   }

   return (passesKinematics && passesIso && passesID);
}

template<typename T> void
LeptonIdentifier::addCommonUserFloats(T& lepton, bool useMINIAODjecs)
{
   double L2L3_SF = 1.;
   pat::Jet matchedJet;
   pat::Jet matchedJetL1;
   double dR = 666.;

   for (const auto &j : jets_) {
      double newDR = helper_.DeltaR(&j, &lepton);
      if (newDR < dR) {
         dR = newDR;
         matchedJet = j;
         matchedJetL1 = j;

         if (useMINIAODjecs)
            matchedJetL1.setP4(j.correctedJet(1).p4());
         L2L3_SF = matchedJet.p4().E() / matchedJetL1.p4().E();
      }
   }

   lepton.addUserFloat("nearestJetDr", min(dR, 0.5)); // no longer used in MVA

   lepton.addUserFloat("nearestJetCsv", 0.);
   lepton.addUserFloat("nearestJetPtRatio", 1.);
   lepton.addUserFloat("nearestJetPtRel", 0.);
   lepton.addUserFloat("nearestJetNDauCharged", 0.);

   // cout << jets_.size() << endl;
   if (jets_.size() > 0) {
      lepton.addUserFloat("nearestJetCsv", max(matchedJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), float(0.0)));

      auto constituents = matchedJet.daughterPtrVector();
      int n_charged_tracks = std::count_if(constituents.begin(), constituents.end(), [&](const reco::CandidatePtr& p) -> bool {
            return (
                  helper_.DeltaR(&lepton, &(p->p4())) <= 0.4
                  and
                  p->charge() != 0
                  and
                  p->bestTrack()
                  and
                  qualityTrack(*(p->bestTrack()), vertex_)
            );
      });
      lepton.addUserFloat("neadestJetNDauCharged", n_charged_tracks);

      if (useMINIAODjecs and (matchedJet.correctedJet(0).p4() - lepton.p4()).Rho() >= 1e-4 && dR <= 0.5) {
         auto lepAwareJetp4 = (matchedJetL1.p4() - lepton.p4()) * L2L3_SF + lepton.p4(); // "lep-aware" JEC
         TLorentzVector muTLV = TLorentzVector(lepton.p4().Px(), lepton.p4().Py(), lepton.p4().Pz(), lepton.p4().E());
         TLorentzVector jetTLV = TLorentzVector(lepAwareJetp4.Px(), lepAwareJetp4.Py(), lepAwareJetp4.Pz(), lepAwareJetp4.E());
         lepton.addUserFloat("nearestJetPtRatio", std::min(lepton.pt() / lepAwareJetp4.pt(), 1.5));
         lepton.addUserFloat("nearestJetPtRel", muTLV.Perp((jetTLV - muTLV).Vect()));
      }
   }

   lepton.addUserFloat("leptonMVA", mva(lepton));

   lepton.addUserFloat("idPreselection", passes(lepton, preselection));
   lepton.addUserFloat("idLooseCut", passes(lepton, looseCut));
   lepton.addUserFloat("idTightCut", passes(lepton, tightCut));

   if (lepton.userFloat("idPreselection") > .5) {
      lepton.addUserFloat("leptonMVA", mva(lepton));
      lepton.addUserFloat("idLooseMVA", passes(lepton, looseMVA));
      lepton.addUserFloat("idTightMVA", passes(lepton, tightMVA));
      lepton.addUserFloat("idFakeable", passes(lepton, fakeable));
      lepton.addUserFloat("idCutBased", passes(lepton, cutbased));
      lepton.addUserFloat("idMVABased", passes(lepton, mvabased));
   } else {
      lepton.addUserFloat("leptonMVA", -666.);
      lepton.addUserFloat("idLooseMVA", -666.);
      lepton.addUserFloat("idTightMVA", -666.);
      lepton.addUserFloat("idFakeable", -666.);
      lepton.addUserFloat("idCutBased", -666.);
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

   edm::Handle<edm::ValueMap<float>> mvaValues;
   edm::Handle<edm::ValueMap<int>> mvaCategories;

   edm::Handle<edm::ValueMap<float>> mvaTrigValues;
   edm::Handle<edm::ValueMap<int>> mvaTrigCategories;

   event.getByToken(rho_token_, rho);
   event.getByToken(packedCand_token_, packedCands);
   event.getByToken(ele_token_, input_ele_raw);
   event.getByToken(jet_token_, input_jet);
   event.getByToken(mu_token_, input_mu);
   event.getByToken(tau_token_, input_tau);
   event.getByToken(vtx_token_, input_vtx);
   event.getByToken(mva_trig_val_token_, mvaTrigValues);
   event.getByToken(mva_trig_cat_token_, mvaTrigCategories);
   event.getByToken(mvaValuesMapToken_, mvaValues);
   event.getByToken(mvaCategoriesMapToken_, mvaCategories);

   const edm::ValueMap<float> ele_mvaValues = (*mvaValues.product());

   helper_.SetRho(*rho);
   helper_.SetPackedCandidates(*packedCands);

   // determine primary vertex
   for (const auto &v : *input_vtx) {
      if (!v.isFake() && v.ndof() >= 4 && abs(v.z()) <= 24. && abs(v.position().Rho()) <= 2.) {
         helper_.SetVertex(v);
         vertex_ = v;
         break;
      }
   }

   auto input_ele = helper_.GetElectronsWithMVAid(input_ele_raw, mvaTrigValues, mvaTrigCategories);

   // auto raw_jets = helper_.GetUncorrectedJets(*input_jet);

   //    cout << "first raw jet energy before corr: " << (*input_jet)[0].correctedJet(0).p4().E() << endl;
   // const JetCorrector* corrector = JetCorrector::getJetCorrector("ak4PFchsL1L2L3", setup);
   // const JetCorrector* corrector = JetCorrector::getJetCorrector("ak4PFCHSL1L2L3Residual", setup);
   // helper_.SetJetCorrector(corrector);
   // auto raw_jets = helper_.GetUncorrectedJets(*input_jet);
   // auto corr_jets = helper_.GetCorrectedJets(raw_jets, event, setup);
   // jets_ = helper_.GetSelectedJets(corr_jets, 5., 2.4, jetID::none, '-');

   bool useMINIAODjecs = true;
   // if(!input_jet.jecSetsAvailable()) useMINIAODjecs = false;

   jets_ = helper_.GetSelectedJets(*input_jet, 5., 2.4, jetID::none, '-'); // already corrected (?)

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

      std::map<std::string, double> miniIso_calculation_params;

      mu.addUserFloat("relIso", helper_.GetMuonRelIso(mu, coneSize::R03, corrType::rhoEA));
      mu.addUserFloat("miniIso", helper_.GetMuonRelIso(mu, coneSize::miniIso, corrType::rhoEA, miniIso_calculation_params));
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

      addCommonUserFloats(mu, useMINIAODjecs);

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
                       helper_.GetElectronRelIso(ele, coneSize::miniIso, corrType::rhoEA, effAreaType::spring15, miniIso_calculation_params));
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
      ele.addUserFloat("eleMvaId", ele_mvaValues.get(ele_index_for_mva - 1));

      ele.addUserFloat("idLooseLJ", helper_.isGoodElectron(ele, 15., 2.4, electronID::electronEndOf15MVA80iso0p15) ? 1. : -666.);
      ele.addUserFloat("idTightLJ", helper_.isGoodElectron(ele, 30., 2.1, electronID::electronEndOf15MVA80iso0p15) ? 1. : -666.);

      addCommonUserFloats(ele, useMINIAODjecs);

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
      tau.addUserFloat("dxy", -666.);
      tau.addUserFloat("dz", -666.);
      tau.addUserFloat("idNonIsolated", -666.);
      tau.addUserFloat("idPreselection", -666.);

      if (tau.leadChargedHadrCand().isAvailable()) {
         auto track = tau.leadChargedHadrCand()->bestTrack();

         if (track) {
            tau.addUserFloat("dxy", track->dxy(vertex_.position()));
            tau.addUserFloat("dz", track->dz(vertex_.position()));
            tau.addUserFloat("idNonIsolated", passes(tau, nonIsolated));
            tau.addUserFloat("idPreselection", passes(tau, preselection));
         }
      }

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
