import FWCore.ParameterSet.Config as cms

ttHLeptons = cms.EDProducer(
    'LeptonIdentifier',
    rhoParam=cms.InputTag('fixedGridRhoFastjetAll'),
    electrons=cms.InputTag('slimmedElectrons'),
    electronMinPt=cms.double(7.0),
    jets=cms.InputTag('slimmedJets'),
    muons=cms.InputTag('slimmedMuons'),
    muonMinPt=cms.double(5.0),
    taus=cms.InputTag('slimmedTaus'),
    tauMinPt=cms.double(20.0),
    LooseCSVWP=cms.double(0.46),
    MediumCSVWP=cms.double(0.80),
    IsHIPSafe=cms.bool(False),
    JECTag=cms.string('')
)
