import FWCore.ParameterSet.Config as cms

ttHLeptons = cms.EDProducer('LeptonIdentifier',
        electrons = cms.InputTag('slimmedElectrons'),
        electronMinPt = cms.double(7.0),
        jets = cms.InputTag('slimmedJets'),
        muons = cms.InputTag('slimmedMuons'),
        muonMinPt = cms.double(5.0)
)
