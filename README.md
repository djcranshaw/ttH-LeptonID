# Lepton ID for Multi-Lepton Analysis

Installation:

    cmsrel CMSSW_7_4_7
    cd CMSSW_7_4_7/src
    cmsenv

    git cms-init
    git remote add cmg-central https://github.com/CERN-PH-CMG/cmg-cmssw.git
    git fetch cmg-central
    cat <<EOF >.git/info/sparse-checkout
    /.gitignore/
    /CMGTools/TTHAnalysis/data/
    EOF
    git checkout -b CMGTools-from-CMSSW_7_4_7 cmg-central/CMGTools-from-CMSSW_7_4_7
    git clone git@github.com:cms-ttH/MiniAOD.git
    git clone git@github.com:cms-ttH/ttH-LeptonID.git ttH/LeptonID

    scram b -j 8

Then include the following configuration in your parameter set:

    from RecoJets.Configuration.RecoJets_cff import *
    from RecoJets.Configuration.RecoPFJets_cff import *
    from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
    from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
    from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

    process.ak4PFCHSL1Fastjet = cms.ESProducer(
            'L1FastjetCorrectionESProducer',
            level       = cms.string('L1FastJet'),
            algorithm   = cms.string('AK4PFchs'),
            srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' ),
            useCondDB = cms.untracked.bool(True)
            )
    process.ak4PFchsL2Relative   =  ak5PFL2Relative.clone( algorithm = 'AK4PFchs' )
    process.ak4PFchsL3Absolute   =  ak5PFL3Absolute.clone( algorithm = 'AK4PFchs' )
    process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
            correctors = cms.vstring(
                'ak4PFCHSL1Fastjet',
                'ak4PFchsL2Relative',
                'ak4PFchsL3Absolute'),
            useCondDB = cms.untracked.bool(True)
    )

    process.load("ttH.LeptonID.ttHLeptons_cfi")

and change the input tags of electrons and muons to `ttHLeptons`.  Note
that only electrons with transverse momentum of 7 GeV are kept (5 GeV for
muons.)  To access the lepton MVA values and ID variables use:

    mu.userFloat("leptonMVA")
    mu.userFloat("idPreselection") > .5
    mu.userFloat("idLooseCut") > .5
    mu.userFloat("idTightCut") > .5
    mu.userFloat("idLooseMVA") > .5
    mu.userFloat("idTightMVA") > .5

and similar for electrons.
