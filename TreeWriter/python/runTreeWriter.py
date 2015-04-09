import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
        # Just a handful of files from the dataset are listed below, for testing
        # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/101611CC-026E-E411-B8D7-00266CFFBF88.root',
        'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1024D6DB-7D6F-E411-AE1D-00266CFF0608.root',
        # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/107B7861-7C6F-E411-974E-00266CFFC80C.root',
        # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12FB6345-C96C-E411-85C9-00266CFFC4D4.root',
        # muons
        # "/store/relval/CMSSW_7_2_0/SingleMu/MINIAOD/PRE_R_72_V10A_RelVal_mu2012D-v2/00000/12AF52D0-945A-E411-A94D-0025905A48C0.root",
        # electrons
        # '/store/relval/CMSSW_7_2_0/SingleElectron/MINIAOD/PRE_R_72_V10A_RelVal_electron2012D-v2/00000/10721F98-F459-E411-B363-0025905B8610.root'
        # Run from a local file
        # 'file:/afs/cern.ch/user/j/jolange/private/data/minoAOD/singleElectron/964E27EF-9B5A-E411-8E30-0026189438FD.root'
        ))

#
# Run some stuff to produce value maps needed for IDs
#
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # physics objects
                                    photons = cms.InputTag("slimmedPhotons"),
                                    jets = cms.InputTag("slimmedJets"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    # photon id related values
                                    full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                    full5x5SigmaIEtaIPhiMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIPhi"),
                                    full5x5E1x3Map      = cms.InputTag("photonIDValueMapProducer:phoFull5x5E1x3"),
                                    full5x5E2x2Map      = cms.InputTag("photonIDValueMapProducer:phoFull5x5E2x2"),
                                    full5x5E2x5MaxMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5E2x5Max"),
                                    full5x5E5x5Map      = cms.InputTag("photonIDValueMapProducer:phoFull5x5E5x5"),
                                    esEffSigmaRRMap     = cms.InputTag("photonIDValueMapProducer:phoESEffSigmaRR"),
                                    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    # electron IDs
                                    electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-veto"),
                                    electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-loose"),
                                    electronMediumIdMap= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-medium"),
                                    electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-tight"),
                                    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('photon_ntuple_mva_mini.root')
                                   )


####################
#     ELECTRONS    #
####################

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V1_miniAOD_cff']
#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

####################
#     RUN          #
####################

process.p = cms.Path(process.egmGsfElectronIDSequence * process.photonIDValueMapProducer * process.TreeWriter)
