import FWCore.ParameterSet.Config as cms

process = cms.Process("TestPhotons")

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
   #
   # Just a handful of files from the dataset are listed below, for testing
   #
       # '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/101611CC-026E-E411-B8D7-00266CFFBF88.root',
       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1024D6DB-7D6F-E411-AE1D-00266CFF0608.root',
       # '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/107B7861-7C6F-E411-974E-00266CFFC80C.root',
       # '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12FB6345-C96C-E411-85C9-00266CFFC4D4.root',
   # Run from a local file
#      'file:/afs/cern.ch/user/i/ikrav/workspace/releases-git/CMSSW_7_2_0/src/ElectronWork/ElectronNtupler/test/00A074A5-BF72-E411-B455-003048F02CBE_mini.root'
 )
)

#
# Run some stuff to produce value maps needed for IDs
#
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

#
# Configure an example module for user analysis of electrons
#

process.ntupler = cms.EDAnalyzer('TreeWriter',
                                 photons = cms.InputTag("slimmedPhotons"),
                                 jets = cms.InputTag("slimmedJets"),
                                 mets = cms.InputTag("slimmedMETs"),
                                 rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 prunedGenParticles = cms.InputTag("prunedGenParticles"),
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
                                 phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('photon_ntuple_mva_mini.root')
                                   )

process.p = cms.Path( process.photonIDValueMapProducer * process.ntupler)
