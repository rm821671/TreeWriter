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
    'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/101611CC-026E-E411-B8D7-00266CFFBF88.root',
    # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
    # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1024D6DB-7D6F-E411-AE1D-00266CFF0608.root',
    # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/107B7861-7C6F-E411-974E-00266CFFC80C.root',
    # 'root://xrootd-cms.infn.it///store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12FB6345-C96C-E411-85C9-00266CFFC4D4.root',
    # muons
    # 'root://xrootd-cms.infn.it///store/relval/CMSSW_7_2_0/SingleMu/MINIAOD/PRE_R_72_V10A_RelVal_mu2012D-v2/00000/12AF52D0-945A-E411-A94D-0025905A48C0.root',
    # electrons
    # 'root://xrootd-cms.infn.it///store/relval/CMSSW_7_2_0/SingleElectron/MINIAOD/PRE_R_72_V10A_RelVal_electron2012D-v2/00000/10721F98-F459-E411-B363-0025905B8610.root'
    # Run from a local file
    # 'file:/afs/cern.ch/user/j/jolange/private/data/minoAOD/singleElectron/964E27EF-9B5A-E411-8E30-0026189438FD.root'
))


process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(200.),
                                    photon_pT_cut=cms.untracked.double(90.),
                                    # physics objects
                                    photons = cms.InputTag("slimmedPhotons"),
                                    jets = cms.InputTag("slimmedJets"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    beamSpot = cms.InputTag('offlineBeamSpot'),
                                    conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
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
                                    # photon IDs
                                    photonLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                    photonMediumIdMap= cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                                    photonTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
                                    # electron IDs
                                    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                                    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
                                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
                                    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('photonTree.root')

)



######################
# PHOTONS, ELECTRONS #
######################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD

# turn on VID producer, indicate data format to be DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

## TODO: should not be necessary any longer, once the 74X MVA-ID recipe is available
# Run some stuff to produce value maps needed for IDs
process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")

####################
#     RUN          #
####################

process.p = cms.Path(process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.TreeWriter)
