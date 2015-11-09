import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import re

import datetime
import calendar

options = VarParsing ('analysis')
options.register ('dataset',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name of the dataset, used to do further settings")
options.register ('user',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name the user. If not set by crab, this script will determine it.")
# get and parse the command line arguments
options.parseArguments()

# determine user if not set by crab
user={}
if not options.user:
    import getpass
    user["name"]=getpass.getuser()
else:
    user["name"]=options.user
# user settings
if user["name"]=="kiesel":
    user["HT_cut"]=500.
    user["photon_pT_cut"]=90.
    user["isolatedPhotons"]=False
elif user["name"]=="lange":
    user["HT_cut"]=0.
    user["photon_pT_cut"]=20.
    user["isolatedPhotons"]=True
else:
    print "you shall not pass!"
    print "(unkown user '%s')"%options.user
    exit()

# setup any defaults you want
#options.inputFiles = 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/027A951D-4103-E511-8B6B-A0040420FE80.root'
options.inputFiles = "file:/home/home4/institut_1b/rmeyer/cms_sw/MiniAOD/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/00759690-D16E-E511-B29E-00261894382D.root"

timestamp = str(calendar.timegm(datetime.datetime.utcnow().timetuple()))
filename = '/home/home4/institut_1b/rmeyer/cms_sw/myTrees/mcDY_' + timestamp + '.root'

options.outputFile = filename
options.maxEvents = -1

# determine if Data or Simulation
isRealData=(not options.dataset.endswith("SIM"))

# the actual TreeWriter module
process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")

# determine global tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
gtName = "auto:run2_data" if isRealData else "auto:run2_mc"

# for further global tags, see here:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
if re.match( "/.*/Run2015.*05Oct2015.*/MINIAOD", options.dataset ): gtName = "74X_dataRun2_reMiniAOD_v0,frontier://FrontierArc/CMS_CONDITIONS_EJ06"
if re.match( "/.*/Run2015D-PromptReco-v4/MINIAOD", options.dataset ): gtName = "74X_dataRun2_Prompt_v4,frontier://FrontierArc/CMS_CONDITIONS_EJ06"
if re.match( "/.*/.*RunIISpring15MiniAODv2.*/MINIAODSIM", options.dataset ): gtName = "74X_mcRun2_asymptotic_v2,frontier://FrontierArc/CMS_CONDITIONS_EJ06"
process.GlobalTag = GlobalTag(process.GlobalTag, gtName, '')

hardPUveto=True if options.dataset.startswith("/QCD_HT100to200") else False

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
    options.inputFiles
))


###############################
# Define MET Filters to apply #
###############################
# See https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=39
applyMetFilters=cms.untracked.vstring(
    "Flag_CSCTightHaloFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"
)
# HBHE has to be manually re-run for early data.
# This is not applied as EDFilter, as suggested, but manually
# checked in TreeWriter (otherwise the "initial" event count
# is wrong)
# TODO: remove, when fixed upstream
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(user["HT_cut"]),
                                    photon_pT_cut=cms.untracked.double(user["photon_pT_cut"]),
                                    isolatedPhotons=cms.untracked.bool(user["isolatedPhotons"]),
                                    # physics objects
                                    photons = cms.InputTag("slimmedPhotons"),
                                    jets = cms.InputTag("slimmedJets"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    	# packed pf candidates
                                    packedCandidates = cms.InputTag("packedPFCandidates"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    # electron IDs
                                    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                    # photon IDs
                                    # TODO: update to 25ns recommendations
                                    photonLooseIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
                                    photonMediumIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
                                    photonTightIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight"),
                                    photonMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
                                    # met filters to apply
                                    metFilterNames=applyMetFilters,
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string( "pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU" ),
                                    HBHENoiseFilterResult = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
                                    HBHEIsoNoiseFilterResult = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
                                    hardPUveto=cms.untracked.bool(hardPUveto),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(
                                        "HLT_Photon90_CaloIdL_PFHT500_v",
                                        "HLT_Photon90_v",
                                        "HLT_PFHT600_v",
                                        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
                                        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v",
                                        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
                                        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v",
                                        "HLT_Photon22_v",
                                        "HLT_Photon30_v",
                                        "HLT_Photon36_v",
                                        "HLT_Photon50_v",
                                        "HLT_Photon165_R9Id90_HE10_IsoM_v",
                                        "HLT_Photon36_R9Id90_HE10_IsoM_v",
                                        "HLT_Photon135_PFMET100_v",
                                        "HLT_Photon175_v",
                                        "HLT_Photon500_v",
                                        "HLT_PFMET170_v",
                                        "HLT_IsoMu18_v",
                                        "HLT_IsoMu20_v",
                                        "HLT_Mu20_v",
                                        "HLT_Mu50_v",
                                    )
)

if re.match( "/SMS-.*/.*/USER", options.dataset ):
    # signal scan
    process.TreeWriter.metFilterNames = [] # no met filters for fastsim
    process.TreeWriter.pileUpSummary = "addPileupInfo" # for miniaod v1

process.TFileService = cms.Service("TFileService",fileName = cms.string(options.outputFile))


######################
# PHOTONS, ELECTRONS #
######################
from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import *

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD

# turn on VID producer, indicate data format to be DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer  (process, dataFormat)

# define which IDs we want to produce
el_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
ph_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff',
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in el_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

####################
#     RUN          #
####################

process.p = cms.Path(
    process.photonIDValueMapProducer
    *process.egmGsfElectronIDSequence
    *process.egmPhotonIDSequence
    *process.HBHENoiseFilterResultProducer #produces HBHE bools (applied in TreeWriter manually)
    *process.TreeWriter
)
