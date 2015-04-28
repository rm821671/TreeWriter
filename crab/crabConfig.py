from WMCore.Configuration import Configuration
import os

cmssw_src=os.environ['CMSSW_BASE']+'/src/'

config = Configuration()

config.section_("General")
config.General.requestName   = 'CRAB3-test1'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName    = cmssw_src+'TreeWriter/TreeWriter/python/runTreeWriter.py'
config.JobType.inputFiles  = [cmssw_src+'TreeWriter/PUreweighting/puWeights.root']

config.section_("Data")
config.Data.inputDataset = '/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
# config.Data.inputDataset = '/SingleElectron/CMSSW_7_2_0-PRE_R_72_V10A_RelVal_electron2012D-v2/MINIAOD'
# config.Data.inputDataset = '/SingleMu/CMSSW_7_2_0-PRE_R_72_V10A_RelVal_mu2012D-v2/MINIAOD'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.publication = False
# This string is used to construct the output dataset name
config.Data.publishDataName = 'CRAB3-test'
config.Data.outLFN = "/store/user/jolange/data/Test/Test"

# These values only make sense for processing data
#    Select input data based on a lumi mask
# config.Data.lumiMask = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#    Select input data based on run-ranges
# config.Data.runRange = '190456-194076'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_DE_RWTH'
