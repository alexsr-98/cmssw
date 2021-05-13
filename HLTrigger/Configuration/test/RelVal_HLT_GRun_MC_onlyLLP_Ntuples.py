# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RelVal --step=HLT:GRun --conditions=auto:run2_mc_GRun --filein=file:RelVal_Raw_GRun_MC.root --custom_conditions= --fileout=RelVal_HLT_GRun_MC.root --number=100 --mc --no_exec --datatier SIM-DIGI-RAW-HLTDEBUG --eventcontent=FEVTDEBUGHLT --customise=HLTrigger/Configuration/CustomConfigs.L1THLT --era=Run2_2018 --customise= --scenario=pp --python_filename=RelVal_HLT_GRun_MC.py --processName=HLT
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('HLT2',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_GRun_onlyLLP_SingleMu_cff_v2')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("HLTrigger.Configuration.MuonMatcher_cfg")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:///eos/cms/store/user/folguera/MuonHLT/Run3/DisplacedMuons_Pt30to100_Dxy0to3000-pythia8-gun_Run3Winter20DRMiniAOD_NoPU_110X_mcRun3_2021_realistic_v6-v3_GEN-SIM-RAW.root'
#        '/store/mc/Run3Winter20DRMiniAOD/DisplacedMuons_Pt30to100_Dxy0to3000-pythia8-gun/GEN-SIM-RAW/NoPU_110X_mcRun3_2021_realistic_v6-v3/10000/3C7EE72E-B3EB-2948-8DDA-8E33A9A7CDCE.root'
#        '/store/mc/Run3Winter20DRPremixMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8_HCAL/GEN-SIM-DIGI-RAW/110X_mcRun3_2021_realistic_v6-v2/280000/01F2B624-56C1-8940-A735-14F363072EBA.root'
        '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/13625078-168A-9A4F-BD90-9D105B09C04E.root'
        #"/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW"
),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('SIM-DIGI-RAW-HLTDEBUG'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('RelVal_HLT_GRun_MC.root'),
    ###outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = cms.untracked.vstring('drop *',
    'keep *_genParticles_*_*',
    'keep recoGenParticles_*_*_*',
    'keep *_hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q_*_*',
    'keep *_hltL1fL1sMu22L1Filtered0_*_*',
    'keep *_hltL3fL1f0L2NoVtx15Filtered24Displaced_*_*',
    'keep *_hltL2fL1f0L2NoVtx24QChaCosmicSeed_*_*',
    'keep *_hltL3fL1f0L2NVf16L3NoFiltersNoVtxFiltered24Displaced_*_*',
    'keep *_hltIter3IterL3MuonL2Candidates_*_*',
    'keep *_hltIterL3MuonCandidates_*_*',
    'keep *_hltIterL3OIGlbDisplacedMuonCandidates_*_*',
    'keep *_hltIterL3OIL3MuonCandidates_*_*',
    'keep *_hltL2MuonCandidates_*_*',
    'keep *_hltL2MuonCandidatesNoVtx_*_*',
    'keep *_hltL2MuonCandidatesNoVtxMeanTimerCosmicSeed_*_*',
    'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*',
    'keep l1tMuonBXVector_*_*_*'),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition  
##SF process.muonNtuples = cms.EDAnalyzer("MuonNtuples", 
##SF                                      offlineVtx = cms.InputTag(""), 
##SF                                      offlineMuons = cms.InputTag(""), 
##SF                                      triggerResult = cms.untracked.InputTag("TriggerResults::HLT2"), 
##SF                                      triggerSummary = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT2"), 
##SF                                      tagTriggerResult = cms.untracked.InputTag("TriggerResults::HLT"), 
##SF                                      tagTriggerSummary = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
##SF                                      L3Candidates = cms.untracked.InputTag("hltIterL3OIGlbDisplacedMuonCandidates"), 
##SF                                      L2Candidates = cms.untracked.InputTag("hltL2MuonCandidatesNoVtxMeanTimerCosmicSeed"), 
##SF                                      L1Candidates = cms.untracked.InputTag('hltGtStage2Digis','Muon'), 
##SF                                      TkMuCandidates = cms.untracked.InputTag(""), 
##SF                                      L3OIMuCandidates = cms.untracked.InputTag("hltIterL3OIGlbDisplacedMuonCandidates"), 
##SF                                      L3IOMuCandidates = cms.untracked.InputTag(""),
##SF                                      theTrackOI = cms.untracked.InputTag("hltIterL3OIGlbDisplacedMuonTrackSelectionLoose"), 
##SF                                      theTrackIOL2 = cms.untracked.InputTag(""), 
##SF                                      theTrackIOL1 = cms.untracked.InputTag(""), 
##SF                                      lumiScalerTag = cms.untracked.InputTag("scalersRawToDigi"), 
##SF                                      puInfoTag = cms.untracked.InputTag("addPileupInfo"), 
##SF                                      genParticlesTag = cms.untracked.InputTag("genParticles"), 
##SF                                      doOffline = cms.untracked.bool(False) 
##SF )

##SF process.TFileService = cms.Service("TFileService", 
##SF                                    fileName = cms.string("muonNtuple.root"), 
##SF                                    closeFileFast = cms.untracked.bool(False) 
##SF ) 
##SF process.HLTValidation = cms.EndPath( process.muonNtuples )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.HLTMuonSeq = cms.Sequence(process.dispGenEta)
process.HLTMuonPath = cms.Path(process.HLTMuonSeq)
# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.HLTMuonPath])
process.schedule.extend([process.endjob_step,
                         process.FEVTDEBUGHLToutput_step,
                     #    process.HLTValidation
                     ])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.CustomConfigs
from HLTrigger.Configuration.CustomConfigs import L1THLT 

#call to customisation function L1THLT imported from HLTrigger.Configuration.CustomConfigs
process = L1THLT(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



