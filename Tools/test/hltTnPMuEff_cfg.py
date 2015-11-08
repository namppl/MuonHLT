import FWCore.ParameterSet.Config as cms

process = cms.Process("EFFICIENCY")

from MuonHLT.Tools.paths25ns_cff import *

# Configuration parameters ##########

fileNames = cms.untracked.vstring(
''
)

secondaryFileNames = cms.untracked.vstring()

#####################################

process.source = cms.Source("PoolSource",
                            fileNames = fileNames,
                            secondaryFileNames = secondaryFileNames
                            )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v2') 

# Configuration parameters ##########

pathToStudy1 = Mu45
pathToStudy2 = Mu50
tagTriggerName1 = "HLT_IsoMu20_v"
tagTriggerName2 = "HLT_IsoMu20_v"

process.tnpMu45 =cms.EDAnalyzer("MuonTriggerEfficiencyAnalyzer",
                           vertexes = cms.InputTag("offlinePrimaryVertices"),
                           muons = cms.InputTag("muons"),
                           triggerProcess = cms.string("HLT"), #"TEST"
                           tagTriggerName = cms.string(tagTriggerName1),
                           triggerName = cms.string("HLT_"+pathToStudy1.TRIGNAME+"_v"),
                           probeFilterL1 = cms.string(pathToStudy1.PROBEL1),
                           probeFilterL2 = cms.string(pathToStudy1.PROBEL2),
                           probeFilterL3 = cms.string(pathToStudy1.PROBEL3),
                           probeFilterL3Iso = cms.string(pathToStudy1.PROBEL3ISO),
                           useRerun = cms.bool(False),
                           tagID = cms.string("Tight"),
                           probeID = cms.string("Tight"),
                           useIso = cms.bool(False),
                           tagPtCut = cms.double(21),
                           probePtCut = cms.double(46),
                           tagEtaCut = cms.double(2.4),
                           probeEtaCut = cms.double(2.1),
                           isolationType = cms.string("TrkIso"),
                           isolationCut = cms.double(0.1),
                           maxNumberMuons = cms.untracked.uint32(999999)
                           )

process.tnpMu50 =cms.EDAnalyzer("MuonTriggerEfficiencyAnalyzer",
                           vertexes = cms.InputTag("offlinePrimaryVertices"),
                           muons = cms.InputTag("muons"),
                           triggerProcess = cms.string("HLT"), #"TEST"
                           tagTriggerName = cms.string(tagTriggerName2),
                           triggerName = cms.string("HLT_"+pathToStudy2.TRIGNAME+"_v"),
                           probeFilterL1 = cms.string(pathToStudy2.PROBEL1),
                           probeFilterL2 = cms.string(pathToStudy2.PROBEL2),
                           probeFilterL3 = cms.string(pathToStudy2.PROBEL3),
                           probeFilterL3Iso = cms.string(pathToStudy2.PROBEL3ISO),
                           useRerun = cms.bool(False),
                           tagID = cms.string("Tight"),
                           probeID = cms.string("Tight"),
                           useIso = cms.bool(False),
                           tagPtCut = cms.double(21),
                           probePtCut = cms.double(51),
                           tagEtaCut = cms.double(2.4),
                           probeEtaCut = cms.double(2.4),
                           isolationType = cms.string("TrkIso"),
                           isolationCut = cms.double(0.1),
                           maxNumberMuons = cms.untracked.uint32(999999)
                           )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("tnp_testpt.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.outPath = cms.EndPath(process.tnpMu45 + process.tnpMu50)
