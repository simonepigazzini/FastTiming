import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.maxEvents = -1
options.parseArguments()

process = cms.Process("FTLDumpJets")
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# import of standard configurations
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Geometry
process.load('Configuration.Geometry.GeometryExtended2023D19Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D19_cff')

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(INPUTFILELIST)
)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            
process.load('PrecisionTiming.FTLAnalysis.FTLDumpJets_cfi')
FTLDumper = process.FTLDumpJets

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:OUTPUTFILE"))

process.path = cms.Path(FTLDumper)

process.schedule = cms.Schedule(process.path)
