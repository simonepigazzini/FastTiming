import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.register('crysLayout',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "crystal layout (tile, barphi, barz)")
options.maxEvents = -1
options.parseArguments()

process = cms.Process("FTLDumpHits")
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
if 'tile' in options.crysLayout:
    process.load('Configuration.Geometry.GeometryExtended2023D24Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D24_cff')
if 'barphi' in options.crysLayout:
    process.load('Configuration.Geometry.GeometryExtended2023D25Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D25_cff')
if 'barz' in options.crysLayout:
    process.load('Configuration.Geometry.GeometryExtended2023D33Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D33_cff')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')

process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")
process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")
process.mtdGeometry.applyAlignment = cms.bool(False)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(INPUTFILELIST)
    )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            
process.load('PrecisionTiming.FTLAnalysis.FTLDumpHits_cfi')
FTLDumper = process.FTLDumpHits

if 'tile' in options.crysLayout:
    FTLDumper.crysLayout = cms.untracked.int32(1)
if 'barphi' in options.crysLayout:
    FTLDumper.crysLayout = cms.untracked.int32(2)
if 'barz' in options.crysLayout:
    FTLDumper.crysLayout = cms.untracked.int32(3)

# Output TFile
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("OUTPUTFILE")
    )

process.path = cms.Path(FTLDumper)

process.schedule = cms.Schedule(process.path)
