import FWCore.ParameterSet.Config as cms

FTLDumpHits = cms.EDAnalyzer('FTLDumpHits',
                             genParticlesTag = cms.untracked.InputTag("genParticles"),
                             simHitsTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
                             recHitsTag = cms.untracked.InputTag("mtdRecHits:FTLBarrel"),
                             tracksTag = cms.untracked.InputTag("generalTracks"),
                             crysLayout = cms.untracked.int32(0),
                             track_hit_DRMax = cms.double(0.05),
                             track_hit_distMax = cms.double(99999.),
                             treeName = cms.untracked.string("DumpHits"),
                             verbosity = cms.bool(True)
                             )
