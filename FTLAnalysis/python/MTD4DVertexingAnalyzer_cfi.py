import FWCore.ParameterSet.Config as cms

MTD4DVertexingAnalyzer = cms.EDAnalyzer('MTD4DVertexingAnalyzer',
                                        genParticlesTag = cms.untracked.InputTag("genParticles"),
                                        trackingParticlesTag = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                                        trackAndTrackingParticlesAssociatorMapTag = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
                                        generalTracksTag = cms.untracked.InputTag("generalTracks"),
                                        extendedTracksTag = cms.untracked.InputTag("trackExtenderWithMTD"),
                                        extTracksPathLengthTag = cms.untracked.InputTag("trackExtenderWithMTD", "pathLength", "RECO"),
                                        extTracksMTDtimeTag = cms.untracked.InputTag("trackExtenderWithMTD", "tmtd", "RECO"),
                                        t0TOFPIDTag = cms.untracked.InputTag("tofPID", "t0", "RECO"),
                                        sigmat0TOFPIDTag = cms.untracked.InputTag("tofPID", "sigmat0", "RECO"),
                                        probPiTOFPIDTag = cms.untracked.InputTag("tofPID", "probPi", "RECO"),
                                        probPTOFPIDTag = cms.untracked.InputTag("tofPID", "probP", "RECO"),
                                        probKTOFPIDTag = cms.untracked.InputTag("tofPID", "probK", "RECO"),
                                        genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
                                        genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
                                        simVtxTag = cms.untracked.InputTag("mix", "InitialVertices", "HLT"),
                                        vtx3DTag = cms.untracked.InputTag("offlinePrimaryVertices3D", "", "RECO"),
                                        vtx4DTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),                                        
                                        vtx4DNoPIDTag = cms.untracked.InputTag("offlinePrimaryVertices4DnoPID", "", "RECO"),                                        
                                        vtxsTreeName = cms.untracked.string("vtxs_tree"),
                                        trksTreeName = cms.untracked.string("trks_tree")
)
