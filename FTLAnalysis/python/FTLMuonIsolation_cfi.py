import FWCore.ParameterSet.Config as cms

FTLMuonIsolation = cms.EDAnalyzer(
    "FTLMuonIsolation",
    ###---Input tags
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),	    
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("particleFlow", "", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIMI"),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---Target time resolution (assumes sample were made with 30ps track t resolution)
    #targetResolutions = cms.untracked.vdouble(0.03, 0.05, 0.07, 0.09, 0.15),
    targetResolutions = cms.untracked.vdouble(0.03),
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Vtx choice option
    useMCTruthPV = cms.untracked.bool(False),
    ###---precessing type w/ or w/o timing
    isTimingSample = cms.untracked.bool(False),
    ###---save tracks info
    saveTracksInfo = cms.untracked.bool(False),
    ###---no ETL, get time in endcap from parametrized HGC response
    HGCToySim = cms.untracked.bool(False),
    ###---Trk-vtx dz cut
    dzCut = cms.untracked.VPSet(
        cms.PSet(
            absEtaMax = cms.untracked.double(0.5),
            dzCutParams = cms.untracked.vdouble(0.152865, 0.343582, -0.140844)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(1.0),
            dzCutParams = cms.untracked.vdouble(0.178083, 0.359211, -0.16321)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(1.5),
            dzCutParams = cms.untracked.vdouble(0.229793, 0.306094, -0.157445)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(1.75),
            dzCutParams = cms.untracked.vdouble(0.280826, 0.563838, -0.328811)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(2.0),
            dzCutParams = cms.untracked.vdouble(0.332952, 0.947662, -0.768322)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(2.25),
            dzCutParams = cms.untracked.vdouble(0.412667, 1.13675, -1.05305)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(2.5),
            dzCutParams = cms.untracked.vdouble(0.494174, 1.13675, -1.05305)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(2.75),
            dzCutParams = cms.untracked.vdouble(0.52859, 1.13675, -1.05305)
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(3.0),
            dzCutParams = cms.untracked.vdouble(0.707472, 1.13675, -1.05305)
        )
    ),
    ###---minimum track pt
    ptCut = cms.untracked.double(0.4),    
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3)
)

FTLMuonIsolationHGCToy = cms.EDAnalyzer(
    "FTLMuonIsolationHGC",
    ###---Input tags
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),	    
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("particleFlow", "", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---Target time resolution (assumes sample were made with 30ps track t resolution)
    targetResolutions = cms.untracked.vdouble(0.03), #ns
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Vtx choice option
    useMCTruthPV = cms.untracked.bool(False),
    ###---precessing type w/ or w/o timing
    isTimingSample = cms.untracked.bool(False),
    ###---Trk-vtx dz cut
    dzCut = cms.untracked.double(0.1),
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3, 0.4),
    ###---Time resolution eta and pt bins
    absEtaBins = cms.untracked.VPSet(
        cms.PSet(
            absEtaMax = cms.untracked.double(1.5),
            ptBins = cms.untracked.vdouble(2.0, 7000.),
            timeRes = cms.untracked.vdouble(0.07, 0.03) #ns
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(4.0),
            ptBins = cms.untracked.vdouble(1.0, 7000.),
            timeRes = cms.untracked.vdouble(-1, 0.03) #ns
        )
    )
)

