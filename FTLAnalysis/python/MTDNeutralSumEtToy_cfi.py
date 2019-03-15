import FWCore.ParameterSet.Config as cms

MTDNeutralSumEtToy = cms.EDAnalyzer('MTDNeutralSumEtToy',
                                     pfCandidatesTag = cms.untracked.InputTag("particleFlow"),
                                     genXYZTag = cms.untracked.InputTag("genParticles", "xyz0"),
                                     genT0Tag = cms.untracked.InputTag("genParticles", "t0"),
                                     timeResolutions = cms.untracked.vdouble(0.035, 0.04, 0.07),                             
                                     outTreeName = cms.untracked.string("sumet_tree")
)
