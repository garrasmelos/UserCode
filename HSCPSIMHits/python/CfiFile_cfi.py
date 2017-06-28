import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HSCPSIMHits',
                      hitsLabel = cms.InputTag("g4SimHits","MuonRPCHits"),
                      genParticleLabel = cms.InputTag("genParticles")
                      )
