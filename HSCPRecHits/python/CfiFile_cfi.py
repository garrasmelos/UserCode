import FWCore.ParameterSet.Config as cms

demo2 = cms.EDAnalyzer('HSCPRecHits',
							recHitLabel = cms.InputTag("rpcRecHits"),
              simHitLabel = cms.InputTag("g4SimHits:MuonRPCHits"),
              genParticlesLabel = cms.InputTag("genParticles")
							
)	
