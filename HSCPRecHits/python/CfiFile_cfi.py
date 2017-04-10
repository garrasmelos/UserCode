import FWCore.ParameterSet.Config as cms

demo2 = cms.EDAnalyzer('HSCPRecHits',
							recHitLabel = cms.InputTag("rpcRecHits",""),
							genParticlesLabel = cms.InputTag("genParticles","","DIGI")
							
)	
