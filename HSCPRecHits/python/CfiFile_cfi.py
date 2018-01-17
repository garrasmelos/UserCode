import FWCore.ParameterSet.Config as cms

demo2 = cms.EDAnalyzer('HSCPRecHits',
							recHitLabel = cms.InputTag("rpcRecHits"),
              simHitLabel = cms.InputTag("g4SimHits:MuonRPCHits"),
              genParticlesLabel = cms.InputTag("genParticles","","HLT"),
              triggerLabel = cms.InputTag("TriggerResults","","HLT"),
              triggNames = cms.vstring(
               'HLT_L1SingleMuOpen_v3'
              ),
              particleId = cms.int32(1000015)
              #particleId = cms.int32(13)

)	
