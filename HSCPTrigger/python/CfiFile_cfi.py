import FWCore.ParameterSet.Config as cms

demo2 = cms.EDAnalyzer('HSCPTrigger',
                     digisLabel = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"), #edm::DetSetVector<RPCDigiSimLink>  "simMuonRPCDigis" "RPCDigiSimLink"  "DIGI2RAW"
                     digis2Label = cms.InputTag("simMuonRPCDigis",""),
                     recHitLabel = cms.InputTag("rpcRecHits",""),
                     genEventInfo = cms.InputTag("generator"),
                     hepEventInfo = cms.InputTag("generator"),
                     triggerLabel = cms.InputTag("TriggerResults","","HLT"),
                     genParticlesLabel = cms.InputTag("genParticles","","SIM"),
                     triggNames = cms.vstring(
                     			'HLT_L1SingleMuOpen_v3'
                     			)
)
