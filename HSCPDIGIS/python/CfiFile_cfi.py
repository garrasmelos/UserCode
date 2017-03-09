import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HSCPDIGIS',
      digisLabel = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"), #edm::DetSetVector<RPCDigiSimLink>     "simMuonRPCDigis"           "RPCDigiSimLink"   "DIGI2RAW"
		digis2Label = cms.InputTag("simMuonRPCDigis",""),
		recHitLabel = cms.InputTag("rpcRecHits","")
)
