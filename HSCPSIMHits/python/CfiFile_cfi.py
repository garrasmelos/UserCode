import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HSCPSIMHits',
                      hitsLabel = cms.InputTag("g4SimHits","MuonRPCHits"),
                      genEventInfo = cms.InputTag("generator"),
                      hepEventInfo = cms.InputTag("generator")
                      )
