import FWCore.ParameterSet.Config as cms

simHitsAnlzr = cms.EDAnalyzer('muSimHits',
    simHitRPCLabel = cms.InputTag("g4SimHits:MuonRPCHits"),
    simHitCSCLabel = cms.InputTag("g4SimHits:MuonCSCHits"),
    simHitDTLabel = cms.InputTag("g4SimHits:MuonDTHits"),
    simHitGEMLabel = cms.InputTag("g4SimHits:MuonGEMHits"),
    genParticleLabel = cms.InputTag("genParticles"),
    particleId = cms.int32(13)
)
