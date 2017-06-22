import FWCore.ParameterSet.Config as cms

process = cms.Process("SimHitsAnlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/CMSSW_9_1_1_patch1/src/SingleMuPt10_pythia8_cfi_py_GEN_SIM_DIGI.root'
    )
)
process.load('HSCPAnalysis.muSimHits.muSimHits_cfi')

#process.simHitsAnlzr = cms.EDAnalyzer('muSimHits'
#)
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("mu_PGun_simHits.root")
)

process.p = cms.Path(process.simHitsAnlzr)
