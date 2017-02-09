import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         'file:/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/generation/histos/CMSSW_7_4_14/src/histograms/step1.root' 
   )
)
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("hTest.root")

)
process.demo = cms.EDAnalyzer('HSCPanalyzer'
)


process.p = cms.Path(process.demo)
