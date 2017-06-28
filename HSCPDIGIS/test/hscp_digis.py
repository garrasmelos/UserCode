import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
          'file:/afs/cern.ch/user/g/garamire/work/private/HSCP_stau-m1599_NoPU/HSCPppstau_M_1599_NoPU.root',
    )
)
process.load('HSCPAnalysis.HSCPDIGIS.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   fileName = cms.string("hTest.root")
                   )
#process.demo = cms.EDAnalyzer('HSCPDIGIS'
#)


process.p = cms.Path(process.demo)
