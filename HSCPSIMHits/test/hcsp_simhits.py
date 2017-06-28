import FWCore.ParameterSet.Config as cms
#import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Demo")

#mylist = FileUtils.loadListFromFile('list.txt')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D12_cff')
process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #     '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_1.root',
    #fileNames = cms.untracked.vstring (*mylist)
          'file:/afs/cern.ch/user/g/garamire/work/private/HSCP_stau-m1599_NoPU/HSCPppstau_M_1599_NoPU.root',
    )
)

process.load('HSCPAnalysis.HSCPSIMHits.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
             fileName = cms.string("hTest.root")
)
#process.demo = cms.EDAnalyzer('HSCPSIMHits'
#         hitsLabel = cms.InputTag("g4SimHits","MuonRPCHits")
#)

process.p = cms.Path(process.demo)
