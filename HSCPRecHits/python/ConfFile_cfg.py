import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Demo2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")

mylist = FileUtils.loadListFromFile('list.txt')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
        	fileNames = cms.untracked.vstring (*mylist)
)

process.load('HSCPAnalysis.HSCPRecHits.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
      	fileName = cms.string("m247_Gate1_BX_ToF_RecHits.root")
)

process.p = cms.Path(process.rpcRecHits*process.demo2)
