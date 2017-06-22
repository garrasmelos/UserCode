import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("demo2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

#process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
#from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

mylist = FileUtils.loadListFromFile('list.txt')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring ('file:/eos/cms/store/user/federica/HSCP_Simulation/HSCP_mass1599_AODSIM/step2_1_1.root')
    fileNames = cms.untracked.vstring (*mylist)
)
process.load('HSCPAnalysis.HSCPTrigger.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   	fileName = cms.string("m1599_MuTrigger.root")
							)

#rpcRecHits.rpcDigiLabel = "simMuonRPCDigis"
#process.p = cms.Path(process.rpcRecHits*process.demo2)
process.p = cms.Path(process.demo2)
