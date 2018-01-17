import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import sys, os

from Configuration.StandardSequences.Eras import eras

process = cms.Process("demo2", eras.Phase2C1)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D12_cff')
process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')

#mylist = FileUtils.loadListFromFile('list_test.txt');
#mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/CMSSW_9_1_1_patch1/src/HSCPAnalysis/HSCPRecHits/test/list.txt');
#mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/CMSSW_9_1_1_patch1/src/HSCPAnalysis/HSCPRecHits/test/list_1218.txt');
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/CMSSW_9_1_1_patch1/src/HSCPAnalysis/HSCPRecHits/test/list_final.txt');
#mylist = FileUtils.loadListFromFile('list_HSCP_Official.txt');
#mylist= FileUtils.loadListFromFile('list_ZMM_PU200.txt');

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
#from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(*mylist)
)
process.load('HSCPAnalysis.HSCPRecHits.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   	#fileName = cms.string("ZMM_RecHits.root")
                    fileName = cms.string("HSCP_1599_PU200_final.root")
							)

#rpcRecHits.rpcDigiLabel = "simMuonRPCDigis"
#process.p = cms.Path(process.rpcRecHits*process.demo2)
process.p = cms.Path(process.demo2)
