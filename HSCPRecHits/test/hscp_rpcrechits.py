import FWCore.ParameterSet.Config as cms
import sys, os

from Configuration.StandardSequences.Eras import eras

process = cms.Process("demo2", eras.Phase2C1)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D12_cff')
process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
#from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring (
          #'/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCPppstau_M_1599_NoPU.root'
          'file:/afs/cern.ch/user/g/garamire/work/private/HSCP_stau-m1599_NoPU/HSCPppstau_M_1599_NoPU.root',
          )
)
process.load('test.HSCPRecHits.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   	fileName = cms.string("m1599_Gate1_BX_ToF_RecHits_test.root")
							)

#rpcRecHits.rpcDigiLabel = "simMuonRPCDigis"
#process.p = cms.Path(process.rpcRecHits*process.demo2)
process.p = cms.Path(process.demo2)
