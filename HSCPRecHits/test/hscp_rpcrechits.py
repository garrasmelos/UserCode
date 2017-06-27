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
          #'/store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D12PU200-v1/10000/0C535B13-8433-E711-ADF6-0CC47A4C8E16.root',
          #'/store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D12PU200-v1/10000/16FCA1D6-8533-E711-A16C-0025905A60D0.root',
          #'/store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D12PU200-v1/10000/1CE54D15-8E33-E711-AE39-0CC47A7C356A.root',
         # ' /store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D12PU200-v1/10000/2EA69035-7133-E711-894D-0CC47A4C8ECA.root',
          #'/store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D12PU200-v1/10000/30A61AF5-B633-E711-9BCC-0CC47A4D7662.root'
          
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/14A91939-1D3F-E711-A09A-0025905A610A.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/1AF43477-223F-E711-89EA-0025905A610C.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/4C4B7246-1E3F-E711-949B-0025905A60DA.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/722FED43-1E3F-E711-81B6-0CC47A78A3EC.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/7671FDAA-1F3F-E711-B7DE-0025905B860C.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/8AFE90A8-1F3F-E711-869F-0CC47A7C34E6.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/96206C74-223F-E711-9EBC-0CC47A7C3458.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/AC5F5FF7-1A3F-E711-9F95-0CC47A4D75EC.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/ACDADD48-203F-E711-ABC1-0CC47A4D75EC.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/B09C3A31-1F3F-E711-B70F-0CC47A78A3D8.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/DCC47538-1F3F-E711-888A-0CC47A78A33E.root',
#          '/store/relval/CMSSW_9_1_1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/DEAB0FD1-1C3F-E711-83A4-0025905A607A.root'

#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/0008B5BF-B34A-E711-820B-0CC47A4D762E.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/24733ABB-B24A-E711-8063-0CC47A4D7600.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/282BF5FE-B14A-E711-95ED-0025905A60DE.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/3E78F52F-B94A-E711-BEE4-0025905B859A.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/407AC7E2-B74A-E711-8E32-0CC47A4D76A0.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/46C6691B-B34A-E711-8ED2-0CC47A4D7614.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/54128E24-B34A-E711-AAEB-0025905B85AE.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/6278CB22-B34A-E711-B663-0025905A6084.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/76A3C5EF-B64A-E711-90F1-0025905A60EE.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/8055EF21-B34A-E711-B4CE-0025905A6122.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/D4047787-B34A-E711-B04D-0CC47A78A2EC.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/D8B866CB-B34A-E711-B9AF-0025905B8580.root',
#         '/store/relval/CMSSW_9_1_1_patch1/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/ECD40D26-B94A-E711-BD53-0CC47A7AB7A0.root'
          #'/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCPppstau_M_1599_NoPU.root'
          
          'file:/afs/cern.ch/user/g/garamire/work/private/HSCP_stau-m1599_NoPU/HSCPppstau_M_1599_NoPU.root',
          
          #'file:/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/data-HSCPppstau_M_1599_14TeV_GEN-SIM-RECO-DQMIO_pythia6_LG25-Full/HSCPstau-Pythia6-14GeV-LBUpgrade.root',
          #'file:/afs/cern.ch/user/g/garamire/work/private/RPCserviceWork/timingStudies/borisSimulation/CMSSW_9_1_0/src/mcGen/step3/step3.root'
          #'/store/mc/PhaseIIFall16DR82/HSCPppstau_M_871_TuneCUETP8M1_14TeV_pythia8/AODSIM/PU200_90X_upgrade2023_realistic_v1-v2/90000/02E4DD49-D9FD-E611-A063-02163E01449D.root'
          )
)
process.load('HSCPAnalysis.HSCPRecHits.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   	fileName = cms.string("HSCP_RecHits.root")
                    #fileName = cms.string(HSCP_MuTrigger_RecHits.root")
                    #fileName = cms.string("HSCP_RecHits.root")
							)

#rpcRecHits.rpcDigiLabel = "simMuonRPCDigis"
#process.p = cms.Path(process.rpcRecHits*process.demo2)
process.p = cms.Path(process.demo2)
