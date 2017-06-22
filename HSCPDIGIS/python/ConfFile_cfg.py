import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_1.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_2.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_3.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_4.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_5.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_6.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_7.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_8.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_9.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_10.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_11.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_12.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_13.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_14.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_15.root'
    )
)
process.load('HSCPAnalysis.HSCPDIGIS.CfiFile_cfi')
process.TFileService = cms.Service("TFileService",
                   fileName = cms.string("hTest.root")
                   )
#process.demo = cms.EDAnalyzer('HSCPDIGIS'
#)


process.p = cms.Path(process.demo)
