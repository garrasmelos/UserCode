import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
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
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_15.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_16.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_17.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_18.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_19.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_20.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_21.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_22.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_23.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_24.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_25.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_26.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_27.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_28.root',
         '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Jan2017/HSCP_MC_GEN-SIM-Jan2017/170124_092652/0000/step1_29.root'
   )
)
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("hTest.root")

)
process.demo = cms.EDAnalyzer('HSCPanalyzer'
)


process.p = cms.Path(process.demo)
