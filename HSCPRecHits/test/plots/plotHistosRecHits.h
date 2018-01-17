#ifndef plotHistosRecHits_h
#define plotHistosRecHits_h
#include <iostream>
#include "styleTemplate.cc"
#include "CMS_lumi.cc"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TPaveStats.h"
#include "TPoint.h"
#include "TCanvas.h"
using namespace std;  


//TString infile = "m1599_RecHits";
TString infile = "rootFiles/HSCP_Rechits.root";
TString infile2 = "rootFiles/HSCP_MuTrigger_RecHits.root";
TString infile3 = "rootFiles/ZMM_RecHits.root";
  
TFile* fin = new TFile(infile);
TFile* fin2 = new TFile(infile2);
TFile* fin3 = new TFile(infile3);

TH1D *fHBeta_pas = new TH1D("fHBeta_pas","fHBeta_pas",50,0.,1.);
TH1D *fHBeta_tot =  new TH1D("fHBeta_tot","fHBeta_tot",50,0.,1.);
TH1D *fHBeta_pas_BetaErrorCut = new TH1D("fHBeta_pas_BetaErrorCut","fHBeta_pas_BetaErrorCut",50,0.,1.);
TH1D *fHBeta_pas_SlopeCut = new TH1D("fHBeta_pas_SlopeCut","fHBeta_pas_SlopeCut",50,0.,1.); 
TH1D *fHBeta_resolution = new TH1D("fHBeta_resolution","fHBeta_resolution", 60, -3.,3.); 
TH1D *fHt0 = new TH1D("fHt0","fHt0",50,-25., 25.); 
TH1D *fHtime = new TH1D("fHtime","fHtime", 60, -10.,50.);
TH1D *fHtimeNewFit = new TH1D("fHtimeNewFit","fHtimeNewFit",50, -25.,25.);
TH1D *fHEta_pas = new TH1D("fHEta_pas", "fHEta_pas", 50,-3.15,3.15); 
TH1D *fHEta_tot = new TH1D("fHEta_tot","fHEta_tot", 50,-3.15,3.15); 
TH1D *fHbx = new TH1D("fHbx","fHbx", 300, -150., 150.); 

TH1D * fHBeta_MuTrig_pas = (TH1D*) fin2->Get("demo2/fHbeta_pas");
TH1D * fHBeta_MuTrig_tot = (TH1D*) fin2->Get("demo2/fHbeta_tot");

TH1D *fHBeta_Zmm_pas = new TH1D("fHBeta_Zmm_pas","fHBeta_Zmm_pas",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_pas");
TH1D *fHBeta_Zmm_tot = new TH1D("fHbeta_Zmm_tot","fHBeta_Zmm_tot",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_tot");
TH1D *fHtime_Zmm = new TH1D("fHtime_Zmm","fHtime_Zmm",60,-10.,50.); //(TH1D*)fin3->Get("demo2/fHrpcHitsHitTime");

unsigned short rpcHits_n, rpcHits_n_mu;
Double_t rpcHitTime[20], rpcHitTime_mu[20],rpcHitTimeErr[20],rpcHitPos[20],rpcHitPosErr[20];
double rpcBeta, rpcBetaErr, t0, t0_bx, genBeta, genEta,fitSlope;
double rpcBeta_mu, rpcBetaErr_mu, t0_mu, t0_bx_mu, genBeta_mu, genEta_mu,fitSlope_mu;
void setBranches()
{
TTree *tree = (TTree*)fin->Get("demo2/tree");
TTree *tree_mu = (TTree*)fin3->Get("demo2/tree");
  tree_mu->SetBranchAddress("rpcHits_n",&rpcHits_n_mu);
  tree_mu->SetBranchAddress("rpcHitTime", rpcHitTime_mu);
  tree_mu->SetBranchAddress("rpcBeta",&rpcBeta_mu); 
  tree_mu->SetBranchAddress("rpcBetaErr",&rpcBetaErr_mu);
  tree_mu->SetBranchAddress("t0",&t0_mu);
  tree_mu->SetBranchAddress("t0_bx",&t0_bx_mu);
  tree_mu->SetBranchAddress("genBeta",&genBeta_mu);
  tree_mu->SetBranchAddress("genEta",&genEta_mu);
  tree_mu->SetBranchAddress("fitSlope",&fitSlope_mu);
  
  tree->SetBranchAddress("rpcHits_n",&rpcHits_n);
  tree->SetBranchAddress("rpcHitTime", rpcHitTime);
  tree->SetBranchAddress("rpcHitTimeErr", rpcHitTimeErr);
  tree->SetBranchAddress("rpcHitPos", rpcHitPos);
  tree->SetBranchAddress("rpcHitPosErr", rpcHitPosErr);
  tree->SetBranchAddress("rpcBeta",&rpcBeta); 
  tree->SetBranchAddress("rpcBetaErr",&rpcBetaErr);
  tree->SetBranchAddress("t0",&t0);
  tree->SetBranchAddress("t0_bx",&t0_bx);
  tree->SetBranchAddress("genBeta",&genBeta);
  tree->SetBranchAddress("genEta",&genEta);
  tree->SetBranchAddress("fitSlope",&fitSlope);
}
#endif
