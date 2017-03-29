#include "styleTemplate.cc"
#include "CMS_lumi.cc"
#include "Utils.cc"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TPaveStats.h"
#include "TPoint.h"
#include "TCanvas.h"


TH1F* fHBetaTot = new TH1F("fHBetaTot","fHBetaTot",100,0.,1.); 
TH1F* fHBetaSel = new TH1F("fHBetaSel","fHBetaSel",100,0.,1.);
TH1F* fHBetaSimRes = new TH1F("fHBetaSimRes","fHBetaSimRes",100,-3.,3.);
TH1F* fHBetaGenRes = new TH1F("fHBetaGenRes","fHBetaGenRes",100,-3.,3.);
TH2F* fHist2DBx_Station = new TH2F("fHist2DBx_Station", "fHist2DBx_Station", 6, 1, 7, 20, 0, 20);
TH2F* fHist2DBetaVsEta = new TH2F("fHist2DBetaVsEta", "fHist2DBetaVsEta", 40, 0., 1., 40, -3.1416, 3.1416);

TString infile = "m1599_Gate1_BX_ToF_RecHits";
TString intreeRecHit = "demo2/rechitTree";
TString intreeEff = "demo2/effTree";
//TString outName = "Bx";

Int_t bunchX=0;
UInt_t stationhit=0;

Float_t beta=0;
Float_t betaSim=0;
Float_t betaGen=0;
Float_t etaGen=0;
UInt_t isChosen=0;


Int_t wCanvas=600;
Int_t hCanvas=600;
