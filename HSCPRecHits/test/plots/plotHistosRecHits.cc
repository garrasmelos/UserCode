#include "plotHistosRecHits.h"
#include "TDRStyle.h"
#include <fstream>
const double relError = 100.;
const double maxEta = 2.4;
const double maxBeta = 0.7;//0.7;
const int nMinHits = 3; 
const double minSlope = 0.;

void drawCMS(TPad *pad);
void drawPU(TPad *pad);
void plotEff( TEfficiency * eff)
{
	setStyleTemplate();
  int iPeriod = 100;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  int iPos    = 11;
  
  TCanvas * c1Eff = new TCanvas("c1Eff", "c1Eff", 1000, 1000);
  TString lumi_13TeV_MC = "MC Pythia6 - #tilde{#tau}(M1599 GeV) of GMSB model" ;
  eff->SetTitle(";#beta_{GEN};Efficiency");
  eff->SetMarkerColor(kGreen+2);
  eff->SetLineColor(kGreen+2);
  eff->SetName("trigEff");
  eff->Draw();
  TLegend *legEff = new TLegend(0.4,0.25,0.65,0.32);
  legEff->AddEntry(eff->GetName(),"lep");
  legEff->SetTextSize(0.025);
  legEff->SetLineColor(1);
  legEff->SetBorderSize(0);
  legEff->Draw();
  CMS_lumi( c1Eff, iPeriod, iPos );
  //drawCMS(c1);
  //drawPU(c1);
  c1Eff->SaveAs("trigEff-HSCPtrigger-OnlySlopeCut.pdf");
  c1Eff->SaveAs("trigEff-HSCPtrigger-OnlySlopeCut.png");
}

vector<double > fitWithErrors(int nhits, double hitTime[20],double hitTimeErr[20],double hitPos[20],double hitPosErr[20])
{
  vector <double> results;
  double hitPos_tmp[nhits],hitPosErr_tmp[nhits],hitTime_tmp[nhits],hitTimeErr_tmp[nhits];
  int nHits=0;
  for(int j = 0; j< nhits; j++)
  {
    bool repeated=false;
    bool causal=true;
    for(int k =0; k<j; k++)
    { 
      if(hitPos[j]==hitPos[k] ) repeated = true;
      if(hitTime[k]>hitTime[j] ) causal=false;
    }
    if(!repeated && causal) 
    {
      hitPos_tmp[nHits]=hitPos[j];
      hitTime_tmp[nHits]=hitTime[j];
      hitPosErr_tmp[nHits]=hitPosErr[j];
      hitTimeErr_tmp[nHits]=hitTimeErr[j]; 
      nHits++;
    }
  }  
  TF1 * poly1 = new TF1("poly1","[0]+[1]*x", 4., 15.);
  poly1->SetParLimits(0,-950.,100.);
  poly1->SetParLimits(1,-10., 10.);
  poly1->SetParameters(0.0, 0.1);
  TGraphErrors *gr = new TGraphErrors(nhits,hitPos,hitTime,hitPosErr,hitTimeErr); 

  gr->Fit("poly1");
  double p0 = poly1->GetParameter(0);
  double p0Err = poly1->GetParError(0);
  double p1 = poly1->GetParameter(1);
  double p1Err = poly1->GetParError(1);
  
  results.push_back(p0);
  results.push_back(p0Err);
  results.push_back(p1);
  results.push_back(p1Err);
  //if(p0 < -950)
  //{
  //  gr2->Fit("poly1","W");
  //  gr2->Draw();
  //  c1->SaveAs("FitPlot.pdf");
  //  
  //  return 0;
  //}
  return results; 
}
int plotHistosRecHits()
{
	setStyleTemplate();
  setBranches();
  ofstream outfile("file.txt");
 // styleTemplate->SetOptStat("nemruo");
styleTemplate->SetOptFit(0);
TString infile = "rootFiles/HSCP_Rechits.root";
TString infile2 = "rootFiles/HSCP_MuTrigger_RecHits.root";
TString infile3 = "rootFiles/ZMM_RecHits.root";
//  TString infile4 = "/afs/cern.ch/work/j/jozobec/public/run2HSCP/ResultsICHEP16/Type2/Histos.root";
TString infile4 = "rootFiles/Histos.root";
TString infile5 = "rootFiles/HSCP_RecHits.root";
double tmin= -100.; 
double tmax= 100.;
TFile* fin = new TFile(infile);
TFile* fin2 = new TFile(infile2);
TFile* fin3 = new TFile(infile3);
TFile* fin4 = new TFile(infile4);
TFile* fin5 = new TFile(infile5);

int W = 800;
int H = 600;
float T = 0.08*H;
float B = 0.12*H; 
float L = 0.12*W;
float R = 0.04*W;
TCanvas * c1 = new TCanvas("c1", "c1",50,50,W,H);
TH1D *fHMass2 = new TH1D("fHMass2","fHMass2",300,0.,3000.);
TH1D *fHBeta_pas = new TH1D("fHBeta_pas","fHBeta_pas",50,0.,1.);
TH1D *fHBeta_tot =  new TH1D("fHBeta_tot","fHBeta_tot",50,0.,1.);
TH1D *fHBeta_pas_BetaErrorCut = new TH1D("fHBeta_pas_BetaErrorCut","fHBeta_pas_BetaErrorCut",50,0.,1.);
TH1D *fHBeta_pas_SlopeCut = new TH1D("fHBeta_pas_SlopeCut","fHBeta_pas_SlopeCut",50,0.,1.); 
TH1D *fHBeta_resolution = new TH1D("fHBeta_resolution","fHBeta_resolution", 60, -3.,3.); 
TH1D *fHBeta_resolution_bx = new TH1D("fHBeta_resolution_bx","fHBeta_resolution_bx", 60, -3.,3.); 
TH1D *fHt0 = new TH1D("fHt0","fHt0",50,tmin,tmax); 
TH1D *fHtime = new TH1D("fHtime","fHtime",  200 , -100., 100.);//60, -10.,50.);
TH1D *fHtimeNewFit = new TH1D("fHtimeNewFit","fHtimeNewFit",50, tmin,tmax);
TH1D *fHtimeNewFit2 = new TH1D("fHtimeNewFit2","fHtimeNewFit2",50, tmin,tmax);
TH1D *fHEta_pas = new TH1D("fHEta_pas", "fHEta_pas", 50,-3.15,3.15); 
TH1D *fHEta_tot = new TH1D("fHEta_tot","fHEta_tot", 50,-3.15,3.15); 
TH1D *fHbx = new TH1D("fHbx","fHbx", 300, -150., 150.); 

TH2D *fHrpcBetaVsGenEta_pas = new TH2D("fHrpcBetaVsGenEta_pas","fHrpcBetaVsGenEta_pas",20, 0., 1.,13, 0.,2.0);
TH2D *fHrpcBetaVsGenEta_tot = new TH2D("fHrpcBetaVsGenEta_tot","fHrpcBetaVsGenEta_tot",20, 0., 1.,13, 0.,2.0);
TH2D *fHrpcBetaVsGenEta_eff = new TH2D("fHrpcBetaVsGenEta_eff","fHrpcBetaVsGenEta_eff",20, 0., 1.,13, 0.,2.0);

TH1D * fHBeta_MuTrig_pas = (TH1D*) fin2->Get("demo2/fHbeta_pas");
TH1D * fHBeta_MuTrig_tot = (TH1D*) fin2->Get("demo2/fHbeta_tot");

TH1D *fHBeta_Zmm_pas = new TH1D("fHBeta_Zmm_pas","fHBeta_Zmm_pas",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_pas");
TH1D *fHBeta_Zmm_tot = new TH1D("fHbeta_Zmm_tot","fHBeta_Zmm_tot",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_tot");
TH1D *fHtime_Zmm = new TH1D("fHtime_Zmm","fHtime_Zmm", 200, -100.,100.);//60,-10.,50.); //(TH1D*)fin3->Get("demo2/fHrpcHitsHitTime");

TH1D *fHBeta_pas_nopu = new TH1D("fHBeta_pas_nopu","fHBeta_pas_nopu",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_pas");
TH1D *fHBeta_tot_nopu = new TH1D("fHbeta_tot_nopu","fHBeta_tot_nopu",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_tot");
TH1D *fHtime_nopu = new TH1D("fHtime_nopu","fHtime_nopu", 200, -100.,100.);//60,-10.,50.); //(TH1D*)fin3->Get("demo2/fHrpcHitsHitTime");

unsigned short rpcHits_n, rpcHits_n_mu, rpcHits_n_nopu;
Double_t rpcHitTime[100], rpcHitTime_mu[100], rpcHitTime_nopu[100], rpcHitTimeErr[100],rpcHitPos[100],rpcHitPosErr[100];
double rpcBeta, rpcBetaErr,rpcBeta_bx, t0, t0_bx, genBeta, genEta,fitSlope,mass;
double rpcBeta_mu, rpcBetaErr_mu, t0_mu, t0_bx_mu, genBeta_mu, genEta_mu,fitSlope_mu;
double rpcBeta_nopu, rpcBetaErr_nopu, t0_nopu, t0_bx_nopu, genBeta_nopu, genEta_nopu,fitSlope_nopu;
TTree *tree = (TTree*)fin->Get("demo2/tree");
TTree *tree_mu = (TTree*)fin3->Get("demo2/tree");
TTree *tree_nopu = (TTree*)fin5->Get("demo2/tree");
TH2D *fHMass1_2D = (TH2D*)fin4->Get("PPStau_13TeV16_M1599/MassTOF");
TH1D* fHMass1 = (TH1D*)fHMass1_2D->ProjectionY ("massAfterTightSelection", 300, 300);

tree_mu->SetBranchAddress("rpcHits_n",&rpcHits_n_mu);
tree_mu->SetBranchAddress("rpcHitTime", rpcHitTime_mu);
tree_mu->SetBranchAddress("rpcBeta",&rpcBeta_mu); 
tree_mu->SetBranchAddress("rpcBetaErr",&rpcBetaErr_mu);
tree_mu->SetBranchAddress("t0",&t0_mu);
tree_mu->SetBranchAddress("t0_bx",&t0_bx_mu);
tree_mu->SetBranchAddress("genBeta",&genBeta_mu);
tree_mu->SetBranchAddress("genEta",&genEta_mu);
tree_mu->SetBranchAddress("fitSlope",&fitSlope_mu);

tree_nopu->SetBranchAddress("rpcHits_n",&rpcHits_n_nopu);
tree_nopu->SetBranchAddress("rpcHitTime", rpcHitTime_nopu);
tree_nopu->SetBranchAddress("rpcBeta",&rpcBeta_nopu); 
tree_nopu->SetBranchAddress("rpcBetaErr",&rpcBetaErr_nopu);
tree_nopu->SetBranchAddress("t0",&t0_nopu);
tree_nopu->SetBranchAddress("t0_bx",&t0_bx_nopu);
tree_nopu->SetBranchAddress("genBeta",&genBeta_nopu);
tree_nopu->SetBranchAddress("genEta",&genEta_nopu);
tree_nopu->SetBranchAddress("fitSlope",&fitSlope_nopu);

tree->SetBranchAddress("rpcHits_n",&rpcHits_n);
tree->SetBranchAddress("rpcHitTime", rpcHitTime);
tree->SetBranchAddress("rpcHitTimeErr", rpcHitTimeErr);
tree->SetBranchAddress("rpcHitPos", rpcHitPos);
tree->SetBranchAddress("rpcHitPosErr", rpcHitPosErr);
tree->SetBranchAddress("rpcBeta",&rpcBeta); 
tree->SetBranchAddress("rpcBeta_bx",&rpcBeta_bx); 
tree->SetBranchAddress("rpcBetaErr",&rpcBetaErr);
tree->SetBranchAddress("t0",&t0);
tree->SetBranchAddress("t0_bx",&t0_bx);
tree->SetBranchAddress("genBeta",&genBeta);
tree->SetBranchAddress("mass",&mass);
tree->SetBranchAddress("genEta",&genEta);
tree->SetBranchAddress("fitSlope",&fitSlope);
 
bool writeExtraText = true;       // if extra text
TString extraText  = "Preliminary";  // default extra text is "Preliminary"
TString lumi_13TeV = "2.7 fb^{-1}";
TString lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
TString lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
int iPeriod = 100;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
int iPos    = 0;
int iPosR    = 33;
//TString lumi_13TeV_MC = "MC Pythia6 - #tilde{#tau}(M1599 GeV) of GMSB model" ;

//fin = new TFile(infile);
  //fin2 = new TFile(infile2);
  //fin3 = new TFile(infile3);

Long64_t nEntries = tree->GetEntries();
Long64_t nEntries_mu = tree_mu->GetEntries();
Long64_t nEntries_nopu = tree_nopu->GetEntries();

for(Long64_t i = 0 ; i < nEntries ; i++)
{
  tree->GetEntry(i);
  if(TMath::Abs(genEta)<maxEta)
  {
    //vector<double> fitRes = fitWithErrors(rpcHits_n, rpcHitTime,rpcHitTimeErr,rpcHitPos,rpcHitPosErr);
    fHtimeNewFit->Fill(0);//(fitRes[0]);
    double res = (genBeta-rpcBeta)/genBeta; 
    fHBeta_resolution->Fill(res);
    double res_bx = (genBeta-rpcBeta_bx)/genBeta;
    fHBeta_resolution_bx->Fill(res_bx);
    //if(rpcBeta > 0. && TMath::Abs(rpcBetaErr/rpcBeta) < relError && rpcBeta< maxBeta && rpcHits_n >=nMinHits )
    if(rpcBeta > 0. && rpcBeta< maxBeta && rpcHits_n >=nMinHits && rpcBetaErr < 10e9 )
    {
      fHBeta_pas->Fill(genBeta); 
      fHMass2->Fill(mass-200.);
      if(genBeta > 0.2 and genBeta <0.6)fHEta_pas->Fill(genEta);
      fHrpcBetaVsGenEta_pas->Fill(genBeta, TMath::Abs(genEta));
    }
    fHrpcBetaVsGenEta_tot->Fill(genBeta, TMath::Abs(genEta));
    
    if(rpcBetaErr/rpcBeta < relError  && rpcBeta<maxBeta )
    {
      fHBeta_pas_BetaErrorCut->Fill(genBeta);
    }
    if(fitSlope > minSlope && rpcBeta<maxBeta )
    {
      fHBeta_pas_SlopeCut->Fill(genBeta);
    }
    fHBeta_tot->Fill(genBeta);
    if(genBeta > 0.2 and genBeta <0.6) fHEta_tot->Fill(genEta);
    fHt0->Fill(t0);
    fHbx->Fill(t0_bx);
    //cout << rpcHitPos[1] << " " << rpcHits_n <<  endl;
  }
  for(UShort_t k=0 ; k < rpcHits_n ; k++)
  {
    fHtime->Fill(rpcHitTime[k]);     
  }  
}


for(Long64_t j = 0 ; j < nEntries_mu ; j++)
{
  tree_mu->GetEntry(j);
  fHBeta_Zmm_tot->Fill(genBeta_mu);
  
  if(fitSlope_mu > 0 && rpcBetaErr_mu/rpcBeta_mu < 0.3 && rpcBeta_mu<0.7)
  {
    fHBeta_Zmm_pas->Fill(genBeta_mu);
  }
  for(UShort_t k=0 ; k < rpcHits_n_mu ; k++)
  {
    fHtime_Zmm->Fill(rpcHitTime_mu[k]);
  }
}

for(Long64_t i = 0 ; i < nEntries_nopu ; i++)
{
  tree_nopu->GetEntry(i);
  if(TMath::Abs(genEta_nopu)<maxEta)
  {
    if(fitSlope_nopu > minSlope && TMath::Abs(rpcBetaErr_nopu/rpcBeta_nopu) < relError && rpcBeta_nopu< maxBeta && rpcHits_n_nopu >=nMinHits )
    {
      fHBeta_pas_nopu->Fill(genBeta_nopu); 
    }
    
    fHBeta_tot_nopu->Fill(genBeta_nopu);
  }
}
/*
for(Long64_t j = 0 ; j < nEntries_nopu ; j++)
{
  tree_nopu->GetEntry(j);
  fHBeta_tot_nopu->Fill(genBeta_nopu);
  
  if(fitSlope_nopu > 0 && rpcBetaErr_nopu/rpcBeta_nopu < 0.3 && rpcBeta_nopu<0.7)
  {
    fHBeta_pas_nopu->Fill(genBeta_nopu);
  }
  for(UShort_t k=0 ; k < rpcHits_n_nopu ; k++)
  {
    fHtime_nopu->Fill(rpcHitTime_nopu[k]);
  }
}
*/
//c1->SetFillColor(0);
//c1->SetBorderMode(0);
//c1->SetFrameFillStyle(0);
//c1->SetFrameBorderMode(0);

c1->SetLeftMargin( 0.15 );
c1->SetRightMargin( 0.16 );
c1->SetTopMargin(0.075);
c1->SetBottomMargin(0.125);
c1->SetTicks(1,1);
//c1->SetRightMargin(.12);
styleTemplate->SetOptStat(0);
  styleTemplate->SetOptTitle(0);
  styleTemplate->SetPalette(51);
  
  fHrpcBetaVsGenEta_pas->Sumw2();
  fHrpcBetaVsGenEta_tot->Sumw2();
  fHrpcBetaVsGenEta_pas->Divide(fHrpcBetaVsGenEta_tot);
  fHrpcBetaVsGenEta_pas->GetXaxis()->SetTitle("#beta");
  fHrpcBetaVsGenEta_pas->GetYaxis()->SetTitle("|#eta|");
  fHrpcBetaVsGenEta_pas->GetZaxis()->SetTitle("Efficiency");
  fHrpcBetaVsGenEta_pas->GetZaxis()->SetTitleOffset(0.91);
  fHrpcBetaVsGenEta_pas->SetMaximum(1.0);
  
  //c1->SetRightMargin( .15 ); 
  fHrpcBetaVsGenEta_pas->Draw("colz");
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("rpcBetaVsGenEta.pdf");
  c1->SaveAs("rpcBetaVsGenEta.png");
  //c1->SetRightMargin( R/W );
 
  c1->SetRightMargin( 0.05 );
  fHMass1->Rebin(6);
  fHMass2->Rebin(6);
  double integral1 = fHMass1->Integral();
  double integral2 = fHMass2->Integral();
  if(integral1!=0) fHMass1->Scale(1000./integral1);
  if(integral2!=0) fHMass2->Scale(1000./integral2);
  fHMass2->GetXaxis()->SetTitle("Mass [GeV]");
  fHMass2->GetYaxis()->SetTitle("Events (a.u.)");
  fHMass2->SetTitle("");
  TF1 *g1 = new TF1("g1","gaus",1000.,2000.);
  TF1 *g2 = new TF1("g2","gaus",1000.,2000.);
  fHMass1->Fit(g1,"R");
  fHMass2->Fit(g2,"R");
  fHMass1->GetFunction("g1")->SetLineColor(kBlue);
  Double_t par[6];
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);

  g1->SetLineColor(kBlue);
  cout << par[2] << endl;
  cout << par[5] << endl;

  fHMass1->SetLineColor(kBlue);
  fHMass1->SetMarkerColor(kBlue);
  fHMass2->SetLineColor(kRed);
  fHMass2->SetMarkerColor(kRed);
  formatHisto(fHMass1);
  formatHisto(fHMass2);

  fHMass2->Draw();
  fHMass1->Draw("SAME");
  drawCMS(c1);
  drawPU(c1);
  TLegend *lm = new TLegend(0.18,0.7,0.42,0.9);
  lm->AddEntry(fHMass2,"Phase 2","lep");
  lm->AddEntry((TObject*)0, Form("#sigma= %.0f GeV",par[5]) ,"");
  lm->AddEntry(fHMass1,"Run 2","lep");
  lm->AddEntry((TObject*)0, Form("#sigma= %.0f GeV",par[2]) ,"");
  lm->SetTextFont(42);
  lm->SetTextSize(0.045);
  lm->SetFillColor(0);
  lm->SetFillStyle(0);
  lm->SetLineColor(1);
  lm->SetBorderSize(0);
  lm->Draw();
  /*TLatex latm;
  latm.SetNDC();
  //latm.SetTextFont(62);
  latm.SetTextSize(0.025);
  latm.DrawLatex(0.18, 0.54, Form("#sigma_{1}= %f",par[2]));
  latm.DrawLatex(0.18, 0.44, Form("#sigma_{2}= %f",par[5]));
  */
  c1->SaveAs("mass.pdf");
  c1->SaveAs("mass.png");
  
  c1->SetRightMargin( 0.05 );
  
  fHt0->Draw();
  fHtimeNewFit->SetLineColor(kBlue);
  fHtimeNewFit->Draw("SAME");
  fHtimeNewFit2->SetLineColor(kRed);
  fHtimeNewFit2->Draw("SAME");
  c1->SaveAs("timeNewFit.pdf");
  c1->SaveAs("timeNewFit.png");

  TEfficiency* trigEff = new TEfficiency(*fHBeta_pas,*fHBeta_tot);
  trigEff->SetTitle(";#beta_{GEN};Efficiency");
  trigEff->SetMarkerColor(kRed+2);
  trigEff->SetLineColor(kRed+2);
  trigEff->SetName("trigEff");
  
  TEfficiency* trigEff_nopu = new TEfficiency(*fHBeta_pas_nopu,*fHBeta_tot_nopu);
  trigEff_nopu->SetTitle(";#beta_{GEN};Efficiency");
  trigEff_nopu->SetMarkerColor(kBlack);
  trigEff_nopu->SetLineColor(kBlack);
  trigEff_nopu->SetName("trigEff_nopu");
  trigEff_nopu->Draw();
  gPad->Update();
  
  TEfficiency* muTrigEff = new TEfficiency(*fHBeta_MuTrig_pas,*fHBeta_MuTrig_tot);
  muTrigEff->SetTitle(";#beta_{GEN};Efficiency");
  muTrigEff->SetMarkerColor(kBlue+2);
  muTrigEff->SetLineColor(kBlue+2);
  muTrigEff->SetName("muTrigEff");
  muTrigEff->Draw();
  gPad->Update();
  auto gMu = muTrigEff->GetPaintedGraph();
  gMu->SetMinimum(0.);
  gMu->GetXaxis()->SetTitleSize(0.05);
  gMu->GetYaxis()->SetTitleSize(0.05);
  gMu->GetYaxis()->SetLabelSize(0.04);
  gMu->GetXaxis()->SetLabelSize(0.04);
  gMu->SetMaximum(1.2);
  gMu->GetXaxis()->SetRangeUser(0.,1.);
  gPad->Update();
  TLegend *lg = new TLegend(0.12,0.80,0.42,0.87);
  lg->AddEntry(muTrigEff,"Phase-1 regular muon trigger (L1 Mu Open)","lep");
  lg->SetTextSize(0.03);
  lg->SetFillColor(0);
  lg->SetFillStyle(0);
  lg->SetLineColor(1);
  lg->SetBorderSize(0);
  lg->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEff-HSCPTrigger.pdf");
  c1->SaveAs("trigEff-HSCPTrigger.png");
   
  trigEff->Draw();
  TLegend *lg2 = new TLegend(0.15,0.15,0.45,0.22);
  lg2->AddEntry(trigEff,"RPC-HSCP trigger","lep");
  lg2->SetTextSize(0.03);
  lg2->SetFillColor(0);
  lg2->SetFillStyle(0);
  lg2->SetLineColor(1);
  lg2->SetBorderSize(0);
  lg2->Draw();
  gPad->Update();
  auto g = trigEff->GetPaintedGraph();
  g->SetMinimum(0.);
  g->SetMaximum(1.2);
  g->GetXaxis()->SetRangeUser(0.,1.);
  gPad->Update();

  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEff-MuTrigger.pdf");
  c1->SaveAs("trigEff-MuTrigger.png");



  muTrigEff->Draw();
  trigEff->Draw("SAME");
  trigEff_nopu->Draw("SAME");
  TLegend *lg3 = new TLegend(0.18,0.78,0.6,0.91);
  lg3->AddEntry(trigEff_nopu,"Phase-2 RPC-HSCP trigger (PU=0)","lep");
  lg3->AddEntry(trigEff,"Phase-2 RPC-HSCP trigger (PU=200)","lep");
  lg3->AddEntry(muTrigEff,"Phase-1 Regular #mu trigger (L1 Mu Open)","lep");
  lg3->SetTextFont(42);
  lg3->SetTextSize(0.035);
  lg3->SetFillColor(0);
  //lg3->SetFillStyle(0);
  lg3->SetLineColor(1);
  lg3->SetBorderSize(0);
  lg3->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEff-Mu-HSCPTriggers.pdf");
  c1->SaveAs("trigEff-Mu-HSCPTriggers.png");



  TEfficiency* trigEff2 = new TEfficiency(*fHBeta_pas_BetaErrorCut,*fHBeta_tot);
  trigEff2->SetTitle(";#beta_{GEN};Efficiency");
  trigEff2->SetMarkerColor(kOrange+2);
  trigEff2->SetLineColor(kOrange+2);
  trigEff2->SetName("trigEff2");
  trigEff2->Draw();
  TLegend *legend2 = new TLegend(0.4,0.25,0.65,0.32);
  legend2->AddEntry(trigEff,"RPC-HSCP trigger - Only cut in Fit Error.","lep");
  legend2->SetTextSize(0.025);
  //legend2->SetTextColor();
  legend2->SetLineColor(1);
  legend2->SetBorderSize(0);
  legend2->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEff-HSCPtrigger-OnlyFitErrorCut.pdf"); 
  c1->SaveAs("trigEff-HSCPtrigger-OnlyFitErrorCut.png"); 
  c1->cd();

  TEfficiency* trigEff3 = new TEfficiency(*fHBeta_pas_SlopeCut,*fHBeta_tot);
  plotEff(trigEff3);
  c1->cd();
  fHEta_tot->Draw();
  c1->SaveAs("Eta.pdf");
  
  TEfficiency* trigEffEta = new TEfficiency(*fHEta_pas,*fHEta_tot);
  trigEffEta->SetTitle(";#eta_{GEN};Efficiency");
  trigEffEta->SetMarkerColor(kRed+2);
  trigEffEta->SetLineColor(kRed+2);
  trigEffEta->SetName("trigEffEta");
  trigEffEta->Draw();
  gPad->Update();
  auto graph = trigEffEta->GetPaintedGraph();
  graph->SetMinimum(0.);
  graph->SetMaximum(1.);
  gPad->Update();
  TLegend *legendEta = new TLegend(0.18,0.45,0.45,0.52);
  legendEta->AddEntry(trigEffEta,"Phase-2 RPC-HSCP trigger","lep");
  legendEta->SetTextSize(0.035);
  //legend->SetTextColor();
  legendEta->SetLineColor(1);
  legendEta->SetBorderSize(0);
  legendEta->Draw();
  TLatex rl;
  rl.SetTextSize(0.035);
  rl.DrawLatexNDC(0.25,0.4,"0.2<#beta_{GEN}<0.6");
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEffEta-HSCPtrigger.pdf");
  c1->SaveAs("trigEffEta-HSCPtrigger.png");


	fHBeta_resolution->GetXaxis()->SetTitle("#beta resolution (#beta_{GEN}-#beta_{RPC})/#beta_{GEN}");
  fHBeta_resolution->GetXaxis()->SetTitleSize(0.05);
  fHBeta_resolution->GetYaxis()->SetTitleSize(0.05);
  fHBeta_resolution->GetYaxis()->SetLabelSize(0.04);
  fHBeta_resolution->GetXaxis()->SetLabelSize(0.04);
	fHBeta_resolution->SetLineColor(kGreen+2);
  fHBeta_resolution->SetLineWidth(4);
  fHBeta_resolution->GetYaxis()->SetTitle("Entries (a.u.)");
	fHBeta_resolution->DrawNormalized();
	
	fHBeta_resolution_bx->GetYaxis()->SetTitle("#beta resolution (#beta_{GEN}-#beta_{RPC})/#beta_{GEN}");
	fHBeta_resolution_bx->SetLineColor(kOrange+2);
	fHBeta_resolution_bx->SetLineWidth(4);
	fHBeta_resolution_bx->DrawNormalized("SAME");
  
  TLegend *lgRes = new TLegend(0.15,0.8,0.42,0.89);
  lgRes->AddEntry(fHBeta_resolution_bx,"Phase 1 - 25ns time resolution","f");
  lgRes->AddEntry(fHBeta_resolution,"Phase 2 - 1.5ns time resolution","f");
  lgRes->SetTextFont(42);
  lgRes->SetTextSize(0.035);
  lgRes->SetFillColor(0);
  lgRes->SetFillStyle(0);
  lgRes->SetLineColor(1);
  lgRes->SetBorderSize(0);
  lgRes->Draw();
	
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
	c1->SaveAs("beta_GenRes.pdf");
	c1->SaveAs("beta_GenRes.png");
  
  fHt0->GetXaxis()->SetTitle("Time [ns]");
  fHt0->SetLineColor(kAzure+1);
  fHt0->SetFillStyle(3444);
  //fHt0->SetFillColor(kAzure+2);
  fHt0->SetName("fHt0");
  //fHt0->Draw();
  fHbx->SetLineColor(kRed+1);
  fHbx->SetName("fHbx");
  fHbx->GetXaxis()->SetTitle("Time [ns]");
  fHbx->SetLineWidth(4);
  fHt0->SetLineWidth(4);
  fHbx->DrawNormalized();
  fHt0->DrawNormalized("SAME");
  TLegend *legendTime = new TLegend(0.28,0.7,0.42,0.8);
  legendTime->SetHeader("Production  BX identification","C");
  legendTime->AddEntry(fHbx,"Phase-1","lep");
  legendTime->AddEntry(fHt0,"Phase-2","lep");
  legendTime->SetTextSize(0.025);
  //legend->SetTextColor();
  legendTime->SetLineColor(1);
  legendTime->SetBorderSize(0);
  legendTime->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("timeZeroAndBx.pdf");
  c1->SaveAs("timeZeroAndBx.png");
  c1->SetLogy();
  c1->SaveAs("timeZeroAndBx_log.pdf");
  c1->SaveAs("timeZeroAndBx_log.png");
  
  //c1->SetLogy(0);
  fHtime_Zmm->GetXaxis()->SetTitle("Time [ns]");
  fHtime->SetLineColor(kPink+2);
//  fHtime->SetFillStyle(3444);
//  fHtime->SetFillColor(kPink+2);
  fHtime_Zmm->GetYaxis()->SetTitle("Entries (a.u.)");
  fHtime_Zmm->SetLineColor(kBlue+2);
  fHtime_Zmm->SetLineWidth(4);
  fHtime->SetLineWidth(4);
  //fHtime_Zmm->GetXaxis()->SetTitleSize(0.05);
  //fHtime_Zmm->GetYaxis()->SetTitleSize(0.05);
  //fHtime_Zmm->GetYaxis()->SetLabelSize(0.04);
  //fHtime_Zmm->GetXaxis()->SetLabelSize(0.04);
  //fHtime_Zmm->DrawNormalized();
  fHtime->DrawNormalized();
  TLegend *lgTime = new TLegend(0.33,0.55,0.78,0.64);
  lgTime->AddEntry(fHtime_Zmm,"#mu from ZMM sample Phase-2 geometry","f");
  lgTime->AddEntry(fHtime,"#tilde{#tau} GMSB m=1599GeV","f");
  lgTime->SetTextFont(42);
  lgTime->SetTextSize(0.042);
  lgTime->SetFillColor(0);
  lgTime->SetFillStyle(0);
  lgTime->SetLineColor(1);
  lgTime->SetBorderSize(0);
  //lgTime->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("time.pdf");
  c1->SaveAs("time.png");


  fHBeta_Zmm_pas->Draw();
  c1->SaveAs("beta_Zmm_pas.pdf");
  c1->SaveAs("beta_Zmm_pas.png");
  fHBeta_Zmm_tot->Draw();
  c1->SaveAs("beta_Zmm_tot.pdf");
  c1->SaveAs("beta_Zmm_tot.png");
  TEfficiency* trigEff_Zmm = new TEfficiency(*fHBeta_Zmm_pas,*fHBeta_Zmm_tot);
  trigEff_Zmm->SetTitle(";#beta_{GEN};Efficiency");
  trigEff_Zmm->SetMarkerColor(kBlack+2);
  trigEff_Zmm->SetLineColor(kBlack+2);
  trigEff_Zmm->SetName("trigEff_Zmm");
  trigEff_Zmm->Draw();
  auto graph2 = trigEff_Zmm->GetPaintedGraph();
  graph2->SetMinimum(0.);
  graph2->SetMaximum(1.);
  graph2->GetXaxis()->SetLimits(0.0,1.0);
  gPad->Update();
  lumi_13TeV_MC = "MC Pythia6 - ZMM - NoPU" ;
  TLegend *legend4 = new TLegend(0.4,0.25,0.65,0.32);
  legend4->AddEntry(trigEff_Zmm,"Phase-2 RPC-HSCP trigger","lep");
  legend4->SetTextSize(0.025);
  //legend4->SetTextColor();
  legend4->SetLineColor(1);
  legend4->SetBorderSize(0);
  legend4->Draw();
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("trigEff_ZMM-HSCPtrigger.pdf");
  c1->SaveAs("trigEff_ZMM-HSCPtrigger.png");
  c1->cd();
  c1->SetLogy(0);

  styleTemplate->SetOptStat("nemruo");
  fHtime_Zmm->GetXaxis()->SetTitle("Time [ns]");
  fHtime_Zmm->SetLineColor(kBlue+4);
  fHtime_Zmm->SetFillStyle(3444);
  fHtime_Zmm->SetFillColor(kBlue+4);
  fHtime_Zmm->Fit("gaus");
  fHtime_Zmm->Draw();
  gPad->Update();
  TPaveStats* sb2=(TPaveStats*)(fHtime_Zmm->FindObject("stats"));
  sb2->SetX1NDC(.62);
  sb2->SetX2NDC(.92);
  sb2->SetY1NDC(.7);
  sb2->SetY2NDC(.9);
  sb2->SetTextColor(1);
  gPad->Modified();
  fHtime_Zmm->GetXaxis()->SetRangeUser(-10.5,10.5);
  //CMS_lumi( c1, iPeriod, iPos );
  drawCMS(c1);
  drawPU(c1);
  c1->SaveAs("time_Zmm.pdf");
  c1->SaveAs("time_Zmm.png");
	return 0;
}
void drawCMS(TPad *pad)
{
  TString text="Preliminary";
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(62);
  latex.SetTextSize(0.0625);
  latex.DrawLatex(0.18, 0.94, "CMS");
  latex.SetTextSize(0.05);
  latex.SetTextFont(52);
  latex.DrawLatex(0.285, 0.94, text);
  return;
}
void drawPU(TPad *pad)
{
  TString text="14 TeV";
  //TString text="14 TeV, <PU>=200";
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.06);
  latex.SetTextAlign(31);
  latex.DrawLatex(0.95,0.94,text);
}
