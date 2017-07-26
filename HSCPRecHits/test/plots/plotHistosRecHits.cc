#include "plotHistosRecHits.h"
#include <fstream>
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
  styleTemplate->SetOptStat("nemruo");
  TString infile = "HSCP_RecHits.root";
  TString infile2 = "HSCP_MuTrigger_RecHits.root";
  TString infile3 = "ZMM_RecHits.root";
  
  double tmin= -100.; 
  double tmax= 100.;
  TFile* fin = new TFile(infile);
  TFile* fin2 = new TFile(infile2);
  TFile* fin3 = new TFile(infile3);
  
  int W = 1000;
  int H = 1000;
  float T = 0.08*H;
  float B = 0.12*H; 
  float L = 0.12*W;
  float R = 0.04*W;
  TCanvas * c1 = new TCanvas("c1", "c1",50,50,W,H);
  TH1D *fHBeta_pas = new TH1D("fHBeta_pas","fHBeta_pas",50,0.,1.);
  TH1D *fHBeta_tot =  new TH1D("fHBeta_tot","fHBeta_tot",50,0.,1.);
  TH1D *fHBeta_pas_BetaErrorCut = new TH1D("fHBeta_pas_BetaErrorCut","fHBeta_pas_BetaErrorCut",50,0.,1.);
  TH1D *fHBeta_pas_SlopeCut = new TH1D("fHBeta_pas_SlopeCut","fHBeta_pas_SlopeCut",50,0.,1.); 
  TH1D *fHBeta_resolution = new TH1D("fHBeta_resolution","fHBeta_resolution", 60, -3.,3.); 
  TH1D *fHBeta_resolution_bx = new TH1D("fHBeta_resolution_bx","fHBeta_resolution_bx", 60, -3.,3.); 
  TH1D *fHt0 = new TH1D("fHt0","fHt0",50,tmin,tmax); 
  TH1D *fHtime = new TH1D("fHtime","fHtime", 60, -10.,50.);
  TH1D *fHtimeNewFit = new TH1D("fHtimeNewFit","fHtimeNewFit",50, tmin,tmax);
  TH1D *fHtimeNewFit2 = new TH1D("fHtimeNewFit2","fHtimeNewFit2",50, tmin,tmax);
  TH1D *fHEta_pas = new TH1D("fHEta_pas", "fHEta_pas", 50,-3.15,3.15); 
  TH1D *fHEta_tot = new TH1D("fHEta_tot","fHEta_tot", 50,-3.15,3.15); 
  TH1D *fHbx = new TH1D("fHbx","fHbx", 300, -150., 150.); 

  TH2D *fHrpcBetaVsGenEta_pas = new TH2D("fHrpcBetaVsGenEta_pas","fHrpcBetaVsGenEta_pas",20, 0., 1.,10, 0.,2.0);
  TH2D *fHrpcBetaVsGenEta_tot = new TH2D("fHrpcBetaVsGenEta_tot","fHrpcBetaVsGenEta_tot",20, 0., 1.,10, 0.,2.0);
  TH2D *fHrpcBetaVsGenEta_eff = new TH2D("fHrpcBetaVsGenEta_eff","fHrpcBetaVsGenEta_eff",20, 0., 1.,10, 0.,2.0);
  
  TH1D * fHBeta_MuTrig_pas = (TH1D*) fin2->Get("demo2/fHbeta_pas");
  TH1D * fHBeta_MuTrig_tot = (TH1D*) fin2->Get("demo2/fHbeta_tot");
  
  TH1D *fHBeta_Zmm_pas = new TH1D("fHBeta_Zmm_pas","fHBeta_Zmm_pas",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_pas");
  TH1D *fHBeta_Zmm_tot = new TH1D("fHbeta_Zmm_tot","fHBeta_Zmm_tot",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_tot");
  TH1D *fHtime_Zmm = new TH1D("fHtime_Zmm","fHtime_Zmm",60,-10.,50.); //(TH1D*)fin3->Get("demo2/fHrpcHitsHitTime");
  
  unsigned short rpcHits_n, rpcHits_n_mu;
  Double_t rpcHitTime[20], rpcHitTime_mu[20],rpcHitTimeErr[20],rpcHitPos[20],rpcHitPosErr[20];
  double rpcBeta, rpcBetaErr,rpcBeta_bx, t0, t0_bx, genBeta, genEta,fitSlope;
  double rpcBeta_mu, rpcBetaErr_mu, t0_mu, t0_bx_mu, genBeta_mu, genEta_mu,fitSlope_mu;
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
  tree->SetBranchAddress("rpcBeta_bx",&rpcBeta_bx); 
  tree->SetBranchAddress("rpcBetaErr",&rpcBetaErr);
  tree->SetBranchAddress("t0",&t0);
  tree->SetBranchAddress("t0_bx",&t0_bx);
  tree->SetBranchAddress("genBeta",&genBeta);
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

  for(Long64_t i = 0 ; i < nEntries ; i++)
  {
    tree->GetEntry(i);
    if( fabs(genBeta) < = 2.4 )
    {
      vector<double> fitRes = fitWithErrors(rpcHits_n, rpcHitTime,rpcHitTimeErr,rpcHitPos,rpcHitPosErr);
      fHtimeNewFit->Fill(fitRes[0]);
      double res = (genBeta-rpcBeta)/genBeta; 
      fHBeta_resolution->Fill(res);
      double res_bx = (genBeta-rpcBeta_bx)/genBeta;
      fHBeta_resolution_bx->Fill(res_bx);
      if(fitSlope > 0 && rpcBetaErr/rpcBeta < 0.3 && rpcBeta<0.7 && rpcHits_n >=4)
      {
        fHBeta_pas->Fill(genBeta);
        fHEta_pas->Fill(genEta);
        fHrpcBetaVsGenEta_pas->Fill(genBeta, TMath::Abs(genEta));
      }
      fHrpcBetaVsGenEta_tot->Fill(genBeta, TMath::Abs(genEta));
      
      if(rpcBetaErr/rpcBeta < 0.3 && rpcBeta<0.7)
      {
        fHBeta_pas_BetaErrorCut->Fill(genBeta);
      }
      if(fitSlope > 0 && rpcBeta<0.7)
      {
        fHBeta_pas_SlopeCut->Fill(genBeta);
      }
      fHBeta_tot->Fill(genBeta);
      fHEta_tot->Fill(genEta);
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
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin( L/W );
  c1->SetRightMargin( R/W );
  c1->SetTopMargin( T/H );
  c1->SetBottomMargin( B/H );
  c1->SetTickx(0);
  c1->SetTicky(0);
  //c1->SetRightMargin(.12);
  styleTemplate->SetOptStat(0);
  styleTemplate->SetPalette(51);
  fHrpcBetaVsGenEta_pas->Sumw2();
  fHrpcBetaVsGenEta_tot->Sumw2();
  fHrpcBetaVsGenEta_pas->Divide(fHrpcBetaVsGenEta_tot);
  fHrpcBetaVsGenEta_pas->GetXaxis()->SetTitle("#beta");
  fHrpcBetaVsGenEta_pas->GetYaxis()->SetTitle("|#eta|");
  c1->SetRightMargin( .15 ); 
  fHrpcBetaVsGenEta_pas->Draw("colz");
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("rpcBetaVsGenEta.pdf");
  c1->SaveAs("rpcBetaVsGenEta.png");
  c1->SetRightMargin( R/W );
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
  
  TEfficiency* muTrigEff = new TEfficiency(*fHBeta_MuTrig_pas,*fHBeta_MuTrig_tot);
  muTrigEff->SetTitle(";#beta_{GEN};Efficiency");
  muTrigEff->SetMarkerColor(kBlue+2);
  muTrigEff->SetLineColor(kBlue+2);
  muTrigEff->SetName("muTrigEff");
  muTrigEff->Draw();
  gPad->Update();
  auto gMu = muTrigEff->GetPaintedGraph();
  gMu->SetMinimum(0.);
  gMu->SetMaximum(1.);
  gMu->GetXaxis()->SetRangeUser(0.,1.);
  gPad->Update();
  TLegend *lg = new TLegend(0.12,0.80,0.42,0.87);
  lg->AddEntry(muTrigEff,"HLT_L1SingleMuOpen_v3 trigger","lep");
  lg->SetTextSize(0.03);
  lg->SetFillColor(0);
  lg->SetFillStyle(0);
  lg->SetLineColor(1);
  lg->SetBorderSize(0);
  lg->Draw();
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("trigEff-HSCPTrigger.pdf");
  c1->SaveAs("trigEff-HSCPTrigger.png");
   
  trigEff->Draw();
  TLegend *lg2 = new TLegend(0.15,0.15,0.45,0.22);
  lg2->AddEntry(trigEff,"HSCP trigger","lep");
  lg2->SetTextSize(0.03);
  lg2->SetFillColor(0);
  lg2->SetFillStyle(0);
  lg2->SetLineColor(1);
  lg2->SetBorderSize(0);
  lg2->Draw();
  gPad->Update();
  auto g = trigEff->GetPaintedGraph();
  g->SetMinimum(0.);
  g->SetMaximum(1.);
  g->GetXaxis()->SetRangeUser(0.,1.);
  gPad->Update();

  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("trigEff-MuTrigger.pdf");
  c1->SaveAs("trigEff-MuTrigger.png");



  muTrigEff->Draw();
  trigEff->Draw("SAME");
  TLegend *lg3 = new TLegend(0.13,0.65,0.4,0.72);
  lg3->AddEntry(trigEff,"HSCP trigger","lep");
  lg3->AddEntry(muTrigEff,"HLT_L1SingleMuOpen_v3 trigger","lep");
  lg3->SetTextSize(0.03);
  lg3->SetFillColor(0);
  //lg3->SetFillStyle(0);
  lg3->SetLineColor(1);
  lg3->SetBorderSize(0);
  lg3->Draw();
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("trigEff-Mu-HSCPTriggers.pdf");
  c1->SaveAs("trigEff-Mu-HSCPTriggers.png");



   TEfficiency* trigEff2 = new TEfficiency(*fHBeta_pas_BetaErrorCut,*fHBeta_tot);
   trigEff2->SetTitle(";#beta_{GEN};Efficiency");
   trigEff2->SetMarkerColor(kOrange+2);
   trigEff2->SetLineColor(kOrange+2);
   trigEff2->SetName("trigEff2");
   trigEff2->Draw();
   TLegend *legend2 = new TLegend(0.4,0.25,0.65,0.32);
   legend2->AddEntry(trigEff,"HSCP trigger - Only cut in Fit Error.","lep");
   legend2->SetTextSize(0.025);
   //legend2->SetTextColor();
   legend2->SetLineColor(1);
   legend2->SetBorderSize(0);
   legend2->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEff-HSCPtrigger-OnlyFitErrorCut.pdf"); 
   c1->SaveAs("trigEff-HSCPtrigger-OnlyFitErrorCut.png"); 
   c1->cd();

   TEfficiency* trigEff3 = new TEfficiency(*fHBeta_pas_SlopeCut,*fHBeta_tot);
   plotEff(trigEff3);
   c1->cd();

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
   TLegend *legendEta = new TLegend(0.4,0.25,0.65,0.32);
   legendEta->AddEntry(trigEffEta,"HSCP trigger","lep");
   legendEta->SetTextSize(0.025);
   //legend->SetTextColor();
   legendEta->SetLineColor(1);
   legendEta->SetBorderSize(0);
   legendEta->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEffEta-HSCPtrigger.pdf");
   c1->SaveAs("trigEffEta-HSCPtrigger.png");


	fHBeta_resolution->GetXaxis()->SetTitle("#beta resolution (#beta_{GEN}-#beta_{RPC})/#beta_{GEN}");
	fHBeta_resolution->SetLineColor(kGreen+2);
  fHBeta_resolution->SetLineWidth(4);
	fHBeta_resolution->DrawNormalized();
	
	fHBeta_resolution_bx->GetXaxis()->SetTitle("#beta resolution (#beta_{GEN}-#beta_{RPC})/#beta_{GEN}");
	fHBeta_resolution_bx->SetLineColor(kOrange+2);
	fHBeta_resolution_bx->SetLineWidth(4);
	fHBeta_resolution_bx->DrawNormalized("SAME");
  
  TLegend *lgRes = new TLegend(0.15,0.55,0.4,0.64);
  lgRes->AddEntry(fHBeta_resolution_bx,"Phase 1 - 25ns t res","f");
  lgRes->AddEntry(fHBeta_resolution,"Phase 2 - 1ns t res","f");
  lgRes->SetTextSize(0.03);
  lgRes->SetFillColor(0);
  lgRes->SetFillStyle(0);
  lgRes->SetLineColor(1);
  lgRes->SetBorderSize(0);
  lgRes->Draw();
	
  CMS_lumi( c1, iPeriod, iPos );
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
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("timeZeroAndBx.pdf");
  c1->SaveAs("timeZeroAndBx.png");
  c1->SetLogy();
  c1->SaveAs("timeZeroAndBx_log.pdf");
  c1->SaveAs("timeZeroAndBx_log.png");
  
  c1->SetLogy(0);
  fHtime_Zmm->GetXaxis()->SetTitle("Time [ns]");
  fHtime->SetLineColor(kPink+2);
//  fHtime->SetFillStyle(3444);
//  fHtime->SetFillColor(kPink+2);
  fHtime_Zmm->SetLineColor(kBlue+2);
  fHtime_Zmm->SetLineWidth(4);
  fHtime->SetLineWidth(4);
  fHtime_Zmm->DrawNormalized();
  fHtime->DrawNormalized("SAME");
  TLegend *lgTime = new TLegend(0.55,0.55,0.9,0.64);
  lgTime->AddEntry(fHtime_Zmm,"Muon","f");
  lgTime->AddEntry(fHtime,"#tilde{#tau}","f");
  lgTime->SetTextSize(0.03);
  lgTime->SetFillColor(0);
  lgTime->SetFillStyle(0);
  lgTime->SetLineColor(1);
  lgTime->SetBorderSize(0);
  lgTime->Draw();
  CMS_lumi( c1, iPeriod, iPos );
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
   gPad->Update();
   auto graph2 = trigEff_Zmm->GetPaintedGraph();
   graph2->SetMinimum(0.);
   graph2->SetMaximum(1.);
   graph2->GetXaxis()->SetLimits(0.0,1.0);
   gPad->Update();
   lumi_13TeV_MC = "MC Pythia6 - ZMM - NoPU" ;
   TLegend *legend4 = new TLegend(0.4,0.25,0.65,0.32);
   legend4->AddEntry(trigEff_Zmm,"HSCP trigger","lep");
   legend4->SetTextSize(0.025);
   //legend4->SetTextColor();
   legend4->SetLineColor(1);
   legend4->SetBorderSize(0);
   legend4->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEff_ZMM-HSCPtrigger.pdf");
   c1->SaveAs("trigEff_ZMM-HSCPtrigger.png");
   c1->cd();
  c1->SetLogy(0);
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
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("time_Zmm.pdf");
  c1->SaveAs("time_Zmm.png");
	return 0;
}
