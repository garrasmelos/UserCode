#include "plotHistosRecHits.h"


int plotHistosRecHits()
{
	setStyleTemplate();
  styleTemplate->SetOptStat(1);
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_13TeV = "2.7 fb^{-1}";
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_13TeV_MC = "MC Pythia6 - #tilde{#tau}(M1599 GeV) of GMSB model" ;
  	
  int iPeriod = 100;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  int iPos    = 11;
	TString infile = "HSCP_RecHits.root";
  TString infile2 = "HSCP_MuTrigger_RecHits.root";
  TString infile3 = "ZMM_RecHits.root";
   
  TFile* fin = new TFile(infile);
  TFile* fin2 = new TFile(infile2);
  TFile* fin3 = new TFile(infile3);

  TCanvas * c1 = new TCanvas("c1", "c1", 1000, 1000);
  TTree *tree = (TTree*)fin->Get("demo2/tree");
  TTree *tree_mu = (TTree*)fin3->Get("demo2/tree");

  TH1D *fHBeta_pas = new TH1D("fHBeta_pas","fHBeta_pas",50,0.,1.);
  TH1D *fHBeta_tot =  new TH1D("fHBeta_tot","fHBeta_tot",50,0.,1.);
  TH1D *fHBeta_pas_BetaErrorCut = new TH1D("fHBeta_pas_BetaErrorCut","fHBeta_pas_BetaErrorCut",50,0.,1.); //(TH1D*)fin->Get("demo2/fHbeta_pas_BetaErrorCut");
  TH1D *fHBeta_pas_SlopeCut = new TH1D("fHBeta_pas_SlopeCut","fHBeta_pas_SlopeCut",50,0.,1.); //(TH1D*)fin->Get("demo2/fHbeta_pas_SlopeCut");
  TH1D *fHBeta_resolution = new TH1D("fHBeta_resolution","fHBeta_resolution", 60, -3.,3.); //(TH1D*)fin->Get("demo2/fHres");
  TH1D *fHt0 = new TH1D("fHt0","fHt0",300,-150., 150.); //(TH1D*)fin->Get("demo2/fHt0");
  TH1D *fHtime = new TH1D("fHtime","fHtime", 60, -10.,50.); //(TH1D*)fin->Get("demo2/fHrpcHitsTime");
  TH1D *fHEta_pas = new TH1D("fHEta_pas", "fHEta_pas", 50,-3.15,3.15); //(TH1D*)fin->Get("demo2/fHeta_pas");
  TH1D *fHEta_tot = new TH1D("fHEta_tot","fHEta_tot", 50,-3.15,3.15); //(TH1D*)fin->Get("demo2/fHeta_tot");
  TH1D *fHbx = new TH1D("fHbx","fHbx", 300, -150., 150.); //(TH1D*)fin->Get("demo2/fHt0_bx");

  TH1D * fHBeta_MuTrig_pas = (TH1D*) fin2->Get("demo2/fHbeta_pas");
  TH1D * fHBeta_MuTrig_tot = (TH1D*) fin2->Get("demo2/fHbeta_tot");

  TH1D *fHBeta_Zmm_pas = new TH1D("fHBeta_Zmm_pas","fHBeta_Zmm_pas",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_pas");
  TH1D *fHBeta_Zmm_tot = new TH1D("fHbeta_Zmm_tot","fHBeta_Zmm_tot",50,0.,1.); //(TH1D*)fin3->Get("demo2/fHbeta_tot");
  TH1D *fHtime_Zmm = new TH1D("fHtime_Zmm","fHtime_Zmm",60,-10.,50.); //(TH1D*)fin3->Get("demo2/fHrpcHitsTime");

  Long64_t nEntries = tree->GetEntries();
  Long64_t nEntries_mu = tree_mu->GetEntries();
  unsigned short rpcHits_n, rpcHits_n_mu;
  double rpcTime[20], rpcTime_mu[20];
  double rpcBeta, rpcBetaErr, t0, t0_bx, genBeta, genEta,fitSlope;
  double rpcBeta_mu, rpcBetaErr_mu, t0_mu, t0_bx_mu, genBeta_mu, genEta_mu,fitSlope_mu;

  tree_mu->SetBranchAddress("rpcHits_n",&rpcHits_n_mu);
  tree_mu->SetBranchAddress("rpcTime", rpcTime_mu);
  tree_mu->SetBranchAddress("rpcBeta",&rpcBeta_mu); 
  tree_mu->SetBranchAddress("rpcBetaErr",&rpcBetaErr_mu);
  tree_mu->SetBranchAddress("t0",&t0_mu);
  tree_mu->SetBranchAddress("t0_bx",&t0_bx_mu);
  tree_mu->SetBranchAddress("genBeta",&genBeta_mu);
  tree_mu->SetBranchAddress("genEta",&genEta_mu);
  tree_mu->SetBranchAddress("fitSlope",&fitSlope_mu);

  tree->SetBranchAddress("rpcHits_n",&rpcHits_n);
  tree->SetBranchAddress("rpcTime", rpcTime);
  tree->SetBranchAddress("rpcBeta",&rpcBeta); 
  tree->SetBranchAddress("rpcBetaErr",&rpcBetaErr);
  tree->SetBranchAddress("t0",&t0);
  tree->SetBranchAddress("t0_bx",&t0_bx);
  tree->SetBranchAddress("genBeta",&genBeta);
  tree->SetBranchAddress("genEta",&genEta);
  tree->SetBranchAddress("fitSlope",&fitSlope);

  for(Long64_t i = 0 ; i < nEntries ; i++)
  {
    tree->GetEntry(i);
    if(rpcHits_n >2)
    {
      if(fitSlope > 0 && rpcBetaErr/rpcBeta < 0.3 && rpcBeta<0.7)
      {
        double res = (genBeta-rpcBeta)/genBeta; 
        fHBeta_resolution->Fill(res);
        fHBeta_pas->Fill(genBeta);
        fHEta_pas->Fill(genEta);
      }
      
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
    }
    for(UShort_t k=0 ; k < rpcHits_n ; k++)
    {
      fHtime->Fill(rpcTime[k]);     
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
      fHtime_Zmm->Fill(rpcTime_mu[k]);
    }
  }
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
   trigEff->Draw("SAME");

   TLegend *legend = new TLegend(0.45,0.15,0.85,0.3);
   legend->AddEntry(trigEff,"HSCP trigger","lep");
   legend->AddEntry(muTrigEff,"HLT_L1SingleMuOpen_v3 trigger","lep");
   legend->SetTextSize(0.025);
   //legend->SetTextColor();
   legend->SetLineColor(1);
   legend->SetBorderSize(0);
   legend->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEff-HSCPandMuTriggers.pdf");

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
   c1->cd();

   TEfficiency* trigEff3 = new TEfficiency(*fHBeta_pas_SlopeCut,*fHBeta_tot);
   trigEff3->SetTitle(";#beta_{GEN};Efficiency");
   trigEff3->SetMarkerColor(kYellow+2);
   trigEff3->SetLineColor(kYellow+2);
   trigEff3->SetName("trigEff3");
   trigEff3->Draw();
   TLegend *legend3 = new TLegend(0.4,0.25,0.65,0.32);
   legend3->AddEntry(trigEff,"HSCP trigger - Only slope cut.","lep");
   legend3->SetTextSize(0.025);
   //legend3->SetTextColor();
   legend3->SetLineColor(1);
   legend3->SetBorderSize(0);
   legend3->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEff-HSCPtrigger-OnlySlopeCut.pdf");
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


	fHBeta_resolution->GetXaxis()->SetTitle("#beta resolution (#beta_{GEN}-#beta_{RPC})/#beta_{GEN}");
	fHBeta_resolution->SetLineColor(kGreen+2);
	fHBeta_resolution->SetFillStyle(3444);
	fHBeta_resolution->SetFillColor(kGreen+2);
	fHBeta_resolution->Draw();
	
	CMS_lumi( c1, iPeriod, iPos );
	c1->SaveAs("beta_GenRes.pdf");
  
  fHt0->GetXaxis()->SetTitle("Time [ns]");
  fHt0->SetLineColor(kAzure+1);
  fHt0->SetFillStyle(3444);
  //fHt0->SetFillColor(kAzure+2);
  fHt0->SetName("fHt0");
  //fHt0->Draw();
  fHbx->SetLineColor(kRed+1);
  fHbx->SetName("fHbx");
  fHbx->GetXaxis()->SetTitle("Time [ns]");
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
  c1->SetLogy();
  c1->SaveAs("timeZeroAndBx_log.pdf");
  
  c1->SetLogy(0);
  fHtime->GetXaxis()->SetTitle("Time [ns]");
  fHtime->SetLineColor(kPink+2);
  fHtime->SetFillStyle(3444);
  fHtime->SetFillColor(kPink+2);
  fHtime->Draw();
  CMS_lumi( c1, iPeriod, iPos );
  c1->SaveAs("time.pdf");
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
  /*
	fHBetaSel->Draw();
	fHBetaSel->SaveAs("histograms/"+infile+"_BetaSel.root");
	fHBetaSel->GetXaxis()->SetTitle("#beta(sTau)");
	c1->SaveAs("images/"+infile+"_BetaSel.pdf");
	
	c1->cd();
	gPad->SetRightMargin(0.05);
	gPad->SetLeftMargin(0.15);
	
	eff->SetTitle(";#beta;Efficiency");
	eff->Draw();
	eff->SaveAs("histograms/"+infile+"_eff.root");
	c1->SaveAs("images/"+infile+"_eff.pdf");
	
	
	
	
*/	
	return 0;
}
