#include "plotHistosRecHits.h"


int plotHistosRecHits()
{
	setStyleTemplate();
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
   TH1D *fHBeta_pas = (TH1D*)fin->Get("demo2/fHbeta_pas");
   TH1D *fHBeta_tot = (TH1D*)fin->Get("demo2/fHbeta_tot");
   TH1D *fHBeta_pas_BetaErrorCut = (TH1D*)fin->Get("demo2/fHbeta_pas_BetaErrorCut");
   TH1D *fHBeta_pas_SlopeCut = (TH1D*)fin->Get("demo2/fHbeta_pas_SlopeCut");
   TH1D *fHBeta_resolution = (TH1D*)fin->Get("demo2/fHres");
   TH1D *fHt0 = (TH1D*)fin->Get("demo2/fHt0");
   TH1D *fHtime = (TH1D*)fin->Get("demo2/fHrpcHitsTime");
   TH1D *fHEta_pas = (TH1D*)fin->Get("demo2/fHeta_pas");
   TH1D *fHEta_tot = (TH1D*)fin->Get("demo2/fHeta_tot");
   TH1D *fHbx = (TH1D*)fin->Get("demo2/fHt0_bx");

   TH1D * fHBeta_MuTrig_pas = (TH1D*) fin2->Get("demo2/fHbeta_pas");
   TH1D * fHBeta_MuTrig_tot = (TH1D*) fin2->Get("demo2/fHbeta_tot");

   TH1D *fHBeta_Zmm_pas = (TH1D*)fin3->Get("demo2/fHbeta_pas");
   TH1D *fHBeta_Zmm_tot = (TH1D*)fin3->Get("demo2/fHbeta_tot");
   TH1D *fHtime_Zmm = (TH1D*)fin3->Get("demo2/fHrpcHitsTime");

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


	fHBeta_resolution->SaveAs("beta_GenRes.root");
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
  c1->SetLogy();

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
