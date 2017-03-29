#include "plotHistosRecHits.h"


int plotHistosRecHits()
{
	setStyleTemplate();
   writeExtraText = true;       // if extra text
  	extraText  = "Preliminary";  // default extra text is "Preliminary"
  	lumi_13TeV = "2.7 fb^{-1}";
  	lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  	lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  	
  	int iPeriod = 100;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  	int iPos    = 33;
  	
  	TString canvName = "c1";
  	TCanvas* c1 = CreateCanvas(canvName);
	
	TFile* fin = new TFile("rootfiles/"+infile+".root");

	TTree* treeRecHit = (TTree*)fin->Get(intreeRecHit);	
	treeRecHit->SetBranchAddress("bunchX",			&bunchX);
	treeRecHit->SetBranchAddress("stationhit",	&stationhit);

	Int_t nRecHits = treeRecHit->GetEntries();
	for(Int_t j = 0; j < nRecHits; j++)
	{
		treeRecHit->GetEntry(j);
		fHist2DBx_Station->Fill(stationhit,bunchX);
	}
	
	TTree* treeEff = (TTree*)fin->Get(intreeEff);
	treeEff->SetBranchAddress("isChosen", &isChosen);
	treeEff->SetBranchAddress("beta", &beta);
	treeEff->SetBranchAddress("betaSim", &betaSim);
	treeEff->SetBranchAddress("betaGen", &betaGen);
	treeEff->SetBranchAddress("etaGen", &etaGen);
	Int_t ntks = treeEff->GetEntries();
	for(Int_t k=0; k < ntks ; k++)
	{
		treeEff->GetEntry(k);
		int chosen= isChosen;
		double b = beta;
		double bSim = betaSim;
		double bGen = betaGen;
		if(b>0 && b <1.)
		{
			fHBetaSimRes->Fill((b-bSim)/bSim);
			fHBetaGenRes->Fill((b-bGen)/bGen);
			cout << bSim << "   " << b << endl;
			fHBetaTot->Fill(b);
			if(chosen) fHBetaSel->Fill(b);
			fHist2DBetaVsEta->Fill(betaGen,etaGen);
		}
		
	}
	
	TEfficiency* eff = 0;
	if(TEfficiency::CheckConsistency(*fHBetaSel,*fHBetaTot))
	{
  		eff = new TEfficiency(*fHBetaSel,*fHBetaTot);
 	}
	
   
	c1 = CreateCanvas(canvName);
	gStyle->SetOptStat(0);
   gPad->SetRightMargin(0.15);
   c1->cd();
	fHist2DBx_Station->Draw("colz");
	c1->SaveAs("images/"+infile+"_BxVsStation.pdf");
	
   gPad->SetRightMargin(0.15);
	c1->cd();
	fHist2DBetaVsEta->SaveAs("histograms/"+infile+"_BetaVsEta.root");
	fHist2DBetaVsEta->GetXaxis()->SetTitle("#beta(sTau)");
	fHist2DBetaVsEta->GetYaxis()->SetTitle("#eta(sTau)");
	fHist2DBetaVsEta->Draw("colz");
	c1->SaveAs("images/"+infile+"_BetaVsEta.pdf");
	
	
	gStyle->SetOptStat(1);
	c1->cd();
	fHBetaTot->Draw();
	fHBetaTot->SaveAs("histograms/"+infile+"_BetaTot.root");
	fHBetaTot->GetXaxis()->SetTitle("#beta(sTau)");
	c1->SaveAs("images/"+infile+"_BetaTot.pdf");
	
	c1->cd();
	fHBetaSel->Draw();
	fHBetaSel->SaveAs("histograms/"+infile+"_BetaSel.root");
	fHBetaSel->GetXaxis()->SetTitle("#beta(sTau)");
	c1->SaveAs("images/"+infile+"_BetaSel.pdf");
	
	c1->cd();
	//eff->GetYaxis()->SetTitle("Efficiency");
	eff->Draw();
	eff->SaveAs("histograms/"+infile+"_eff.root");
	c1->SaveAs("images/"+infile+"_eff.pdf");
	
	c1->cd();
	fHBetaSimRes->SaveAs("histograms/"+infile+"_SimRes.root");
	fHBetaSimRes->GetXaxis()->SetTitle("#beta resolution (#beta - #beta_{SIM})/#beta_{SIM}");
	fHBetaSimRes->Draw();
	c1->SaveAs("images/"+infile+"_SimRes.pdf");
	
	c1->cd();
	fHBetaGenRes->SaveAs("histograms/"+infile+"_GenRes.root");
	fHBetaGenRes->GetXaxis()->SetTitle("#beta resolution (#beta - #beta_{GEN})/#beta_{GEN}");
	fHBetaGenRes->Draw();
	c1->SaveAs("images/"+infile+"_GenRes.pdf");
	
	
	return 0;
}
