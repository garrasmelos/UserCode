#include "CMS_lumi.cc"
#include "styleTemplate.cc"
int plot_eff()
{
	setStyleTemplate();
   writeExtraText = true;       // if extra text
  	extraText  = "Preliminary";  // default extra text is "Preliminary"
  	lumi_13TeV = "2.7 fb^{-1}";
  	lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  	lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  	lumi_13TeV_MC = "MC Pythia6 - #tilde{#tau} of GMSB model" ;
  	
  	int iPeriod = 100;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  	int iPos    = 11;
  	
   TString infile = "m1599_eff.root";
   TFile* fin = new TFile(infile);
   TCanvas * c1 = new TCanvas("c1", "c1", 1000, 1000);
   TH1D *fHBeta_pas = (TH1D*)fin->Get("demo2/fHbeta_pas");
   TH1D *fHBeta_tot = (TH1D*)fin->Get("demo2/fHbeta_tot");
   TEfficiency* trigEff = new TEfficiency(*fHBeta_pas,*fHBeta_tot);
   trigEff->SetTitle(";#beta_{GEN};Efficiency");
   trigEff->SetMarkerColor(kAzure+2);
   trigEff->SetLineColor(kAzure+2);
   trigEff->SetName("trigEff");
   trigEff->Draw();
   TLegend *legend = new TLegend(0.2,0.75,0.45,0.82);
   legend->AddEntry(trigEff,"HLT_L1SingleMuOpen_v3","lep");
   legend->SetTextSize(0.025);
   //legend->SetTextColor();
   legend->SetLineColor(1);
   legend->SetBorderSize(0);
   legend->Draw();
   CMS_lumi( c1, iPeriod, iPos );
   c1->SaveAs("trigEff.pdf");
   return 0;

}
