int plot_eff()
{
   TString infile = "m1599_eff.root";
   TFile* fin = new TFile(infile);
   TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
   TH1D *fHBeta_pas = (TH1D*)fin->Get("demo2/fHbeta_pas");
   TH1D *fHBeta_tot = (TH1D*)fin->Get("demo2/fHbeta_tot");
   TEfficiency* trigEff = new TEfficiency(*fHBeta_pas,*fHBeta_tot);
   trigEff->Draw();
   c1->SaveAs("trigEff.pdf");
   return 0;

}
