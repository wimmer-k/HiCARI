{
  TFile *f = new TFile("~/rootfiles/test.root");
  TCanvas* cc = new TCanvas("cc","cc",800,800);
  cc->Divide(2,2);
  TCanvas* cs = new TCanvas("cs","cs",800,800);
  cs->Divide(4,4);
  TH1F *hc[4];
  TH2F *hs[4][4];
  for(int m=0;m<4;m++){
    cc->cd(m+1);
    //gPad->SetLogz();
    hc[m] = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m+6));
    hc[m]->Draw();
    for(int c=0;c<4;c++){
      cs->cd(m*4+c+1);
      gPad->SetLogz();
      hs[m][c] = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m+6,c));
      hs[m][c]->Draw("colz");
    }
  }
  


}
