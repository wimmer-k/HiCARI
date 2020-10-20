void showClover(int run){
  TFile *f = new TFile(Form("~/rootfiles/run%04d.root",run));
  TCanvas* cc = new TCanvas("cc","cc",800,800);
  cc->Divide(4,4);
  TCanvas* cs = new TCanvas("cs","cs",800,800);
  cs->Divide(4,4);
  TH1F *hc[4][4];
  TH2F *hs[4][4];
  for(int m=0;m<4;m++){
    for(int c=0;c<4;c++){
      cc->cd(m*4+c+1);
      hc[m][c] = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m+6,c));
      if(hc[m][c]!=NULL)
	hc[m][c]->Draw();
      cs->cd(m*4+c+1);
      gPad->SetLogz();
      hs[m][c] = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m+6,c));
      if(hs[m][c]!=NULL)
	hs[m][c]->Draw("colz");
    }
  }
}
void showMiniball(int run){
  TFile *f = new TFile(Form("~/rootfiles/run%04d.root",run));
  TCanvas* cc = new TCanvas("cc","cc",800,800);
  cc->Divide(3,6);
  TCanvas* cs = new TCanvas("cs","cs",800,800);
  cs->Divide(3,6);
  TH1F *hc[6][3];
  TH2F *hs[6][3];
  for(int m=0;m<6;m++){
    for(int c=0;c<3;c++){
      cc->cd(m*3+c+1);
      hc[m][c] = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
      if(hc[m][c]!=NULL)
	hc[m][c]->Draw();
      cs->cd(m*3+c+1);
      gPad->SetLogz();
      hs[m][c] = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
      if(hs[m][c]!=NULL)
	hs[m][c]->Draw("colz");
    }
  }
}
