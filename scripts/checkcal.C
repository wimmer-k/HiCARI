TFile *f = new TFile("~/rootfiles/cal0394bkp.root");
TCanvas *c = NULL;
TH1F* h;
TH2F* h2;
void checkcal(int u, int y){
  if(c==NULL){
    gStyle->SetOptStat(0);
    c = new TCanvas("c","c",1900,800);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetLogy();
    c->cd(2);
    gPad->SetLogz();
  }
  h = NULL;
  h2 = NULL;
  h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",u,y));
  h2 = (TH2F*)f->Get(Form("h_segen_vs_nr_clus%02d_crys%02d",u,y));

  if(h!=NULL){
    c->cd(1);    
    h->Draw();
  }
  if(h2!=NULL){
    c->cd(2);
    h2->Draw("colz");
    h2->GetYaxis()->SetRangeUser(900,1500);
  }

}
