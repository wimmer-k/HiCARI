TFile *f = new TFile("~/rootfiles/cal0399.root");
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
void checkcores(int run){
  f = new TFile(Form("./hist/hcal%04d.root",run));
  h2 = (TH2F*)f->Get("h_en_summary");
  if(h2!=NULL){
    h2->Draw("colz");
    h2->GetYaxis()->SetRangeUser(0,1500);
  }
  
}
void plotall(){
  THStack *hs = new THStack("hs","hs");
  for(int r=448;r<483;r++){
    f = NULL;
    f = new TFile(Form("./hist/hcal%04d.root",r));
    if(f==NULL)
      continue;
    h2 = NULL;
    h2 = (TH2F*)f->Get("h_en_summary");
    if(h2==NULL)
      continue;
    h = (TH1F*)h2->ProjectionY(Form("hp_%04d",r));
    h->GetXaxis()->SetRangeUser(100,1500);
    h->Scale(1./h->Integral(1900,3900));
    h->SetLineColor(r%6+1);
    hs->Add(h);
  }
  hs->Draw("nostack");
}
