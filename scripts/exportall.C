void exportall(){
  TFile *f = NULL;
  TH2F* h2;
  TH1F* h;
  for(int r=447; r<483;r++){
    f = new TFile(Form("hist/hallcal%04d.root",r));
    if(!f->IsOpen()){
      continue;
    }
    h2 = (TH2F*)f->Get("h_en_summary");
    h = (TH1F*)h2->ProjectionY(Form("h_%04d",r));
    hist1D(h,"python/data",Form("egamall_%04d",r));
    f->Close();
  }
}
