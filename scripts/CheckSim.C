#include <fstream>
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "fit.h"

TH1F* loaddata(int run){
  TFile *f = new TFile(Form("hist/hcal%04d.root",run));
  TH2F* h2d = (TH2F*)f->Get("h_en_summary_fine");
  TH1F* h = (TH1F*)h2d->ProjectionY(Form("h%04d",run));
  return h;
}
TH1F* loadsim(int source, double time, double act){
  TFile *f;
  if(source == 0)
    f = new TFile("simulation/hist/mysourceCo60.root");
  else
    f = new TFile("simulation/hist/mysourceEu152.root");
  TH2F* h2d = (TH2F*)f->Get("egam_summary_fine");
  TH1F* h = (TH1F*)h2d->ProjectionY(Form("hs%d",source));
  if(time>1 && act >1)
    h->Scale(time*act/1e6);
  if(source == 0)
    h->Scale(10);
  return h;

}
void CheckSim(){
  TCanvas* c = new TCanvas("c","c",900,400);
  c->Divide(2,1);

  TH1F* bg[2];
  bg[0] = loaddata(470);
  //bg[0]->Scale(1./5294.5);
  bg[1] = (TH1F*)bg[0]->Clone("bgcopy");

  TH1F* d[2];
  d[0] = loaddata(475);
  d[1] = loaddata(455);

  double t[2] = { 3735.0,  3629.2};
  double a[2] = {23484.2,  6349.4};

  TH1F* s[2];
  TH1F* r[2];
  double p[2][2];
  TF1* fu[2][2];
  for(int i=0;i<2;i++){
    cout <<"------------------------------"<<endl;
    c->cd(1+i);
    d[i]->Draw();
    d[i]->GetXaxis()->SetRangeUser(10,1550);
    fu[i][0] = new TF1(Form("f%s",d[i]->GetName()),f2gammagaussbg,1450,1475,9);
    fu[i][0]->SetLineWidth(1);
    fu[i][0]->SetParameter(0,0);//bg const
    fu[i][0]->SetParameter(1,0);//bg slope
    fu[i][0]->SetParameter(2,d[i]->Integral(1450,1470)/2);//norm
    fu[i][0]->SetParameter(3,1460);//mean
    fu[i][0]->SetParLimits(3,1450,1470);//mean
    fu[i][0]->SetParameter(4,1);//sigma
    fu[i][0]->SetParLimits(4,0.1,3);//sigma
    fu[i][0]->SetParameter(5,d[i]->Integral(1450,1470)/2);//norm
    fu[i][0]->SetParameter(6,1460);//mean
    fu[i][0]->SetParLimits(6,1450,1470);//mean
    fu[i][0]->SetParameter(7,4);//sigma
    fu[i][0]->SetParLimits(7,2,7);//sigma
    fu[i][0]->SetParameter(8,d[i]->GetBinContent(d[i]->FindBin(1450)));//step
    d[i]->Fit(fu[i][0],"Rn");
    fu[i][0]->SetLineColor(kGray);
    fu[i][0]->Draw("same");
    p[i][0] = fu[i][0]->GetParameter(2)+ fu[i][0]->GetParameter(5);
    
    fu[i][1] = new TF1(Form("f%s",bg[i]->GetName()),f2gammagaussbg,1450,1475,9);
    fu[i][1]->SetLineWidth(1);
    fu[i][1]->SetParameter(0,0);//bg const
    fu[i][1]->SetParameter(1,0);//bg slope
    fu[i][1]->SetParameter(2,bg[i]->Integral(1450,1470)/2);//norm
    fu[i][1]->SetParameter(3,1460);//mean
    fu[i][1]->SetParLimits(3,1450,1470);//mean
    fu[i][1]->SetParameter(4,1);//sigma
    fu[i][1]->SetParLimits(4,0.1,3);//sigma
    fu[i][1]->SetParameter(5,bg[i]->Integral(1450,1470)/2);//norm
    fu[i][1]->SetParameter(6,1460);//mean
    fu[i][1]->SetParLimits(6,1450,1470);//mean
    fu[i][1]->SetParameter(7,4);//sigma
    fu[i][1]->SetParLimits(7,2,7);//sigma
    fu[i][1]->SetParameter(8,bg[i]->GetBinContent(bg[i]->FindBin(1450)));//step
    bg[i]->Fit(fu[i][1],"Rn");
    fu[i][1]->SetLineColor(kRed-8);
    fu[i][1]->Draw("same");
    p[i][1] = fu[i][1]->GetParameter(2)+ fu[i][1]->GetParameter(5);
    cout << p[i][0] << "\t" << p[i][1] << endl;
    bg[i]->Scale(p[i][0]/p[i][1]);
    bg[i]->SetLineColor(2);
    bg[i]->Draw("same");
    
    //s[i] = loadsim(i,t[i],a[i]);
    s[i] = loadsim(i,1,1);
    s[i]->SetLineColor(4);
    s[i]->Draw("same");

    r[i] = (TH1F*)s[i]->Clone(Form("r%d",i));
    r[i]->Add(bg[i]);
    r[i]->SetLineColor(3);
    r[i]->Draw("same");
    
    gPad->SetLogy();
    cout << d[i]->GetName() << endl;
    ofstream fout;
    fout.open(Form("python/data/%s.dat",d[i]->GetName()));
    fout << d[i]->GetNbinsX() << "\t" << d[i]->GetBinLowEdge(1)<< "\t" << d[i]->GetBinLowEdge(d[i]->GetNbinsX()+1) << endl;
    for(int b=1;b<d[i]->GetNbinsX();b++)
      fout << d[i]->GetBinCenter(b) << "\t" << d[i]->GetBinContent(b) << "\t" << s[i]->GetBinContent(b) << "\t" << r[i]->GetBinContent(b) << "\t" << r[i]->GetBinContent(b) << "\t" << d[i]->GetBinContent(b)-s[i]->GetBinContent(b) << endl;
    fout.close();
  }
}
