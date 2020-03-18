#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TEnv.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

using namespace std;
vector<double> defreturn;
const int banks = 20;
const int slots = 12;
const int chans = 10;
const int corec = 9;
double resolution = 2;
double frange = 5000;
double range[2] = {150000,450000};
TCanvas *ca;
char* filen = (char*)"/home/gamma20/rootfiles/hist0268.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
vector<double> fit(TH1F* h, bool bg = true, bool draw = false);

void SetRange(double low, double hig){
  range[0] = low;
  range[1] = hig;
}
void test(int b=11, int s = 3, int c =9){
  TFile *f = new TFile(filen);
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_bank%02d_slot%02d_chan%02d",b,s,c));
  fit(h,1,1);
}
void core(int m=0, int c =0){
  TFile *f = new TFile(filen);
  TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  if(h==NULL)
    return;
  fit(h,1,1);
}
void segment(int m=0, int c =0, int s =0){
  resolution = 5;
  TFile *f = new TFile(filen);
  TH2F* h2 = (TH2F*)f->Get(Form("h_segen_vs_nr_clus%02d_crys%02d",m,c));
  TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
  if(h==NULL || h->Integral()<100)
    return;
  fit(h,0,1);
}
vector<double> fit(TH1F* h, bool bg, bool draw){
  defreturn.push_back(-1);
  h->GetXaxis()->SetRangeUser(range[0],range[1]);
  if(h->Integral() < 10)
    return defreturn;
  TSpectrum *sp = new TSpectrum(3,resolution);
  sp->SetResolution(resolution);
  Int_t nfound = 0;
  if(draw)
    nfound = sp->Search(h,resolution,"nobackground",0.4);
  else
    nfound = sp->Search(h,resolution,"nobackgroundgoff",0.4);
  if(nfound!=2){
    cout << "Found " << nfound << " peaks in spectrum, not 2, try again" << endl;
    TH1F* hc = (TH1F*)h->Clone("testclone");
    hc->Rebin(2);
    hc->GetXaxis()->SetRangeUser(range[0],range[1]);
    nfound = 0;
    if(draw)
      nfound = sp->Search(hc,resolution,"nobackground",0.4);
    else
      nfound = sp->Search(hc,resolution,"nobackgroundgoff",0.4);
      
    if(nfound!=2){
      cout << "Found " << nfound << " peaks in spectrum, not 2, aborting" << endl;
      if(draw)
	hc->DrawCopy();
      return defreturn;
    }      
  }
  Float_t *xpeaks = sp->GetPositionX();
  Float_t *ypeaks = sp->GetPositionY();

  if(draw){
    for(int p=0;p<nfound;p++){
      cout << xpeaks[p] << "\t" << ypeaks[p] << endl;
    }
    ca = new TCanvas("ca","ca",600,600);
    ca->Divide(2,2);
    ca->cd(1);
    h->DrawCopy();
  }
  //check if first peak is lower in energy, otherwise swap them
  if(xpeaks[0]>xpeaks[1]){
    Float_t temp = xpeaks[1];
    xpeaks[1] = xpeaks[0];
    xpeaks[0] = temp;
    temp = ypeaks[1];
    ypeaks[1] = ypeaks[0];
    ypeaks[0] = temp;
	
  }
  //h->Rebin(4);
  TF1 *fu[2];
  TF1 *fus[2][3];
  for(int p=0;p<nfound;p++){
    h->GetXaxis()->SetRangeUser(xpeaks[p]-frange,xpeaks[p]+frange);
    fu[p] = new TF1(Form("f%s_p%d",h->GetName(),p),fgammagaussbg,xpeaks[p]-frange,xpeaks[p]+frange,6);
    fu[p]->SetLineColor(3);
    fu[p]->SetLineWidth(1);
    fu[p]->SetParameter(0,0);//bg const
    fu[p]->SetParameter(1,0);//bg slope
    if(!bg){
      fu[p]->FixParameter(0,0);//bg const
      fu[p]->FixParameter(1,0);//bg slope
    }
    fu[p]->SetParameter(2,h->Integral(xpeaks[p]-frange,xpeaks[p]+frange));//norm
    fu[p]->SetParameter(3,xpeaks[p]);//mean
    fu[p]->SetParLimits(3,xpeaks[p]-500,xpeaks[p]+500);//mean
    fu[p]->SetParameter(4,200);//sigma
    fu[p]->SetParLimits(4,100,1000);//sigma
    fu[p]->SetParameter(5,h->GetBinContent(h->FindBin(xpeaks[p]-frange)));//step
    if(draw)
      h->Fit(fu[p],"Rn");
    else
      h->Fit(fu[p],"Rqn");

    //draw results.
    if(draw){
      ca->cd(2+p);
      h->DrawCopy();
      fu[p]->Draw("same");
      fus[p][0] = new TF1(Form("f%s_p%d_bg",h->GetName(),p),fgammabg,xpeaks[p]-frange,xpeaks[p]+frange,6);
      fus[p][1] = new TF1(Form("f%s_p%d_st",h->GetName(),p),fgammastep,xpeaks[p]-frange,xpeaks[p]+frange,6);
      fus[p][2] = new TF1(Form("f%s_p%d_ga",h->GetName(),p),fgammagaus,xpeaks[p]-frange,xpeaks[p]+frange,6);
	

      fus[p][0]->SetLineColor(5);
      fus[p][1]->SetLineColor(4);
      fus[p][2]->SetLineColor(2);
      for(int k=0;k<3;k++){
	fus[p][k]->SetLineWidth(1);
	for(int l=0;l<6;l++)
	  fus[p][k]->SetParameter(l,fu[p]->GetParameter(l));
	fus[p][k]->Draw("same");
      }

    }// draw
  }// peaks
  if(draw){
    ca->cd(nfound+2);
    double y[2] = {1173.2,1332.5};
    double x[2] = {fu[0]->GetParameter(3), fu[1]->GetParameter(3)};
    double e[2] = {fu[0]->GetParError(3), fu[1]->GetParError(3)};
    TGraphErrors *g = new TGraphErrors(2,x,y,e);
    g->Draw("AP");
    g->Fit("pol1");
  }
  
  //cout << fu[0]->GetParameter(3) << "\t" << fu[1]->GetParameter(3) << endl;
  double a =  (1332.5-1173.2)/(fu[1]->GetParameter(3)-fu[0]->GetParameter(3));
  double b = (1332.5+1173.2 -a*(fu[1]->GetParameter(3)+fu[0]->GetParameter(3)))/2;
  cout << a << "\t" << b << "\t" <<  a*fu[0]->GetParameter(4) << "\t" << a*fu[1]->GetParameter(4) <<"\t" << h->GetName() << endl;
  vector<double> rv;
  rv.push_back(a);
  rv.push_back(b);
  rv.push_back(a*fu[0]->GetParameter(4));
  rv.push_back(a*fu[1]->GetParameter(4));
  return rv;
}
void CalibGe(){
  frange = 3000;
  TEnv *cf = new TEnv(Form("%s.cal",filen));
  TFile *f = new TFile(filen);
  vector <double> gain;
  vector <double> offs;
  vector <double> res0;
  vector <double> res1;
  
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      //cout << clu << "\t" << cry << endl;
      TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	continue;
      vector<double> r = fit(h,1);
      //cout << " fitted " << r[0]<< endl;
      if(r[0]<0)
	continue;
      gain.push_back(r[0]);
      offs.push_back(r[1]);
      res0.push_back(r[2]);
      res1.push_back(r[3]);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),r[0]);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),r[1]);
    }
  }
  cout << "finished" << endl;
  ca = new TCanvas("ca","ca",900,400);
  ca->Divide(2,2);
  ca->cd(1);
  //cout << gain.size() << endl;
  TGraph* g = new TGraph(gain.size(),&gain[0],&offs[0]);
  g->Draw("AP*");
  
  ca->cd(2);
  TGraph* gr = new TGraph(res0.size(),&res0[0],&res1[0]);
  gr->Draw("AP*");
  
  gain.clear();
  offs.clear();
  res0.clear();
  res1.clear();
  resolution = 2;
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      //cout << clu << "\t" << cry << endl;
      TH2F* h2 = (TH2F*)f->Get(Form("h_segen_vs_nr_clus%02d_crys%02d",clu,cry));
      if(h2==NULL)
	continue;
      for(int s=0;s<40;s++){
	//cout << s << endl;
	TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
	if(h==NULL || h->Integral()<100)
	  continue;
	frange = 3000;
	SetRange(150000,450000);
	if(clu==2 && cry==1 && s==0)
	  SetRange(250000,350000);
	if(clu==3 && cry==2 && s==4)
	  SetRange(280000,380000);
	if(clu==4 && cry==1 && s==2)
	  SetRange(280000,380000);
	if(clu==1 && cry==1 && s==2)
       	  frange = 2000;
	if(clu==3 && cry==0 && s==3)
       	  frange = 2000;
	vector<double> r = fit(h,0);
	//cout << " fitted " << r[0]<< endl;
	if(r[0]<0)
	  continue;
	gain.push_back(r[0]);
	offs.push_back(r[1]);
	res0.push_back(r[2]);
	res1.push_back(r[3]);
	//if(r[2]>3)
	//  return;
	//	if(r[2]<0.8)
	//	  return;
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Gain",clu,cry,s),r[0]);
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Offset",clu,cry,s),r[1]);
      }//segments
    }//crystals
  }//clusters
  cout << "finished segments" << endl;
  ca->cd(3);
  //cout << gain.size() << endl;
  TGraph* gs = new TGraph(gain.size(),&gain[0],&offs[0]);
  gs->Draw("AP*");
  
  ca->cd(4);
  TGraph* grs = new TGraph(res0.size(),&res0[0],&res1[0]);
  grs->Draw("AP*");
  
  cf->SaveLevel(kEnvLocal);
}
void CalibMode3(){
  TFile *f = new TFile(filen);
  vector <double> gain;
  vector <double> offs;
  vector <double> res0;
  vector <double> res1;
  
  for(int b=0;b<banks;b++){
    for(int s=0;s<slots;s++){
      for(int c=0;c<chans;c++){
	TH1F* h = (TH1F*)f->Get(Form("hraw_en_bank%02d_slot%02d_chan%02d",b,s,c));
	if(c>5 && c<9)
	  continue;
	if(!h)
	  continue;
	if(c==9)
	  SetRange(100000,500000);
	else if( (s==0 && c==3) || (s==2 && c==2))
	  SetRange(280000,500000);
	else
	  SetRange(250000,500000);
	//cout << h->GetName() << "\t"<< h->GetEntries() << endl;
	vector<double> r = fit(h);
	//cout << " fitted " << r[0]<< endl;
	if(r[0]<0)
	  continue;
	gain.push_back(r[0]);
	offs.push_back(r[1]);
	res0.push_back(r[2]);
	res1.push_back(r[3]);
	
      }
    }
  }
  cout << "finished" << endl;
  ca = new TCanvas("c","c",900,400);
  ca->Divide(2,1);
  ca->cd(1);
  //cout << gain.size() << endl;
  TGraph* g = new TGraph(gain.size(),&gain[0],&offs[0]);
  g->Draw("AP*");
  
  ca->cd(2);
  TGraph* gr = new TGraph(res0.size(),&res0[0],&res1[0]);
  gr->Draw("AP*");
  
  
}

Double_t fgammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t fgammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  //  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];
  arg = (x[0]-mean)/(sqrt2*sigma);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];


  arg = (x[0]-mean)/(sqrt2*sigma);
  Double_t result = 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  return result;

}
