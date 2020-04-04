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
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TMarker.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
char* fileEu = (char*)"/home/gamma20/rootfiles/run0455.root";
char* fileBG = (char*)"/home/gamma20/rootfiles/run0470.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t f2gammagaussbg(Double_t *x, Double_t *par);
Double_t f2gammabg(Double_t *x, Double_t *par);
Double_t f2gammastep(Double_t *x, Double_t *par);
Double_t f2gammagaus(Double_t *x, Double_t *par);
vector <vector<double> > fitEu(TH1F* h, bool draw);
double resolution = 2;
double frange = 15;
TCanvas *ca;
vector <vector<double> > defreturn;

void core(int clu, int cry){
  TFile *f = new TFile(fileEu);
  TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",clu,cry));
  if(h==NULL)
    return;
  fitEu(h,1);
}
void fullfile(char* filen = NULL){
  TFile *f;
  if(filen == NULL)
    f = new TFile(fileEu);
  else
    f = new TFile(filen);
  TH2F* h2 = (TH2F*)f->Get("h_en_summary");
  TH1F* h = (TH1F*)h2->ProjectionY("hp");
  if(h==NULL)
    return;
  fitEu(h,1);
}
void effEu(double runtime=0, double act=0){
  TFile *f = new TFile(fileEu);
  TH1F* h =NULL;
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  TCanvas* ca2 = new TCanvas("ca2","ca2",1200,800);
  TCanvas* ca3 = new TCanvas("ca3","ca3",1200,800);
  TMultiGraph *mg = new TMultiGraph();
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      h = NULL;
      h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	continue;
      vector<vector<double> > r = fitEu(h,0);
      if(runtime>0 && act >0){
	for(int i=0;i<r.at(0).size();i++){
	  r.at(6)[i] = r.at(6)[i]/runtime/act; 
	  r.at(7)[i] = r.at(7)[i]/runtime/act; 
	}
      }
      ca->cd(cry*12+1+clu);
      TGraphErrors *ge = new TGraphErrors(r.at(0).size(),&r.at(5)[0],&r.at(6)[0],0,&r.at(7)[0]);
      ge->Draw("AP*");
      mg->Add(ge,"P");

    }
  }
  ca2->cd();
  mg->Draw("a");
  h = NULL;
  TH2F* h2 = (TH2F*)f->Get("h_en_summary");
  h = (TH1F*)h2->ProjectionY("hp");
  vector<vector<double> > r = fitEu(h,0);
  if(runtime>0 && act >0){
    for(int i=0;i<r.at(0).size();i++){
      r.at(6)[i] = r.at(6)[i]/runtime/act; 
      r.at(7)[i] = r.at(7)[i]/runtime/act; 
    }
  }
  ca3->cd();
  TGraphErrors *ge = new TGraphErrors(r.at(0).size(),&r.at(5)[0],&r.at(6)[0],0,&r.at(7)[0]);
  ge->Draw("AP*");
  
}


vector <vector<double> > fitEu(TH1F* h, bool draw){
  defreturn.resize(4);
  defreturn[0].push_back(-1);
  if(h->Integral() < 10)
    return defreturn;
  TSpectrum *sp = new TSpectrum(2,resolution);
  sp->SetResolution(resolution);
  const int n = 10;
  ifstream intensity;
  intensity.open("/home/gamma20/HiCARI/scripts/Eudecay.dat");
  intensity.ignore(1000,'\n');
  double en[n], inten[n];//, eff[n];
  if(draw){
    ca = new TCanvas("ca","ca",600,600);
    ca->Divide(4,3);
  }
  vector< vector<double> > peaks;
  peaks.resize(8);
  for(int i=0;i<8;i++)
    peaks.at(i).clear();
  for(int p=0; p<n; p++){
    intensity >> en[p] >> inten[p];
    cout << en[p] << "\t" << inten[p] << endl;
    h->GetXaxis()->SetRangeUser(en[p]-frange,en[p]+frange);
    Int_t nfound = 0;
    if(draw){
      ca->cd(1+p);
      nfound = sp->Search(h,resolution,"nobackground",0.7);
    }
    else
      nfound = sp->Search(h,resolution,"nobackgroundgoff",0.7);
    if(draw){
      ca->cd(1+p);
      h->DrawCopy();
    }
    if(nfound!=1){
      cout << "Found " << nfound << " peaks in spectrum, not 1" << endl;
      continue;
    }
    Float_t *xpeaks = sp->GetPositionX();
    TF1 *fu;
    TF1 *fus[3];
    h->GetXaxis()->SetRangeUser(xpeaks[0]-frange,xpeaks[0]+frange);
    fu = new TF1(Form("f%s_p%d",h->GetName(),p),fgammagaussbg,xpeaks[0]-frange,xpeaks[0]+frange,6);
    fu->SetLineColor(3);
    fu->SetLineWidth(1);
    fu->SetParameter(0,0);//bg const
    fu->SetParameter(1,0);//bg slope
    // if(!bg){
    //   fu->FixParameter(0,0);//bg const
    //   fu->FixParameter(1,0);//bg slope
    // }
    fu->SetParameter(2,h->Integral(xpeaks[0]-frange,xpeaks[0]+frange));//norm
    fu->SetParameter(3,xpeaks[0]);//mean
    fu->SetParLimits(3,xpeaks[0]-50,xpeaks[0]+50);//mean
    fu->SetParameter(4,1);//sigma
    fu->SetParLimits(4,0.1,10);//sigma
    fu->SetParameter(5,h->GetBinContent(h->FindBin(xpeaks[0]-frange)));//step
    if(draw)
      h->Fit(fu,"Rn");
    else
      h->Fit(fu,"Rqn");

    //draw results.
    if(draw){
      ca->cd(1+p);
      h->DrawCopy();
      fu->Draw("same");
      fus[0] = new TF1(Form("f%s_p%d_bg",h->GetName(),p),fgammabg,xpeaks[0]-frange,xpeaks[0]+frange,6);
      fus[1] = new TF1(Form("f%s_p%d_st",h->GetName(),p),fgammastep,xpeaks[0]-frange,xpeaks[0]+frange,6);
      fus[2] = new TF1(Form("f%s_p%d_ga",h->GetName(),p),fgammagaus,xpeaks[0]-frange,xpeaks[0]+frange,6);
	

      fus[0]->SetLineColor(5);
      fus[1]->SetLineColor(4);
      fus[2]->SetLineColor(2);
      for(int k=0;k<3;k++){
	fus[k]->SetLineWidth(1);
	for(int l=0;l<6;l++)
	  fus[k]->SetParameter(l,fu->GetParameter(l));
	fus[k]->Draw("same");
      }
    }// draw
    //cout << (peaks.at(0)).size() << "\t" << fu->GetParameter(3) << endl;
    (peaks.at(0)).push_back(fu->GetParameter(3)-en[p]); // mean
    (peaks.at(1)).push_back(fu->GetParError(3));  // error mean
    (peaks.at(2)).push_back(fu->GetParameter(4)); // sigma
    (peaks.at(3)).push_back(fu->GetParameter(2)/h->GetBinWidth(1)); // content
    (peaks.at(4)).push_back(fu->GetParError(2)/h->GetBinWidth(1)); // content error 
    (peaks.at(5)).push_back(en[p]);
    (peaks.at(6)).push_back(fu->GetParameter(2)/h->GetBinWidth(1)/inten[p]);
    (peaks.at(7)).push_back(fu->GetParError(2)/h->GetBinWidth(1)/inten[p]);
  }// peaks
  if(draw){
    //cout << (peaks.at(0)).size() << endl;
    ca->cd(n+1);
    TGraphErrors *g = new TGraphErrors(peaks.at(0).size(),&peaks.at(5)[0],&peaks.at(0)[0],0,&peaks.at(1)[0]);
    g->Draw("AP*");
    //g->Fit("pol1");    
    ca->cd(n+2);
    TGraphErrors *ge = new TGraphErrors(peaks.at(0).size(),&peaks.at(5)[0],&peaks.at(6)[0],0,&peaks.at(7)[0]);
    ge->Draw("AP*");
  }
  return peaks;
}
void findall(bool AB = false, double runtime=0, double act=0, char* filen = NULL){
  TFile *f;
  if(filen == NULL)
    f = new TFile(fileEu);
  else
    f = new TFile(filen);
  
  TH2F* h2 = (TH2F*)f->Get("h_en_summary");
  if(AB)
    h2 = (TH2F*)f->Get("hAB_en_summary");
  TH1F* h = (TH1F*)h2->ProjectionY("hp");
  if(h==NULL)
    return;
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(1,2);
  ca->cd(1);
  h->GetXaxis()->SetRangeUser(10,3010);
  h->DrawCopy();
  TF1* fu = new TF1(Form("f%s",h->GetName()),fgammagaussbg,1450,1475,6);
  fu->SetLineColor(3);
  fu->SetLineWidth(1);
  fu->SetParameter(0,0);//bg const
  fu->SetParameter(1,0);//bg slope
  fu->SetParameter(2,h->Integral(1450,1470));//norm
  fu->SetParameter(3,1460);//mean
  fu->SetParLimits(3,1450,1470);//mean
  fu->SetParameter(4,1);//sigma
  fu->SetParLimits(4,0.1,10);//sigma
  fu->SetParameter(5,h->GetBinContent(h->FindBin(1450)));//step
  h->Fit(fu,"R");
  double p = fu->GetParameter(2);
  f = new TFile(fileBG);
  h2 = (TH2F*)f->Get("h_en_summary");
  if(AB)
    h2 = (TH2F*)f->Get("hAB_en_summary");
  TH1F* hbg = (TH1F*)h2->ProjectionY("hbgp");
  //hbg->DrawCopy("same");
  fu->SetParameter(2,hbg->Integral(1450,1470));//norm
  fu->SetParameter(5,hbg->GetBinContent(h->FindBin(1450)));//step
  hbg->Fit(fu,"Rn");
  fu->Draw("same");
  p /= fu->GetParameter(2);
  TH1F* hbgs = (TH1F*)hbg->Clone("hbgs");
  hbgs->Scale(p);
  hbgs->SetLineColor(2);
  hbgs->DrawCopy("same");

  ca->cd(2);
  TH1F* hs = (TH1F*)h->Clone("hs");
  hs->Add(hbgs,-1);
  hs->GetXaxis()->SetRangeUser(10,3010);
  hs->DrawCopy();
  ifstream intensity;
  intensity.open("/home/gamma20/HiCARI/scripts/Eudecaydetail.dat");
  intensity.ignore(1000,'\n');
  vector<double> en;
  vector<double> in;
  vector<double> fc;
  vector<double> fe;
  TCanvas* ca2 = new TCanvas("ca2","ca2",1200,800);
  ca2->Divide(5,3);
  //                      1 2 3    4    5 6   7 8   9 10    11    12 13    14
  double secpeak[14] = {351,0,0,1085,1292,0,238,0,862, 0, 1091, 1120, 0, 1402}; 
  while(!intensity.eof()){
    double e,i;
    intensity >> e >> i;
    
    if(intensity.eof())
      break;
    intensity.ignore(1000,'\n');
    if(e>100 && i > 1 && i*125000./28.53 > h->GetBinContent(h->FindBin(e))){
      ca->cd(2);
      TLine *l = new TLine(e,0,e,h->GetMaximum());
      l->SetLineColor(2);
      l->Draw();
      TMarker *m = new TMarker(e,i*110000/28.53,20);
      m->SetMarkerColor(2);
      m->Draw();
      cout << e << "\t" << i << endl;
      en.push_back(e);
      in.push_back(i);
      ca2->cd(en.size());
      TF1 *fus[3];
      h->GetXaxis()->SetRangeUser(e-frange,e+frange);
      h->GetYaxis()->UnZoom();
      frange = 12;
      //one peak
      int pars = 6;
      if(secpeak[en.size()-1]>0){
	pars = 9;
	fu = new TF1(Form("f%s_p%d",h->GetName(),(int)e),f2gammagaussbg,e-frange,e+frange,pars);
	fus[0] = new TF1(Form("f%s_p%d_bg",h->GetName(),(int)e),f2gammabg,e-frange,e+frange,pars);
	fus[1] = new TF1(Form("f%s_p%d_st",h->GetName(),(int)e),f2gammastep,e-frange,e+frange,pars);
	fus[2] = new TF1(Form("f%s_p%d_ga",h->GetName(),(int)e),f2gammagaus,e-frange,e+frange,pars);
      }
      else{
	fu = new TF1(Form("f%s_p%d",h->GetName(),(int)e),fgammagaussbg,e-frange,e+frange,pars);
	fus[0] = new TF1(Form("f%s_p%d_bg",h->GetName(),(int)e),fgammabg,e-frange,e+frange,pars);
	fus[1] = new TF1(Form("f%s_p%d_st",h->GetName(),(int)e),fgammastep,e-frange,e+frange,pars);
	fus[2] = new TF1(Form("f%s_p%d_ga",h->GetName(),(int)e),fgammagaus,e-frange,e+frange,pars);
      }
      fu->SetLineColor(3);
      fu->SetLineWidth(1);
      fu->SetParameter(0,0);//bg const
      fu->SetParameter(1,0);//bg slope
      fu->SetParameter(2,h->Integral(e-frange,e+frange));//norm
      fu->SetParameter(3,e);//mean
      fu->SetParLimits(3,e-5,e+5);//mean
      fu->SetParameter(4,1);//sigma
      fu->SetParLimits(4,0.1,10);//sigma
      if(pars==6){
	fu->SetParameter(5,h->GetBinContent(h->FindBin(e-frange))/2);//step
	fu->SetParLimits(5,0.0,h->GetBinContent(h->FindBin(e-frange)));//step
      }
      if(pars==9){
	fu->SetParameter(8,h->GetBinContent(h->FindBin(e-frange))/4);//step
	fu->SetParLimits(8,0.0,h->GetBinContent(h->FindBin(e-frange))/2);//step
	fu->SetParameter(5,h->Integral(e-frange,e+frange)/20);//norm
	fu->SetParameter(6,secpeak[en.size()-1]);//mean
	fu->SetParLimits(6,secpeak[en.size()-1]-1,secpeak[en.size()-1]+1);//mean
	fu->SetParameter(7,1);//sigma
	fu->SetParLimits(7,0.8,1.5);//sigma
      }
      h->Fit(fu,"Rn");
      h->DrawCopy();
      fu->Draw("same");
	

      fus[0]->SetLineColor(5);
      fus[1]->SetLineColor(4);
      fus[2]->SetLineColor(2);
      for(int k=0;k<3;k++){
	fus[k]->SetLineWidth(1);
	for(int j=0;j<pars;j++)
	  fus[k]->SetParameter(j,fu->GetParameter(j));
	fus[k]->Draw("same");
      }
      if(runtime>0 && act >0){
	fc.push_back(fu->GetParameter(2)/h->GetBinWidth(1)/(i/100)/runtime/act);
	fe.push_back(fu->GetParError(2)/h->GetBinWidth(1)/(i/100)/runtime/act);
      }
      else{
	fc.push_back(fu->GetParameter(2)/h->GetBinWidth(1)/i);
	fe.push_back(fu->GetParError(2)/h->GetBinWidth(1)/i);
      }
    }//good energies
  }
  cout << en.size() << endl;
  TGraphErrors* g = new TGraphErrors(en.size(),&en[0],&fc[0],0,&fe[0]);
  ca2->cd(en.size()+1);
  g->Draw("AP");
}
// one peak fitting function
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
// two peak fitting function
Double_t f2gammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];

  Double_t step = par[8];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t f2gammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];

  Double_t step = par[8];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  double result = 0;
  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);

  return result;

}
