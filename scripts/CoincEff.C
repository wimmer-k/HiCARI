#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
char* fileCo = (char*)"./hist/hcal0469.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t f2gammagaussbg(Double_t *x, Double_t *par);
Double_t f2gammabg(Double_t *x, Double_t *par);
Double_t f2gammastep(Double_t *x, Double_t *par);
Double_t f2gammabgstep(Double_t *x, Double_t *par);
Double_t f2gammagaus0(Double_t *x, Double_t *par);
Double_t f2gammagaus1(Double_t *x, Double_t *par);
vector<double> notinthis(int clu, int cry, bool AB = false, bool draw = false);
void eff(int clu, int cry, bool AB =false, bool draw = false);
//double xpeaks[2] = {1173,1332};
double xpeaks[2] = {898,1831};
double grange = 8;
double crange = 7;
double frange[2] = {20,13};
TCanvas *ca;
ofstream fout;

using namespace std;
void hist1D(TH1F* h, char* foldername = "python", char* histname = NULL){
  ofstream fout;
  if(histname==NULL)
    fout.open(Form("%s/%s.dat",foldername,h->GetName()));
  else
    fout.open(Form("%s/%s.dat",foldername,histname));
  fout << h->GetNbinsX() << "\t" << h->GetBinLowEdge(1)<< "\t" << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;
  for(int i=1;i<h->GetNbinsX();i++)
    fout << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << endl;
  fout.close();
}
void func(TF1* f, double from, double to, int steps, char* foldername = "python", char* fname = NULL){
  ofstream fout;
  if(fname==NULL)
    fout.open(Form("%s/%s.dat",foldername,f->GetName()));
  else
    fout.open(Form("%s/%s.dat",foldername,fname));
  fout << steps << endl;
  for(int i=0;i<steps;i++){
    double x = from+i*(to-from)/steps;
    double y= f->Eval(x);
    //cout << i << "\t" << x <<"\t" << y << endl;
    fout << i << "\t" << x <<"\t" << y << endl;

  }
  fout.close();
}

void setCo(){
  xpeaks[0] = 1173;
  xpeaks[1] = 1332;
  fileCo = (char*)"./hist/hcal0475.root";
  crange = 4;
}
void eff(){
  fout.open("python/data/efficienciesY.dat");
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      eff(clu, cry);
    }
  }
  fout.close();
  fout.open("python/data/efficienciesABY.dat");
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      eff(clu, cry, 1);
    }
  }
  fout.close();
}

void eff(int clu, int cry, bool AB, bool draw){
  TFile* f = new TFile(fileCo);
  TH2F* h2;
  if(AB)
    h2 = (TH2F*)f->Get(Form("hAB_engg_TC_clus%02d_crys%02d",clu,cry));
  else
    h2 = (TH2F*)f->Get(Form("h_engg_TC_clus%02d_crys%02d",clu,cry));
  if(h2==NULL)
    return;
  h2->GetXaxis()->SetRangeUser(500,2000);
  h2->GetYaxis()->SetRangeUser(500,2000);
  if(h2->Integral()<1)
    return;
  
  TH1F* h[2];
  if(draw){
    ca = new TCanvas("ca","ca",1200,400);
    ca->Divide(3,1);
    ca->cd(1);
  }
  double c[2];
  double d[2];
  double e[2];
  for(int p=0;p<2;p++){
    h[p] = (TH1F*)h2->ProjectionX(Form("gated_p%d_cry%02d_clu%02d",p,clu,cry),h2->GetYaxis()->FindBin(xpeaks[p]-grange),h2->GetYaxis()->FindBin(xpeaks[p]+grange));
    c[p] = h[p]->Integral(h[p]->FindBin(xpeaks[1-p]-crange), h[p]->FindBin(xpeaks[1-p]+crange));
    e[p] = sqrt(c[p]);
    //d[p] = h[p]->Integral(h[p]->FindBin(xpeaks[p]-grange), h[p]->FindBin(xpeaks[p]+grange)); //fake coincidences
    //e[p] = sqrt(c[p] + d[p]);
    //c[p] -= d[p];
    if(draw){
      h[p]->Draw();
      ca->cd(p+2);
    }
  }
  vector<double> gates = notinthis(clu, cry, AB, draw);
  double ef[2];
  double ee[2];
  fout << clu << "\t" << cry;
  for(int p=0;p<2;p++){
    ef[p] = c[p]/gates.at(p*2);
    ee[p] = ef[p]*sqrt(pow(gates.at(p*2+1)/gates.at(p*2),2) + pow(e[p]/c[p],2));
    cout << xpeaks[p] << "\t" << ef[p] << " +- " << ee[p] << endl;
    fout << "\t" << ef[p] << "\t" << ee[p];
  }
  fout << endl;
}

vector<double> notinthis(int clu, int cry, bool AB, bool draw){
  TFile* f = new TFile(fileCo);
  TH2F* h2;
  if(AB)
    h2 = (TH2F*)f->Get("hAB_en_summary_fine");
  else
    h2 = (TH2F*)f->Get("h_en_summary_fine");
  //h2->ProjectionY()->Draw();
  TH1F* h;
  if(clu==0 && cry == 0)
    h = (TH1F*)h2->ProjectionY(Form("notin_cry%02d_clu%02d",clu,cry),2,48);
  else{
    h = (TH1F*)h2->ProjectionY(Form("notin_cry%02d_clu%02d",clu,cry),1,clu*4+cry);
    h->Add((TH1F*)h2->ProjectionY(Form("notin_cry%02d_clu%02d_part2",clu,cry),clu*4+cry+2,48));
  }
  // h->SetLineColor(2);
  // h->Draw("same");
  // TH1F* ht = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",clu,cry));
  // ht->SetLineColor(3);
  // ht->Draw("same");
  // TH1F* hc = (TH1F*)h->Clone(Form("notin_cry%02d_clu%02d_clone",clu,cry));
  // hc->Add(ht,1);
  // hc->SetLineColor(4);
  // hc->Draw("same");
  TF1 *fu[2];
  TF1 *fus[2][5];
  if(draw){
    h->Draw();
  }
  vector<double> rv;
  rv.clear();
    
  for(int p=0;p<2;p++){
    h->GetXaxis()->SetRangeUser(xpeaks[p]-frange[0],xpeaks[p]+frange[1]);
    fu[p] = new TF1(Form("f%s_p%d",h->GetName(),p),f2gammagaussbg,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
    fu[p]->SetLineColor(3);
    fu[p]->SetLineWidth(1);
    fu[p]->SetParameter(0,0);//bg const
    fu[p]->SetParameter(1,0);//bg slope
    fu[p]->SetParameter(2,h->Integral(xpeaks[p]-frange[0],xpeaks[p]+frange[1])/2);//norm
    fu[p]->SetParameter(3,xpeaks[p]);//mean
    fu[p]->SetParLimits(3,xpeaks[p]-50,xpeaks[p]+50);//mean
    fu[p]->SetParameter(4,1);//sigma
    fu[p]->SetParLimits(4,0.1,2);//sigma
    fu[p]->SetParameter(5,h->Integral(xpeaks[p]-frange[0],xpeaks[p]+frange[1])/2);//norm
    fu[p]->SetParameter(6,xpeaks[p]);//mean
    fu[p]->SetParLimits(6,xpeaks[p]-50,xpeaks[p]+50);//mean
    fu[p]->SetParameter(7,5);//sigma
    fu[p]->SetParLimits(7,2,10);//sigma
    fu[p]->SetParameter(8,h->GetBinContent(h->FindBin(xpeaks[p]-frange[0])));//step
    if(draw)
      h->Fit(fu[p],"Rn");
    else
      h->Fit(fu[p],"Rqn");

    rv.push_back(fu[p]->GetParameter(2)/h->GetBinWidth(1) + fu[p]->GetParameter(5)/h->GetBinWidth(1));
    rv.push_back(sqrt(fu[p]->GetParError(2)*fu[p]->GetParError(2)+fu[p]->GetParError(5)*fu[p]->GetParError(5))/h->GetBinWidth(1));
    //draw results.
    if(draw){
      hist1D(h,"python/data");
      fu[p]->Draw("same");
      if(p==0)
	func(fu[p],xpeaks[p]-frange[0],xpeaks[p]+frange[1],1000,"python/data");
      fus[p][0] = new TF1(Form("f%s_p%d_bg",h->GetName(),p),f2gammabg,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
      fus[p][1] = new TF1(Form("f%s_p%d_st",h->GetName(),p),f2gammastep,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
      fus[p][2] = new TF1(Form("f%s_p%d_g0",h->GetName(),p),f2gammagaus0,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
      fus[p][3] = new TF1(Form("f%s_p%d_g1",h->GetName(),p),f2gammagaus1,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
      fus[p][4] = new TF1(Form("f%s_p%d_bs",h->GetName(),p),f2gammabgstep,xpeaks[p]-frange[0],xpeaks[p]+frange[1],9);
	

      fus[p][0]->SetLineColor(5);
      fus[p][1]->SetLineColor(5);
      fus[p][2]->SetLineColor(4);
      fus[p][3]->SetLineColor(4);
      fus[p][4]->SetLineColor(2);
      for(int k=0;k<5;k++){
	fus[p][k]->SetLineWidth(1);
	for(int l=0;l<9;l++)
	  fus[p][k]->SetParameter(l,fu[p]->GetParameter(l));
	fus[p][k]->Draw("same");
      if(p==0)
	func(fus[p][k],xpeaks[p]-frange[0],xpeaks[p]+frange[1],1000,"python/data");
      }
    }// draw

  }
  if(draw)
    h->GetXaxis()->SetRangeUser(500,2000);

   
  return rv;
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
Double_t f2gammabgstep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
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
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammagaus0(Double_t *x, Double_t *par){
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

  return result;

}
Double_t f2gammagaus1(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  double result = 0;
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);

  return result;

}
