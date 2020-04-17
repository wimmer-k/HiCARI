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
#include "scripts/fit.h"

using namespace std;
TFile* f = new TFile((char*)"./hist/hcal0469.root");
TCanvas *ca;
double* fitone(TH1F* h, bool draw = false);
double range[2] = {10,550};
double frange[2] = {10,150};
void test(int m, int c){
  TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  fitone(h,1);
}
void Thresh(){
  string abc[4] = {"A","B","C","D"};
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  vector<double> e0;
  vector<double> de;
  TEnv *cf = new TEnv("settings/thresholds0416.dat");
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      ca->cd(cry*12+1+clu);
      frange[1]=150;
      if(clu==0)
	frange[1] = 350;
      
      if(clu==4)
	frange[1] = 200;
      
      if((clu==8 || clu==9) && cry==0)
	frange[1] = 250;
      
      TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",clu,cry));
      double *p;
      p = fitone(h);
      if(p[0]<0 || p[1]<0)
	continue;
      e0.push_back(p[0]);
      de.push_back(p[1]);
      cf->SetValue(Form("Detector.%d.Crystal.%d.SimThreshold",clu,cry),1);
      cf->SetValue(Form("Detector.%d.Crystal.%d.E",clu,cry),p[0]);
      cf->SetValue(Form("Detector.%d.Crystal.%d.dE",clu,cry),p[1]);
    }
  }
  TGraph* g = new TGraph(e0.size(),&e0[0],&de[0]);
  TCanvas* ca2 = new TCanvas("ca2","ca2",1200,800);
  ca2->cd();
  g->Draw("AP*");
  cf->SaveLevel(kEnvLocal);
  
}
double* fitone(TH1F* h,  bool draw){
  static double rv[2] = {-1,-1};
  rv[0] = -1;
  rv[1] = -1;
  if(h==NULL)
    return rv;
  h->GetXaxis()->SetRangeUser(range[0],range[1]);
  if(h->Integral() < 10)
    return rv;
  if(draw){
    ca = new TCanvas("ca","ca",600,600);
    ca->cd();
  }
  TF1 *fus;
  TF1 *fu = new TF1(Form("f%s",h->GetName()),fthreshbg,frange[0],frange[1],4);
  fu->SetParameters(90,10,150,50);
  if(draw)
    h->Fit(fu,"Rn");
  else
    h->Fit(fu,"Rqn");
  //draw results.
  
  h->DrawCopy();
  fu->Draw("same");
  fus = new TF1(Form("f%s_bg",h->GetName()),fthresh,frange[0],frange[1],4);

  fus->SetLineColor(4);
  fus->SetLineWidth(1);
  for(int l=0;l<4;l++)
    fus->SetParameter(l,fu->GetParameter(l));
  fus->Draw("same");
 
  rv[0] = fu->GetParameter(0);
  rv[1] = fu->GetParameter(1);
  return rv;
  
}
