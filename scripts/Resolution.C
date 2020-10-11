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
char* fileCo = (char*)"./hist/hcal0475.root";
char* fileEu = (char*)"./hist/hcal0448.root";
char* fileBa = (char*)"./hist/hcal0468.root";
char* fileYy = (char*)"./hist/hcal0469.root";
TCanvas *ca;
TGraph* core(int m, int c, bool draw=false);
double fitonepeak(TH1F* h, double ene, bool draw = false);
double frange = 20;
void test(int m, int c, double e){
  TFile* f = new TFile(fileEu);
  TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  fitonepeak(h,e,1);
}
void ResCore(){
  string abc[4] = {"A","B","C","D"};
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  TEnv *cf = new TEnv("settings/resolutions0416.dat");
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      ca->cd(cry*12+1+clu);
      TGraph *g = core(clu,cry);
      if(!g)
	continue;
      g->SetTitle(Form("DET %d%s Core resolution",clu,(char*)abc[cry].c_str()));
      g->Draw("AP*");
      if(g->GetN()<3)
	continue;
      TF1 *fr = new TF1("fr",freso,0,3000,3);
      fr->SetParameters(1.15,0.0006,0);
      fr->FixParameter(2,0);
      g->Fit(fr,"q");
      double A = fr->GetParameter(0);
      double B = fr->GetParameter(1);
      cf->SetValue(Form("Detector.%d.Crystal.%d.SimResolution",clu,cry),1);
      cf->SetValue(Form("Detector.%d.Crystal.%d.A",clu,cry),A);
      cf->SetValue(Form("Detector.%d.Crystal.%d.B",clu,cry),B);
      
    }
  }
  cf->SaveLevel(kEnvLocal);
}
TGraph* core(int m, int c, bool draw){
  vector<double> ene;
  vector<double> res;
  const int nEu = 11;
  double enEu[nEu] = {121.77, 244.66, 344.28, 411.11, 443.96, 778.90, 867.3, 963.38, 1085.84, 1112.08, 1408.01};
  const int nCo = 2;
  double enCo[nCo] = {1173.23, 1332.49};
  const int nBa = 5;
  double enBa[nBa] = {81.00, 276.40, 302.85, 356.01, 383.85};
  const int nYy = 2;
  double enYy[nYy] = {898.04, 1836.06};

  frange = 20;

  TFile* f = new TFile(fileEu);
  TH1F* h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nEu;i++){
    double fr = fitonepeak(h,enEu[i],0);
    if(fr>0){
      ene.push_back(enEu[i]);
      res.push_back(fr);
    }
  }
  f = new TFile(fileCo);
  h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nCo;i++){
    double fr = fitonepeak(h,enCo[i],0);
    if(fr>0){
      ene.push_back(enCo[i]);
      res.push_back(fr);
    }
  }
  f = new TFile(fileBa);
  h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nBa;i++){
    double fr = fitonepeak(h,enBa[i],0);
    if(fr>0){
      if(m==7 && c==2 && fr>3)
	continue;
      ene.push_back(enBa[i]);
      res.push_back(fr);
    }
  }
  f = new TFile(fileYy);
  if(f->IsOpen()){
    h = (TH1F*)f->Get(Form("h_en_clus%02d_crys%02d",m,c));
    for(int i=0;i<nYy;i++){
      double fr = fitonepeak(h,enYy[i],0);
      if(fr>0){
	ene.push_back(enYy[i]);
	res.push_back(fr);
      }
    }
  }
  TGraph* g = new TGraph(ene.size(),&ene[0],&res[0]);
  TF1 *fr = new TF1("fr",freso,0,3000,3);
  fr->SetParameters(1.15,0.0006,0);
  fr->FixParameter(2,0);
  if(draw){
    ca = new TCanvas("ca","ca",800,800);
    ca->cd();
    g->Draw("AP*");
    g->Fit(fr);
  }
  else
    g->Fit(fr,"qn");

  f->Close();
  

  return g;

}

double fitonepeak(TH1F* h, double ene, bool draw){
  if(h==NULL)
    return -1;
  h->GetXaxis()->SetRangeUser(ene-frange,ene+frange);
  if(h->Integral() < 10)
    return -1;
  if(draw){
    ca = new TCanvas("ca","ca",600,600);
    ca->cd();
    h->DrawCopy();
  }
  TF1 *fu;
  TF1 *fus[3];
  h->GetXaxis()->SetRangeUser(ene-frange,ene+frange);
  fu = new TF1(Form("f%s",h->GetName()),fgammagaussbg,ene-frange,ene+frange,6);
  fu->SetLineColor(3);
  fu->SetLineWidth(1);
  fu->SetParameter(0,0);//bg const
  fu->SetParameter(1,0);//bg slope
  fu->SetParameter(2,h->Integral(ene-frange,ene+frange));//norm
  fu->SetParameter(3,ene);//mean
  fu->SetParLimits(3,ene-5,ene+5);//mean
  fu->SetParameter(4,1);//sigma
  fu->SetParLimits(4,0.1,20);//sigma
  fu->SetParameter(5,h->GetBinContent(h->FindBin(ene-frange)));//step
  if(draw)
    h->Fit(fu,"Rn");
  else
    h->Fit(fu,"Rqn");

  //draw results.
  if(draw){
    fu->Draw("same");
    fus[0] = new TF1(Form("f%s_bg",h->GetName()),fgammabg,ene-frange,ene+frange,6);
    fus[1] = new TF1(Form("f%s_st",h->GetName()),fgammastep,ene-frange,ene+frange,6);
    fus[2] = new TF1(Form("f%s_ga",h->GetName()),fgammagaus,ene-frange,ene+frange,6);
	

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
  return fu->GetParameter(4);
}
