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
char* fileCo = (char*)"./hist/hcal0475.root";
char* fileEu = (char*)"./hist/hcal0448.root";
char* fileBa = (char*)"./hist/hcal0468.root";
char* fileYy = (char*)"./hist/hcal0469.root";
char* fileBG = (char*)"./hist/hcal0470.root";
ofstream fout;

//char* fileEu = (char*)"/home/gamma20/rootfiles/run0455.root";
//char* fileBG = (char*)"/home/gamma20/rootfiles/run0470.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);

Double_t f2gammagaussbg(Double_t *x, Double_t *par);
Double_t f2gammabg(Double_t *x, Double_t *par);
Double_t f2gammastep(Double_t *x, Double_t *par);
Double_t f2gammagaus(Double_t *x, Double_t *par);
Double_t f2gammagaus0(Double_t *x, Double_t *par);
Double_t f2gammagaus1(Double_t *x, Double_t *par);

Double_t f3gammagaussbg(Double_t *x, Double_t *par);
Double_t f3gammabg(Double_t *x, Double_t *par);
Double_t f3gammastep(Double_t *x, Double_t *par);
Double_t f3gammagaus(Double_t *x, Double_t *par);
Double_t f3gammagaus0(Double_t *x, Double_t *par);
Double_t f3gammagaus1(Double_t *x, Double_t *par);
Double_t f3gammagaus2(Double_t *x, Double_t *par);
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
	for(UShort_t i=0;i<r.at(0).size();i++){
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
    for(UShort_t i=0;i<r.at(0).size();i++){
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
  cout << "scaling " << p << endl;
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
  intensity.open("./scripts/Eudecaydetail.dat");
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
    if(e>100 && i > 1 && i*140000./28.53 > h->GetBinContent(h->FindBin(e))){
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
      frange = 15;
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
      double s = (h->GetBinContent(h->FindBin(e-frange))-h->GetBinContent(h->FindBin(e+frange)))/(e-frange - (e+frange));
      fu->SetParameter(1,s);//bg slope
      fu->SetParLimits(1,0.1*s,2*s);//bg slope
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
	fu->SetParLimits(7,0.8,2.5);//sigma
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
void efficiency(int run, int source, double runtime=0, double act=0, bool AB = false){
  TFile *f = new TFile(Form("hist/hcal%04d.root",run));
  if(!f->IsOpen()){
    return;
  }
  
  TH2F* h2 = NULL;
  h2 = (TH2F*)f->Get("h_en_summary_fine");
  if(AB)
    h2 = (TH2F*)f->Get("hAB_en_summary_fine");
  if(h2==NULL)
    return; 
  TH1F* h = NULL;
  h = (TH1F*)h2->ProjectionY("hp");
  if(h==NULL)
    return;
  
  char* infile[4] = {(char*)"./scripts/Codecay_NNDC.dat",(char*)(char*)"./scripts/Ydecay_NNDC.dat",(char*)"./scripts/Eudecay_NNDC.dat",(char*)"./scripts/Badecay_NNDC.dat"};
  
  ifstream intensity;
  
  intensity.open(infile[source]);
  intensity.ignore(1000,'\n');
  vector<double> en;
  vector<double> in;
  vector<double> fc;
  vector<double> fe;
  int ctr = 0;
  TCanvas* ca2 = new TCanvas("ca2","ca2",1200,800);
  if(source==2)
    ca2->Divide(5,3);
  if(source==3)
    ca2->Divide(3,3);
  if(source==0 || source ==1)
    ca2->Divide(2,2);
  while(!intensity.eof()){
    h = (TH1F*)h2->ProjectionY("hp");
    if(h==NULL)
      return;
    double e,i,sp;
    intensity >> e >> i >> sp;
    
    if(intensity.eof())
      break;
    intensity.ignore(1000,'\n');
    cout << e << "\t" << i << endl;
    en.push_back(e);
    in.push_back(i);
    ca2->cd(en.size());
    TF1* fu;
    TF1 *fus[3];
    frange = 15;
    if(source==2 && (e>410 && e<412))
      frange = 20;
    if(source==3 && (e>150 && e<250))
      frange = 7;
    h->GetXaxis()->SetRangeUser(e-frange,e+frange);
    h->GetYaxis()->UnZoom();
    //one peak
    int pars = 6;
    if(sp>0){
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
    //cout <<h->FindBin(e-frange) << "\t" << h->FindBin(e+frange) << "\t" << h->GetBinContent(h->FindBin(e-frange))  << "\t" << h->GetBinContent(h->FindBin(e+frange))  << "\t" << (h->GetBinContent(h->FindBin(e-frange))+h->GetBinContent(h->FindBin(e+frange)))*0.5 << endl;
    fu->SetParameter(0,(h->GetBinContent(h->FindBin(e-frange))+h->GetBinContent(h->FindBin(e+frange)))*0.5);//bg const
    double s = (h->GetBinContent(h->FindBin(e-frange))-h->GetBinContent(h->FindBin(e+frange)))/(e-frange - (e+frange));
    fu->SetParameter(1,0);//bg slope
    fu->SetParLimits(1,0,2*s);//bg slope
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
      fu->SetParameter(6,sp);//mean
      fu->SetParLimits(6,sp-1,sp+1);//mean
      fu->SetParameter(7,1);//sigma
      fu->SetParLimits(7,0.8,2.5);//sigma
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
  }
  for(UShort_t i=0;i<en.size();i++){
    fout << en[i] << "\t" << fc[i] << "\t" << fe[i] << "\t" << source<< endl;
    cout << en[i] << "\t" << fc[i] << "\t" << fe[i] << "\t" << source<< endl;
  }
  //cout << en.size() << endl;
  TGraphErrors* g = new TGraphErrors(en.size(),&en[0],&fc[0],0,&fe[0]);
  ca2->cd(en.size()+1);
  g->Draw("AP");
  f->Close();
}
void efficiency_twocomp(int run, int source, double runtime=0, double act=0, bool AB = false){
  TFile *f = new TFile(Form("hist/hcal%04d.root",run));
  if(!f->IsOpen()){
    return;
  }
  
  TH2F* h2 = NULL;
  h2 = (TH2F*)f->Get("h_en_summary_fine");
  if(AB)
    h2 = (TH2F*)f->Get("hAB_en_summary_fine");
  if(h2==NULL)
    return; 
  TH1F* h = NULL;
  h = (TH1F*)h2->ProjectionY("hp");
  if(h==NULL)
    return;
  
  char* infile[4] = {(char*)"./scripts/Codecay_NNDC.dat",(char*)(char*)"./scripts/Ydecay_NNDC.dat",(char*)"./scripts/Eudecay_NNDC.dat",(char*)"./scripts/Badecay_NNDC.dat"};
  
  ifstream intensity;
  
  intensity.open(infile[source]);
  intensity.ignore(1000,'\n');
  vector<double> en;
  vector<double> in;
  vector<double> fc;
  vector<double> fe;
  int ctr = 0;
  TCanvas* ca2 = new TCanvas("ca2","ca2",1200,800);
  if(source==2)
    ca2->Divide(5,3);
  if(source==3)
    ca2->Divide(3,3);
  if(source==0 || source ==1)
    ca2->Divide(2,2);
  while(!intensity.eof()){
    h = (TH1F*)h2->ProjectionY("hp");
    if(h==NULL)
      return;
    double e,i,sp;
    intensity >> e >> i >> sp;
    
    if(intensity.eof())
      break;
    intensity.ignore(1000,'\n');
    cout << e << "\t" << i << endl;
    en.push_back(e);
    in.push_back(i);
    ca2->cd(en.size());
    TF1* fu;
    TF1 *fus[6];
    frange = 15;
    if(source==2 && (e>110 && e<135))
      frange = 10;
    if(source==2 && (e>410 && e<412))
      frange = 20;
    if(source==2 && (e>770 && e<810))
      frange = 20;
    if(source==3 && (e>150 && e<250))
      frange = 8;
//    if(source==1 && AB)
//      frange = 9;
    cout << "frange " << frange << endl;
    h->GetXaxis()->SetRangeUser(e-frange,e+frange);
    h->GetYaxis()->UnZoom();
    //one peak
    int pars = 9;
    if(sp>0){
      pars = 12;
      fu = new TF1(Form("f%s_p%d",h->GetName(),(int)e),f3gammagaussbg,e-frange,e+frange,pars);
      fus[0] = new TF1(Form("f%s_p%d_bg",h->GetName(),(int)e),f3gammabg,e-frange,e+frange,pars);
      fus[1] = new TF1(Form("f%s_p%d_st",h->GetName(),(int)e),f3gammastep,e-frange,e+frange,pars);
      fus[2] = new TF1(Form("f%s_p%d_ga",h->GetName(),(int)e),f3gammagaus,e-frange,e+frange,pars);
      fus[3] = new TF1(Form("f%s_p%d_g0",h->GetName(),(int)e),f3gammagaus0,e-frange,e+frange,pars);
      fus[4] = new TF1(Form("f%s_p%d_g1",h->GetName(),(int)e),f3gammagaus1,e-frange,e+frange,pars);
      fus[5] = new TF1(Form("f%s_p%d_g2",h->GetName(),(int)e),f3gammagaus2,e-frange,e+frange,pars);
    }
    else{
      fu = new TF1(Form("f%s_p%d",h->GetName(),(int)e),f2gammagaussbg,e-frange,e+frange,pars);
      fus[0] = new TF1(Form("f%s_p%d_bg",h->GetName(),(int)e),f2gammabg,e-frange,e+frange,pars);
      fus[1] = new TF1(Form("f%s_p%d_st",h->GetName(),(int)e),f2gammastep,e-frange,e+frange,pars);
      fus[2] = new TF1(Form("f%s_p%d_ga",h->GetName(),(int)e),f2gammagaus,e-frange,e+frange,pars);
      fus[3] = new TF1(Form("f%s_p%d_g0",h->GetName(),(int)e),f2gammagaus0,e-frange,e+frange,pars);
      fus[4] = new TF1(Form("f%s_p%d_g1",h->GetName(),(int)e),f2gammagaus1,e-frange,e+frange,pars);
    }
    fu->SetLineColor(3);
    fu->SetLineWidth(1);
    //cout <<h->FindBin(e-frange) << "\t" << h->FindBin(e+frange) << "\t" << h->GetBinContent(h->FindBin(e-frange))  << "\t" << h->GetBinContent(h->FindBin(e+frange))  << "\t" << (h->GetBinContent(h->FindBin(e-frange))+h->GetBinContent(h->FindBin(e+frange)))*0.5 << endl;
    fu->SetParameter(0,(h->GetBinContent(h->FindBin(e-frange))+h->GetBinContent(h->FindBin(e+frange)))*0.5);//bg const
    double s = (h->GetBinContent(h->FindBin(e-frange))-h->GetBinContent(h->FindBin(e+frange)))/(e-frange - (e+frange));
    fu->SetParameter(1,0);//bg slope
    fu->SetParLimits(1,-fabs(5*s),fabs(5*s));//bg slope
    fu->SetParameter(2,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange))/2);//norm
    fu->SetParLimits(2,0,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange)));//norm
    fu->SetParameter(3,e);//mean
    fu->SetParLimits(3,e-2,e+2);//mean
    fu->SetParameter(4,1);//sigma
    fu->SetParLimits(4,0.1,3);//sigma
    fu->SetParameter(5,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange))/2);//norm
    fu->SetParLimits(5,0,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange)));//norm
    fu->SetParameter(6,e);//mean
    fu->SetParLimits(6,e-2,e+2);//mean
    if(source==2)
      fu->SetParLimits(6,e-0.2,e+0.2);//mean
    fu->SetParameter(7,5);//sigma
    fu->SetParLimits(7,2,7);//sigma
    if(source==3 && (e>150 && e<250)){
      fu->SetParLimits(4,0,1.1);//sigma
      fu->SetParLimits(7,1,2.5);//sigma
      fu->SetParLimits(6,e-0.2,e+0.2);//mean
      fu->SetParLimits(3,e-0.2,e+0.2);//mean
    }
    if(source==2 && (e>200 && e<250)){
      fu->SetParLimits(4,1,3.1);//sigma
      fu->SetParLimits(7,0,1.1);//sigma
      fu->SetParLimits(6,e-0.2,e+0.2);//mean
    }
    if(e>1000)
      fu->SetParLimits(7,2,10);//sigma
    if(source==2 && AB && e>1100 &&e<1120){
      fu->SetParLimits(7,5,8);//sigma
      //fu->SetParLimits(1,-1e7,0);//sigma
      fu->SetRange(e-15,e+23);
    }
    if(source==2 && AB && e>1200 &&e<1220){
      fu->SetParLimits(1,-1e7,0);//sigma
    }
//    if(source==1 && AB){
//      fu->SetParLimits(6,e-0.2,e+0.2);//mean
//      fu->SetParLimits(3,e-0.2,e+0.2);//mean
//      fu->SetRange(e-9,e+15);
//      if(e<1000)
//	fu->SetRange(e-8,e+15);
//    }
   
    if(pars==9){
      fu->SetParameter(8,h->GetBinContent(h->FindBin(e-frange))/2);//step
      fu->SetParLimits(8,0.0,h->GetBinContent(h->FindBin(e-frange)));//step
    }
    if(pars==12){
      fu->SetParameter(11,h->GetBinContent(h->FindBin(e-frange))/4);//step
      fu->SetParLimits(11,0.0,h->GetBinContent(h->FindBin(e-frange))/2);//step
      fu->SetParameter(8,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange))/20);//norm
      fu->SetParameter(9,sp);//mean
      fu->FixParameter(9,sp);//mean
      //fu->SetParLimits(9,sp-1,sp+1);//mean
      fu->SetParameter(10,1);//sigma
      fu->SetParLimits(10,0.8,2.5);//sigma
      if(source==2 && (e>420 && e<450)){
	fu->SetParLimits(4,0,1.1);//sigma
	fu->SetParLimits(7,1,2.5);//sigma
	fu->SetParLimits(6,442.8-0.2,442.8+0.2);//mean
      }
      if(source==2 && (e>1089 && e<1090)){
	//fu->SetParLimits(4,0,1.5);//sigma
	fu->SetParLimits(1,0,10);//bg slope
	fu->SetParLimits(7,3,3.8);//sigma
	fu->SetParLimits(6,e-0.02,e+0.02);//mean
	fu->SetRange(e-10,e+8);
	//fu->FixParameter(5,0);//mean
      }
    }
    h->Fit(fu,"Rn");
    h->DrawCopy();
    fu->Draw("same");
	

    fus[0]->SetLineColor(5);
    fus[1]->SetLineColor(4);
    fus[2]->SetLineColor(2);
    fus[3]->SetLineColor(6);
    fus[4]->SetLineColor(6);
    int kk = 5;
    if(pars==12)
      kk=6;
    for(int k=0;k<kk;k++){
      fus[k]->SetLineWidth(1);
      for(int j=0;j<pars;j++)
	fus[k]->SetParameter(j,fu->GetParameter(j));
      fus[k]->Draw("same");
    }
    if(runtime>0 && act >0){
      fc.push_back((fu->GetParameter(2) + fu->GetParameter(5))/h->GetBinWidth(1)/(i/100)/runtime/act);
      fe.push_back(sqrt(fu->GetParError(2)*fu->GetParError(2)+fu->GetParError(5)*fu->GetParError(5))/h->GetBinWidth(1)/(i/100)/runtime/act);
      // if(source==2 && e> 1089 && e< 1090){
      // 	cout << "asdasdffsdf  " <<fu->GetParameter(2) << "\t" << fu->GetParameter(2)/h->GetBinWidth(1)/(i/100)/runtime/act << endl;
      // 	fc.back() = fu->GetParameter(2)/h->GetBinWidth(1)/(i/100)/runtime/act;
      // 	fe.back() = fu->GetParError(2)/h->GetBinWidth(1)/(i/100)/runtime/act;
      // }
    }
    else{
      fc.push_back((fu->GetParameter(2) + fu->GetParameter(5))/h->GetBinWidth(1)/i);
      fe.push_back(sqrt(fu->GetParError(2)*fu->GetParError(2)+fu->GetParError(5)*fu->GetParError(5))/h->GetBinWidth(1)/i);
    }
  }
  for(UShort_t i=0;i<en.size();i++){
    fout << en[i] << "\t" << fc[i] << "\t" << fe[i] << "\t" << source<< endl;
    cout << en[i] << "\t" << fc[i] << "\t" << fe[i] << "\t" << source<< endl;
  }  
  //cout << en.size() << endl;
  TGraphErrors* g = new TGraphErrors(en.size(),&en[0],&fc[0],0,&fe[0]);
  ca2->cd(en.size()+1);
  g->Draw("AP");
  f->Close();
}
void checkruns(){
  int r[4] = {475,469,455,468};
  TFile *f;
  TH2F* h2 = NULL;
  TH1F* h = NULL;
  int ctr=0;
  for(int i=0;i<4;i++){
    f = new TFile(Form("hist/hcal%04d.root",r[i]));
    if(!f->IsOpen())
      continue;
  
    h2 = (TH2F*)f->Get("h_en_summary");
    if(h2==NULL)
      continue;
    h = (TH1F*)h2->ProjectionX(Form("hp%d",r[i]),30,4000);
    if(h==NULL)
      continue;
    h->Scale(1./h->Integral());
    h->SetLineColor(1+i);
    if(i==0)
      h->DrawCopy();
    else
      h->DrawCopy("same");
    f->Close();
  }
}
void effcurve(){
  fout.open("python/data/effcurve_twocomp.dat");
  efficiency_twocomp(475, 0, 3735.0, 23484.2);
  efficiency_twocomp(469, 1, 2703.1, 26565.8);
  efficiency_twocomp(455, 2, 3629.2,  6349.4);
  efficiency_twocomp(468, 3, 2060.9, 50773.6);
  fout.close();
  fout.open("python/data/effcurveAB_twocomp.dat");
  efficiency_twocomp(475, 0, 3735.0, 23484.2,1); // Co
  efficiency_twocomp(469, 1, 2703.1, 26565.8,1); // Y
  efficiency_twocomp(455, 2, 3629.2,  6349.4,1); // Eu 
  efficiency_twocomp(468, 3, 2060.9, 50773.6,1); // Ba
  fout.close();
  /*
  fout.open("python/data/effcurve_onecomp.dat");
  efficiency(475, 0, 3735.0, 23484.2);
  efficiency(469, 1, 2703.1, 26565.8);
  efficiency(455, 2, 3629.2,  6349.4);
  efficiency(468, 3, 2060.9, 50773.6);
  fout.close();
  fout.open("python/data/effcurveAB_onecomp.dat");
  efficiency(475, 0, 3735.0, 23484.2,1);
  efficiency(469, 1, 2703.1, 26565.8,1);
  efficiency(455, 2, 3629.2,  6349.4,1);
  efficiency(468, 3, 2060.9, 50773.6,1);
  fout.close();
  */
}
Double_t effcurve(Double_t *x, Double_t *p){
  double xx = log(x[0]);
  return exp(p[0] + p[1]*xx + p[2]*pow(xx,2) + p[3]*pow(xx,3));
}
Double_t effcurve2(Double_t *x, Double_t *p){
  double xx = log(x[0]/100);
  double yy = log(x[0]/1000);
  return exp( (p[0] + p[1]*xx + p[2]*pow(xx,2)) + (p[3]*yy + p[3]*pow(yy,2)));
}
Double_t effcurve3(Double_t *x, Double_t *p){
  double xx = log(x[0]/100);
  double yy = log(x[0]/1000);
  return exp( pow(  pow(p[0] + p[1]*xx + p[2]*pow(xx,2),-p[6]) + pow(p[3] + p[4]*yy + p[5]*pow(yy,2),-p[6]) , -1./p[6]) );
}
void fiteff(){
  ifstream data[2];
  //data[0].open("python/data/effcurve_twocomp.dat");
  //data[1].open("python/data/effcurve_onecomp.dat");
  data[0].open("python/data/effcurveAB_twocomp.dat");
  data[1].open("python/data/bkp/effcurveAB_twocomp.dat");
  vector<double> ve[2],vf[2],vu[2];
  TGraphErrors* g[2];
  TF1 *ff3[2];
  for(int j=0;j<2;j++){
    while(!data[j].eof()){
      double e,f,u;
      data[j] >> e >> f >> u;
      data[j].ignore(100,'\n');
      if(data[j].eof())
	break;
      ve[j].push_back(e);
      vf[j].push_back(f*100);
      vu[j].push_back(u*100);
    }
    g[j] = new TGraphErrors(ve[j].size(),&ve[j][0],&vf[j][0],0,&vu[j][0]);
    g[j]->SetMarkerStyle(20+j);
    if(j==0)
      g[j]->Draw("AP");
    else
      g[j]->Draw("P");
    // TF1 *ff = new TF1("effcurve",effcurve,0,2000,4);
    // ff->SetParameters(2,-0.2,0.02,-0.002);
    // g->Fit(ff,"Rn");
    // ff->Draw("same");
    // TF1 *ff2 = new TF1("effcurve2",effcurve2,0,2000,5);
    // ff2->SetParameters(1,-0.2,0.02,-0.002,0.02,-0.002);
    // g->Fit(ff2,"Rn");
    // ff2->Draw("same");
    ff3[j] = new TF1(Form("effcurve3_%d",j),effcurve3,0,2000,7);
    ff3[j]->SetParameter(0,1);
    ff3[j]->SetParameter(1,1);
    ff3[j]->SetParameter(2,0);
    ff3[j]->SetParameter(3,0.5);
    ff3[j]->SetParameter(4,-0.5);
    ff3[j]->SetParameter(5,0.05);
    ff3[j]->SetParameter(6,10);
    ff3[j]->FixParameter(2,0);
    ff3[j]->FixParameter(6,10);
    g[j]->Fit(ff3[j],"Rn");
    //ff3->Draw();
    ff3[j]->Draw("same");
  }
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
// three peak fitting function
Double_t f3gammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];

  Double_t step = par[11];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);
  arg = (x[0]-mean3)/(sqrt2*sigma3);
  result += 1/(sqrt2pi*sigma3) * norm3 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f3gammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t f3gammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];

  Double_t step = par[11];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean3)/(sqrt2*sigma3);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f3gammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];
  double result = 0;
  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);
  arg = (x[0]-mean3)/(sqrt2*sigma3);
  result += 1/(sqrt2pi*sigma3) * norm3 * exp(-arg*arg);

  return result;

}
Double_t f3gammagaus0(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];
  double result = 0;
  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);

  return result;

}
Double_t f3gammagaus1(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];
  double result = 0;
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);

  return result;

}
Double_t f3gammagaus2(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  Double_t norm3  = par[8];
  Double_t mean3  = par[9];
  Double_t sigma3 = par[10];
  double result = 0;
  arg = (x[0]-mean3)/(sqrt2*sigma3);
  result += 1/(sqrt2pi*sigma3) * norm3 * exp(-arg*arg);

  return result;

}
