#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "scripts/fit.h"
double crange=5;
double frange=25;
TH1F* h;
TFile *f = NULL;
vector<double> eff(int e, bool AB, bool draw = false);
vector<double> eff(int e, int m, int c, bool AB, bool draw = false);
double simevt = 1e6;
void effcurve(bool AB = false){
  vector<double> en;
  vector<double> ef;
  vector<double> ee;
  //frange
  ofstream fout;
  if(!AB)
    fout.open("python/data/simeff_scaled.dat");
  else
    fout.open("python/data/simeffAB_scaled.dat");
  
  for(int e=50; e<500; e+=50){
    f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
    vector<double> fitv = eff(e,AB);
    en.push_back(e*1.0);
    cout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    fout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    ef.push_back(fitv.at(0));
    ee.push_back(fitv.at(1));
  }
  for(int e=500; e<1000; e+=100){
    f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
    vector<double> fitv = eff(e,AB);
    en.push_back(e*1.0);
    cout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    fout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    ef.push_back(fitv.at(0));
    ee.push_back(fitv.at(1));
  }
  for(int e=1000; e<2100; e+=200){
    f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
    vector<double> fitv = eff(e,AB);
    en.push_back(e*1.0);
    cout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    fout << e << "\t" << fitv.at(0) << "\t" << fitv.at(1) << endl;
    ef.push_back(fitv.at(0));
    ee.push_back(fitv.at(1));
  }
  TGraphErrors* g = new TGraphErrors(en.size(),&en[0],&ef[0],0,&ee[0]);
  g->Draw("AP*");
  fout.close();
}
void checksum(int e, bool AB=false){
  double esum = 0;
  double eser = 0;
  double csum = 0;
  f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
  for(int m=0;m<10;m++){
    for(int c=0;c<4;c++){
      crange=5;
      if(m==4)
	crange=10;
      vector<double> fitv = eff(e,m,c,AB);
      if(fitv.at(0)>0){
	cout << m << "\t" << c << "\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2) << endl;
	esum+=fitv.at(0);
	eser+=fitv.at(1)*fitv.at(1);
	csum+=fitv.at(2);
      }
    }
  }
  crange=10;
  vector<double> fitv = eff(e,AB);
  if(fitv.at(0)>0){
    cout << "total\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2) << endl;
  }
  cout << "sum\t" << esum << "\t" << sqrt(eser) << "\t" << csum << endl;

}
void crystals(int e, bool AB=false){
  ofstream fout;
  if(!AB)
    fout.open(Form("python/data/simefficiencies_e%04d.dat",e));
  else
    fout.open(Form("python/data/simefficienciesAB_e%04d.dat",e));
  f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
  for(int m=0;m<10;m++){
    for(int c=0;c<4;c++){
      vector<double> fitv = eff(e,m,c,AB);
      if(fitv.at(0)>0){
	cout << m << "\t" << c << "\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2) << endl;
	fout << m << "\t" << c << "\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2) << endl;
	
      }
    }
  }
  fout.close();
}
void radii(bool AB=false){
  char* fnames[8] = {"source1173keV_MB35_SC34",
		     "source1173keV_MB34_SC33",
		     "source1173keV_MB33_SC32",
		     "source1173keV_MB32_SC31",
		     "source1173keV_MB31_SC30",
		     "source1173keV_MB30_SC29",
		     "source1173keV_MB29_SC28",
		     "source1173keV_MB28_SC27"};
  ofstream fout;
  for(int d=0;d<8;d++){
    if(!AB)
      fout.open(Form("python/data/simeff/%s.dat",fnames[d]));
    else
      fout.open(Form("python/data/simeff/%s_AB.dat",fnames[d]));
    f = new TFile(Form("simulation/hist/%s.root",fnames[d]));
    for(int m=0;m<10;m++){
      for(int c=0;c<4;c++){
	vector<double> fitv = eff(1173,m,c,AB);
	if(fitv.at(0)>0){
	  cout << m << "\t" << c << "\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2)/simevt << endl;
	  fout << m << "\t" << c << "\t" << fitv.at(0) << "\t" << fitv.at(1) << "\t" << fitv.at(2)/simevt << endl;
	
	}
      }
    }
    fout.close();
  }
}
vector<double> eff(int e, bool AB, bool draw){
  if(f==NULL)
    f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
  TH2F *h2;
  if(AB)
    h2 = (TH2F*)f->Get("egamAB_summary_fine");
  else
    h2 = (TH2F*)f->Get("egam_summary_fine");
  h = (TH1F*)h2->ProjectionY(Form("hp%04d",e));
  TF1 *fu;
  TF1 *fus[5];
  if(draw){
    h->Draw();
  }
  vector<double> rv;
  rv.clear();
    
  h->GetXaxis()->SetRangeUser(e-frange*2,e+frange*2);
  fu = new TF1(Form("f%s",h->GetName()),f2gammagaussbg,e-frange,e+frange,9);
  fu->SetLineColor(3);
  fu->SetLineWidth(1);
  fu->SetParameter(0,0);//bg const
  fu->FixParameter(1,0);//bg slope
  fu->SetParameter(2,h->Integral(e-frange,e+frange)/2);//norm
  fu->SetParLimits(2,1,1e6);
  fu->SetParameter(3,e);//mean
  fu->SetParLimits(3,e-1,e+1);//mean
  fu->SetParameter(4,1);//sigma
  fu->SetParLimits(4,0.1,5);//sigma
  fu->SetParameter(5,h->Integral(e-frange,e+frange)/2);//norm
  fu->SetParLimits(5,1,1e6);
  fu->SetParameter(6,e);//mean
  fu->SetParLimits(6,e-1,e+1);//mean
  fu->SetParameter(7,5);//sigma
  fu->SetParLimits(7,2,10);//sigma
  fu->SetParameter(8,h->GetBinContent(h->FindBin(e-frange)));//step
  // fu->FixParameter(0,0);
  // fu->FixParameter(1,0);
  if(e<100){
    fu->FixParameter(0,0);
    fu->FixParameter(1,0);
    fu->FixParameter(8,0);
    fu->FixParameter(5,0);
  }
  if(draw)
    h->Fit(fu,"Rn");
  else
    h->Fit(fu,"Rqn");

  rv.push_back((fu->GetParameter(2) + fu->GetParameter(5))/h->GetBinWidth(1)/simevt);
  rv.push_back(sqrt(fu->GetParError(2)*fu->GetParError(2)+fu->GetParError(5)*fu->GetParError(5))/h->GetBinWidth(1)/simevt);
  rv.push_back(h->Integral(h->FindBin(e-crange),h->FindBin(e+crange)));
  //draw results.
  if(draw){
    fu->Draw("same");
    fus[0] = new TF1(Form("f%s_bg",h->GetName()),f2gammabg,e-frange,e+frange,9);
    fus[1] = new TF1(Form("f%s_st",h->GetName()),f2gammastep,e-frange,e+frange,9);
    fus[2] = new TF1(Form("f%s_g0",h->GetName()),f2gammagaus0,e-frange,e+frange,9);
    fus[3] = new TF1(Form("f%s_g1",h->GetName()),f2gammagaus1,e-frange,e+frange,9);
    fus[4] = new TF1(Form("f%s_bs",h->GetName()),f2gammabgstep,e-frange,e+frange,9);
	

    fus[0]->SetLineColor(5);
    fus[1]->SetLineColor(5);
    fus[2]->SetLineColor(4);
    fus[3]->SetLineColor(4);
    fus[4]->SetLineColor(2);
    for(int k=0;k<5;k++){
      fus[k]->SetLineWidth(1);
      for(int l=0;l<9;l++)
	fus[k]->SetParameter(l,fu->GetParameter(l));
      fus[k]->Draw("same");
    }
    //h->GetXaxis()->SetRangeUser(0,2000);
  }//draw
   
  return rv;
}
vector<double> eff(int e, int m, int c, bool AB, bool draw){
  if(f==NULL)
    f = new TFile(Form("simulation/hist/source%04dkeV.root",e));
  TH2F *h2;
  if(AB)
    h2 = (TH2F*)f->Get("egamAB_summary_fine");
  else
    h2 = (TH2F*)f->Get("egam_summary_fine");
  h = (TH1F*)h2->ProjectionY(Form("hp%04d",e),4*m+c+1,4*m+c+1);
  TF1 *fu;
  TF1 *fus[4];
  if(draw){
    h->Draw();
  }
  vector<double> rv;
  rv.clear();
    
  h->GetXaxis()->SetRangeUser(e-frange*2,e+frange*2);
  fu = new TF1(Form("f%s",h->GetName()),fgammagaussbg,e-frange,e+frange,6);
  fu->SetLineColor(3);
  fu->SetLineWidth(1);
  fu->SetParameter(0,0);//bg const
  fu->FixParameter(1,0);//bg slope
  fu->SetParameter(2,h->Integral(h->FindBin(e-frange),h->FindBin(e+frange)));//norm
  fu->SetParLimits(2,1,1e6);
  fu->SetParameter(3,e);//mean
  fu->SetParLimits(3,e-1,e+1);//mean
  fu->SetParameter(4,1);//sigma
  fu->SetParLimits(4,0.1,5);//sigma
  fu->SetParameter(5,h->GetBinContent(h->FindBin(e-frange)));//step
  //fu->FixParameter(0,0);
  //fu->FixParameter(1,0);
  //fu->FixParameter(5,0);
  if(e<100){
    fu->FixParameter(0,0);
    fu->FixParameter(1,0);
    fu->FixParameter(5,0);
  }
  if(draw)
    h->Fit(fu,"Rn");
  else
    h->Fit(fu,"Rqn");

  rv.push_back(fu->GetParameter(2)/h->GetBinWidth(1)/simevt);
  rv.push_back(fu->GetParError(2)/h->GetBinWidth(1)/simevt);
  rv.push_back(h->Integral(h->FindBin(e-crange),h->FindBin(e+crange)));
  //draw results.
  if(draw){
    fu->Draw("same");
    fus[0] = new TF1(Form("f%s_bg",h->GetName()),fgammabg,e-frange,e+frange,6);
    fus[1] = new TF1(Form("f%s_st",h->GetName()),fgammastep,e-frange,e+frange,6);
    fus[2] = new TF1(Form("f%s_g0",h->GetName()),fgammagaus,e-frange,e+frange,6);
    fus[3] = new TF1(Form("f%s_bs",h->GetName()),fgammabgstep,e-frange,e+frange,6);
	

    fus[0]->SetLineColor(5);
    fus[1]->SetLineColor(5);
    fus[2]->SetLineColor(4);
    fus[3]->SetLineColor(2);
    for(int k=0;k<4;k++){
      fus[k]->SetLineWidth(1);
      for(int l=0;l<6;l++)
	fus[k]->SetParameter(l,fu->GetParameter(l));
      fus[k]->Draw("same");
    }
    //h->GetXaxis()->SetRangeUser(0,2000);
  }//draw
   
  return rv;
}
