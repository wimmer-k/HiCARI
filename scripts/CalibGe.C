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
vector <vector<double> > defreturn2;
const int banks = 20;
const int slots = 12;
const int chans = 10;
const int corec = 9;
double resolution = 2;
double frange = 5000;
double range[2] = {150000,450000};
TCanvas *ca;
//char* fileCo = (char*)"/home/gamma20/rootfiles/hist0268.root";
char* fileCo = (char*)"/home/gamma20/rootfiles/run0420.root";
char* fileEu = (char*)"/home/gamma20/rootfiles/Eu_281_282.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t flinear(Double_t *x, Double_t *par);
vector<double> fitCo(TH1F* h, bool bg = true, bool draw = false);
vector <vector<double> > fitEu(TH1F* h, double roguh, bool bg = true, bool draw = false);

void SetRange(double low, double hig){
  range[0] = low;
  range[1] = hig;
}
void test(int b=11, int s = 3, int c =9){
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hmode3_en_bank%02d_slot%02d_chan%02d",b,s,c));
  fitCo(h,1,1);
}
void core(int m=0, int c =0){
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  if(h==NULL)
    return;
  vector<double> v = fitCo(h,1,1);
  cout << v[0] << "\t" << v[1] << endl;
  f = new TFile(fileEu);
  h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  frange = 3000;
  if(m==4)
    frange = 5000;
  fitEu(h,v[0],1,1);
}
void coreCo(int m=0, int c =0){
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  if(h==NULL)
    return;
  h->Draw();
  vector<double> v = fitCo(h,1,1);
}
void segment(int m=0, int c =0, int s =0){
  resolution = 5;
  TFile *f = new TFile(fileCo);
  TH2F* h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
  TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
  if(h==NULL || h->Integral()<100)
    return;
  fitCo(h,0,1);
}

vector <vector<double> > fitEu(TH1F* h, double rough, bool bg, bool draw){
  defreturn2.resize(4);
  defreturn2[0].push_back(-1);
  if(h->Integral() < 10)
    return defreturn2;
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
    //cout << en[p] << "\t" << en[p]/rough-frange << "\t" << en[p]/rough+frange << endl;
    h->GetXaxis()->SetRangeUser(en[p]/rough-frange,en[p]/rough+frange);
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
    if(!bg){
      fu->FixParameter(0,0);//bg const
      fu->FixParameter(1,0);//bg slope
    }
    fu->SetParameter(2,h->Integral(xpeaks[0]-frange,xpeaks[0]+frange));//norm
    fu->SetParameter(3,xpeaks[0]);//mean
    fu->SetParLimits(3,xpeaks[0]-500,xpeaks[0]+500);//mean
    fu->SetParameter(4,200);//sigma
    fu->SetParLimits(4,100,2000);//sigma
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
    (peaks.at(0)).push_back(fu->GetParameter(3)); // mean
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
    TGraphErrors *g = new TGraphErrors(peaks.at(0).size(),&peaks.at(0)[0],&peaks.at(5)[0],&peaks.at(1)[0]);
    g->Draw("AP*");
    g->Fit("pol1");    
    ca->cd(n+2);
    TGraphErrors *ge = new TGraphErrors(peaks.at(0).size(),&peaks.at(5)[0],&peaks.at(6)[0],0,&peaks.at(7)[0]);
    ge->Draw("AP*");
  }
  return peaks;
}
vector<double> fitCo(TH1F* h, bool bg, bool draw){
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
    fu[p]->SetParLimits(4,100,frange);//sigma
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
  }// peaksfu->GetParameter(2)/h->GetBinWidth(1)
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
void CalibGeEu(){
  string abc[3] = {"A","B","C"};
  frange = 3000;
  TFile *fco = new TFile(fileCo);                                                  
  TFile *feu = new TFile(fileEu);                                                  
  TEnv *cf = new TEnv(Form("%s.cal",fileEu));
  vector <double> gain;
  vector <double> offs;

  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(6,3);
  TCanvas *ca2 = new TCanvas("ca2","ca2",1200,800);
  ca2->Divide(6,3);
  TCanvas *ca3 = new TCanvas("ca3","ca3",1200,800);
  ca3->Divide(6,3);
  for(int clu=0;clu<6;clu++){
    for(int cry=0;cry<3;cry++){
      TH1F* h = (TH1F*)fco->Get(Form("hraw_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	return;
      vector<double> vco = fitCo(h,1,0);
      if(vco[0]<0)
	continue;
      h = (TH1F*)feu->Get(Form("hraw_en_clus%02d_crys%02d",clu,cry));
      frange = 3000;
      if(cry==4)
	frange = 5000;
      vector< vector<double> > veu = fitEu(h,vco[0],1,0);
      cout << clu << "\t" << cry << "\t" << abc[cry] << "\t" << veu[0][0] << endl;
      ca->cd(cry*6+1+clu);
      TGraphErrors* gc = new TGraphErrors(veu.at(0).size(), &veu.at(0)[0], &veu.at(5)[0], &veu.at(1)[0]);
      gc->SetTitle(Form("MB%d%s Core calibration",clu,(char*)abc[cry].c_str()));
      gc->Draw("AP*");
      TF1 *fl = new TF1("fl",flinear,0,500000,2);
      fl->SetParameters(0.004,1);
      gc->Fit(fl);
      ca2->cd(cry+1+clu*3);
      TGraphErrors* ge = new TGraphErrors(veu.at(0).size(), &veu.at(5)[0], &veu.at(6)[0],0, &veu.at(7)[0]);
      ge->SetTitle(Form("MB%d%s Core efficiency",clu,(char*)abc[cry].c_str()));
      ge->Draw("AP*");
      //ge->Fit("pol1");
      ca3->cd(cry+1+clu*3);
      TGraphErrors* gr = new TGraphErrors(veu.at(0).size(), &veu.at(5)[0], &veu.at(2)[0],0, 0);
      gr->SetTitle(Form("MB%d%s Core resolution",clu,(char*)abc[cry].c_str()));
      gr->Draw("AP*");
      
      //pritn calibration parameters
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),fl->GetParameter(0));
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),fl->GetParameter(1));
      gain.push_back(fl->GetParameter(0));
      offs.push_back(fl->GetParameter(1));
    }// crystal
  }//clu
  cf->SaveLevel(kEnvLocal);
  TCanvas *ca4 = new TCanvas("ca4","ca4",600,300);
  ca4->cd();
  TGraph* g = new TGraph(gain.size(),&gain[0],&offs[0]);
  g->Draw("AP*");

}
void CalibGeCo(){
  frange = 3000;
  TEnv *cf = new TEnv(Form("%s.cal",fileCo));
  TFile *f = new TFile(fileCo);
  vector <double> gain;
  vector <double> offs;
  vector <double> res0;
  vector <double> res1;
  
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      //cout << clu << "\t" << cry << endl;
      TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	continue;
      frange = 3000;
      SetRange(150000,450000);
      if(clu==5 && cry==1){
	frange = 10000;
      }
      if(clu==1 && cry==1){
	frange = 4000;
      }
      if(clu==11 && cry==1){
	frange = 30000;
	SetRange(1200000,1700000);
      }
      if(clu==10){
	frange = 5000;
	if(cry==0)
	  frange = 3000;
	if(cry==3)
	  frange = 8000;
	SetRange(400000,600000);
      }
      vector<double> r = fitCo(h,1);
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
      TH2F* h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",clu,cry));
      if(h2==NULL)
	continue;
      for(int s=0;s<40;s++){
	//cout << s << endl;
	TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
	if(h==NULL || h->Integral()<10)
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
	if(clu==11)
	  SetRange(420000,540000);
	if(clu==11&&cry==2&&s ==33)
	  SetRange(580000,740000);
	if(clu==10){
	  SetRange(400000,600000);
	  frange = 5000;
	  if(cry==0){
	    frange = 3000;
	    if(s==6)
	      SetRange(580000,840000);
	    if(s==6)
	      SetRange(580000,840000);
	    if(s==33)
	      SetRange(480000,620000);
	    if(s==34)
	      SetRange(440000,580000);
	  }
	  if(cry==1){
	    frange = 5000;
	  }
	  else
	    continue;
	  if(cry==3)
	    frange = 8000;
	}
	vector<double> r = fitCo(h,0);
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
  TFile *f = new TFile(fileCo);
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
	vector<double> r = fitCo(h);
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
Double_t flinear(Double_t *x, Double_t *par){
  return x[0]*par[0] + par[1];

}
