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
double frange = 50;
double range[2] = {1500,4500};
TCanvas *ca;
char*   fileCo = (char*)"./rootfiles/mode3_0668_0670.root";
char*   fileEu = (char*)"./rootfiles/mode3_0667.root";
char*   fileBa = (char*)"./rootfiles/mode3_0673_0674.root";
char* fileYy = (char*)"./hist/hraw0469.root";
Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t fgammagausnostep(Double_t *x, Double_t *par);
Double_t flinear(Double_t *x, Double_t *par);
vector<double> fitCo(TH1F* h, bool bg = true, bool draw = false);
vector<double> fitCo2(TH1F* h, vector< vector<double> > &CoPar, bool bg = true, bool draw = false);
vector<double> fitCoBa(TH1F* h, double roguh, vector< vector<double> > CoPar, bool bg = true, bool draw = false);
vector <vector<double> > fitEu(TH1F* h, double roguh, bool bg = true, bool draw = false);
TGraph* core(int m, int c, bool draw=false);
TGraph* corerun(int m, int c, int run, bool draw=false);
double fitonepeak(TH1F* h, double low, double hig, bool draw = false);
ofstream fout;

TFile *f_fitEu=new TFile("/home/wimmer/source/hist/FitSpecEu.root", "recreate");

void SetRange(double low, double hig){
  range[0] = low;
  range[1] = hig;
}
void test(int b=11, int s = 3, int c =9){
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hmode3_en_bank%02d_slot%02d_chan%02d",b,s,c));
  fitCo(h,1,1);
}
void coreEu(int m=0, int c =0){
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  if(h==NULL)
    return;
  vector<double> v = fitCo(h,1,1);
  cout << v[0] << "\t" << v[1] << endl;
  f = new TFile(fileEu);
  h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  frange = 20;
  if(m==4 && c==2)
    frange = 30;
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
void segmentCoBa(int m=0, int c =0, int s =0){
  //resolution = 5;
  vector< vector <double> > CoPar;
  TFile *f = new TFile(fileCo);
  TH2F* h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
  TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
  if(h==NULL || h->Integral()<100)
    return;
  vector<double> v = fitCo2(h,CoPar,0,1);
  cout << v[0] << "\t" << v[1] << " CoPar.at(0).size()=" << CoPar.at(0).size() << endl;
  for(int i=0;i<CoPar.at(0).size();i++) cout << "Mean " << i << ": " << CoPar.at(0).at(i) << " ";
  cout << endl;
  f = new TFile(fileBa);
  h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
  h = (TH1F*)h2->ProjectionY(Form("%s_%02d",Form("%s_2",h2->GetName()),s),s+1,s+1);
  frange = 30;
  fitCoBa(h,v[0],CoPar,1,1);
}
void segmentEu(int m=0, int c =0, int s =0){
  //resolution = 5;
  TFile *f = new TFile(fileCo);
  TH2F* h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
  TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
  if(h==NULL || h->Integral()<100)
    return;
  vector<double> v = fitCo(h,0,1);
  cout << v[0] << "\t" << v[1] << endl;
  f = new TFile(fileEu);
  h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",m,c));
  h = (TH1F*)h2->ProjectionY(Form("%s_Eu_%02d",h2->GetName(),s),s+1,s+1);
  frange = 30;
  if(m==4)
    frange = 50;
  fitEu(h,v[0],1,1);
  //fitEu(h,0.415256,1,1);
}
void segment(int m=0, int c =0, int s =0){
  //resolution = 5;
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
  intensity.open("/home/wimmer/progs/HiCARI/scripts/Eudecay.dat");
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
    cout << en[p] << "\t" << en[p]/rough-frange << "\t" << en[p]/rough+frange << endl;
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
    Double_t *xpeaks = sp->GetPositionX();
    TF1 *fu;
    TF1 *fus[3];
    h->GetXaxis()->SetRangeUser(xpeaks[0]-3,xpeaks[0]+3);
    double rms=h->GetRMS();
    //cout << "rms=" << rms << "xpeaks[0]=" << xpeaks[0] << endl;
    h->GetXaxis()->SetRangeUser(xpeaks[0]-frange,xpeaks[0]+frange);
    ////////////////////////////////////////////////////////////////////////////////////////////// Addition BM
    double norm,mean;
    fu = new TF1(Form("f%s_p%d",h->GetName(),p),"gaus",xpeaks[0]-rms,xpeaks[0]+rms);
    fu->SetParameter(0,h->Integral(xpeaks[0]-frange,xpeaks[0]+frange));//norm
    fu->SetParameter(1,xpeaks[0]);//mean
    fu->SetParLimits(1,xpeaks[0]-frange/2.,xpeaks[0]+frange/2.);//mean
    fu->SetParameter(2,rms);//sigma
    fu->SetParLimits(2,rms/2.,rms*2);//sigma
    h->Fit(fu,"Rqn");
    norm=fu->GetParameter(0);
    //cout << "norm=" << norm << endl;
    mean=fu->GetParameter(1);
    rms=fu->GetParameter(2);
 
    delete fu;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    fu = new TF1(Form("f%s_p%d",h->GetName(),p),fgammagausnostep,xpeaks[0]-frange,xpeaks[0]+frange,6);
    fu->SetLineColor(3);
    fu->SetLineWidth(1);
    fu->SetParameter(0,0);//bg const
    fu->SetParameter(1,0);//bg slope
    if(!bg){
      fu->FixParameter(0,0);//bg const
      fu->FixParameter(1,0);//bg slope
    }
    //////////////////////////// Change BM
    fu->SetParameter(2,norm);//norm
    fu->SetParameter(3,mean);//mean
    fu->SetParLimits(3,mean-2,mean+2);//mean
    fu->SetParameter(4,rms);//sigma//*/
    //////////////////////////// 
    /*fu->SetParameter(2,h->Integral(xpeaks[0]-frange,xpeaks[0]+frange));//norm
    fu->SetParameter(3,xpeaks[0]);//mean
    fu->SetParLimits(3,xpeaks[0]-frange/2.,xpeaks[0]+frange/2.);//mean
    fu->FixParameter(4,rms);//sigma//*/
    fu->SetParLimits(4,rms-4,rms+4);//sigma
    //  fu->SetParameter(5,h->GetBinContent(h->FindBin(xpeaks[0]-frange)));//step
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
      // fus[2] = new TF1(Form("f%s_p%d_st",h->GetName(),p),fgammastep,xpeaks[0]-frange,xpeaks[0]+frange,6);
      fus[1] = new TF1(Form("f%s_p%d_ga",h->GetName(),p),fgammagaus,xpeaks[0]-frange,xpeaks[0]+frange,6);
	

      fus[0]->SetLineColor(5);
      // fus[2]->SetLineColor(4);
      fus[1]->SetLineColor(2);
      for(int k=0;k<2;k++){
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
    //ca->SaveAs(Form("/home/gamma20/fall2020/plots/FitEu/fitEu_%s.png",h->GetName()));
    ca->SetName(h->GetName());
    f_fitEu->cd();
    ca->Write();
    //ca->Close();
    ca->WaitPrimitive();
    // delete ca;
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
  Double_t *xpeaks = sp->GetPositionX();
  Double_t *ypeaks = sp->GetPositionY();

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
    Double_t temp = xpeaks[1];
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
    fu[p]->SetParLimits(3,xpeaks[p]-5,xpeaks[p]+5);//mean
    fu[p]->SetParameter(4,2);//sigma
    fu[p]->SetParLimits(4,1,frange);//sigma
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
      cout << "New frange=" << frange << endl;
      //ca->WaitPrimitive();

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
    cout << "New frange=" << frange << endl;
    //    ca->WaitPrimitive();
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
vector<double> fitCo2(TH1F* h, vector< vector<double> > &peaks, bool bg, bool draw){
  peaks.resize(6);
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
  Double_t *xpeaks = sp->GetPositionX();
  Double_t *ypeaks = sp->GetPositionY();

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
    Double_t temp = xpeaks[1];
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
    fu[p]->SetParLimits(3,xpeaks[p]-5,xpeaks[p]+5);//mean
    fu[p]->SetParameter(4,2);//sigma
    fu[p]->SetParLimits(4,1,frange);//sigma
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
      cout << "New frange=" << frange << endl;
      ca->WaitPrimitive();

    }// draw
  }// peaksfu->GetParameter(2)/h->GetBinWidth(1)
  double y[2] = {1173.2,1332.5};
  if(draw){
    ca->cd(nfound+2);
    double x[2] = {fu[0]->GetParameter(3), fu[1]->GetParameter(3)};
    double e[2] = {fu[0]->GetParError(3), fu[1]->GetParError(3)};
 
    TGraphErrors *g = new TGraphErrors(2,x,y,e);
    g->Draw("AP");
    g->Fit("pol1");
    cout << "New frange=" << frange << endl;
    // ca->WaitPrimitive();
  }
  for(int ij=0;ij<2;ij++)
    {
      (peaks.at(0)).push_back(fu[ij]->GetParameter(3)); // mean
      (peaks.at(1)).push_back(fu[ij]->GetParError(3));  // error mean
      (peaks.at(2)).push_back(fu[ij]->GetParameter(4)); // sigma
      (peaks.at(3)).push_back(fu[ij]->GetParameter(2)/h->GetBinWidth(1)); // content
      (peaks.at(4)).push_back(fu[ij]->GetParError(2)/h->GetBinWidth(1)); // content error 
      (peaks.at(5)).push_back(y[ij]);
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
vector<double> fitCoBa(TH1F* h, double rough, vector< vector<double> > CoPar, bool bg, bool draw){
  cout << "CoPar.size()=" << CoPar.size() << endl;
  TSpectrum *sp = new TSpectrum(2,resolution);
  sp->SetResolution(resolution);
  const int n = 3;
  ifstream intensity;
  intensity.open("/home/wimmer/progs/HiCARI/scripts/Badecay.dat");
  intensity.ignore(1000,'\n');
  double en[n], inten[n];//, eff[n];
  if(draw){
    ca = new TCanvas("ca","ca",600,600);
    ca->Divide(4,3);
  }
  vector< vector<double> > peaks;
  peaks.resize(6);
  for(int i=0;i<6;i++)
    peaks.at(i).clear();
  for(int p=0; p<n; p++){
    intensity >> en[p] >> inten[p];
    cout << en[p] << "\t" << en[p]/rough-frange << "\t" << en[p]/rough+frange << endl;

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
    Double_t *xpeaks = sp->GetPositionX();
    TF1 *fu;
    TF1 *fus[3];
    h->GetXaxis()->SetRangeUser(xpeaks[0]-3,xpeaks[0]+3);
    double rms=h->GetRMS();
    //cout << "rms=" << rms << "xpeaks[0]=" << xpeaks[0] << endl;
    h->GetXaxis()->SetRangeUser(xpeaks[0]-frange,xpeaks[0]+frange);
    double norm,mean;
    fu = new TF1(Form("f%s_p%d",h->GetName(),p),"gaus",xpeaks[0]-rms,xpeaks[0]+rms);
    fu->SetParameter(0,h->Integral(xpeaks[0]-frange,xpeaks[0]+frange));//norm
    fu->SetParameter(1,xpeaks[0]);//mean
    fu->SetParLimits(1,xpeaks[0]-frange/2.,xpeaks[0]+frange/2.);//mean
    fu->SetParameter(2,rms);//sigma
    fu->SetParLimits(2,rms/2.,rms*2);//sigma
    h->Fit(fu,"Rqn");
    norm=fu->GetParameter(0);
    //cout << "norm=" << norm << endl;
    mean=fu->GetParameter(1);
    rms=fu->GetParameter(2);
 
    delete fu;
    fu = new TF1(Form("f%s_p%d",h->GetName(),p),fgammagausnostep,xpeaks[0]-frange,xpeaks[0]+frange,6);
    fu->SetLineColor(3);
    fu->SetLineWidth(1);
    fu->SetParameter(0,0);//bg const
    fu->SetParameter(1,0);//bg slope
    if(!bg){
      fu->FixParameter(0,0);//bg const
      fu->FixParameter(1,0);//bg slope
    }
    fu->SetParameter(2,norm);//norm
    fu->SetParameter(3,mean);//mean
    fu->SetParLimits(3,mean-2,mean+2);//mean
    fu->SetParameter(4,rms);//sigma//*/
    fu->SetParLimits(4,rms-4,rms+4);//sigma
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
      fus[1] = new TF1(Form("f%s_p%d_ga",h->GetName(),p),fgammagaus,xpeaks[0]-frange,xpeaks[0]+frange,6);
	

      fus[0]->SetLineColor(5);
      fus[1]->SetLineColor(2);
      for(int k=0;k<2;k++){
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
  }// peaks
  cout << "CoPar.size()=" << CoPar.size() << endl;
  for(int j=0;j<CoPar.size();j++)
    {
      cout << "CoPar.at(0).size()=" << CoPar.at(0).size() << endl;
      for(int i=0;i<CoPar.at(0).size();i++)
	(peaks.at(j)).push_back(CoPar.at(j).at(i));
    }
  
     TGraphErrors *g = new TGraphErrors(peaks.at(0).size(),&peaks.at(0)[0],&peaks.at(5)[0],&peaks.at(1)[0]);
     TF1 *fl = new TF1("fl","pol1",0,5000);
     fl->SetParameters(0.004,1);
     cout << "fl->GetParameter(0)=" << fl->GetParameter(0) << " fl->GetParameter(1)=" << fl->GetParameter(1) << endl;
     g->Fit(fl,"Rn");    
     if(draw){
       //cout << (peaks.at(0)).size() << endl;
       ca->cd(n+1);
       g->Draw("AP*");
       g->Fit(fl);    
       ca->cd(n+2);
       ca->Update();
       //ca->SaveAs(Form("/home/gamma20/fall2020/plots/FitEu/fitCoBa_%s.png",h->GetName()));
       ca->SetName(h->GetName());
       f_fitEu->cd();
       ca->Write();//*/
       //ca->Close();
       ca->WaitPrimitive();
       // delete ca;
     }
     cout << "write parameters peaks.at(0).size()=" << peaks.at(0).size() << endl;
     double a=fl->GetParameter(1);
     double b=fl->GetParameter(0);
     cout << a << "\t" << b << "\t" << peaks.at(0).size() << " " << h->GetName() << endl;
     vector<double> rv;
     rv.push_back(a);
     rv.push_back(b);
     rv.push_back(a*(peaks.at(0)).at(0));
     rv.push_back(a*(peaks.at(0)).at((peaks.at(0)).size()-1));
     return rv;
}
void CalibGeEu(){
  string abc[4] = {"A","B","C","D"};
  frange = 30;
  TFile *fco = new TFile(fileCo);                                                  
  TFile *feu = new TFile(fileEu);                                                  
  TEnv *cf = new TEnv(Form("%s.cal",fileEu));
  vector <double> gain;
  vector <double> offs;

  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  //  TCanvas *ca2 = new TCanvas("ca2","ca2",1200,800);
  //  ca2->Divide(12,4);
  // TCanvas *ca3 = new TCanvas("ca3","ca3",1200,800);
  // ca3->Divide(12,4);
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      TH1F* h = (TH1F*)fco->Get(Form("hraw_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	continue;
      ///////////// Added from CalibGeCo BM
      cout << "clu=" << clu << " cry=" << cry << endl;
      frange = 20;
      SetRange(1800,4500);
      if(clu==3 && cry==2){
	frange = 20; 
	SetRange(2000,3000);
	cout << "done" << endl;
      }
      /*if(clu==1 && cry==1){
       	frange = 20;
	}//*/
      if(clu==11 && cry==1){
       	frange = 100;
       	SetRange(9000,12000);
      }
      if(clu==10){
       	//frange = 15;
       	if(cry>0)
       	  frange = 30;
	//	if(cry==3)
	// frange = 80;
       	SetRange(3200,4600);
       }
      ///////////// End: Added from CalibGeCo
      vector<double> vco = fitCo(h,1,0);
      if(vco[0]<0)
	continue;
      h = (TH1F*)feu->Get(Form("hraw_en_clus%02d_crys%02d",clu,cry));
      if(h==NULL)
	continue;
      frange=25;
      if(cry==4)
	frange = 50;
      vector< vector<double> > veu = fitEu(h,vco[0],1,0);
      cout << "clu " << clu << "\tcry " << cry << "\t" << abc[cry] << "\t" << veu[0][0] << endl;
      ca->cd(cry*12+1+clu);
      TGraphErrors* gc = new TGraphErrors(veu.at(0).size(), &veu.at(0)[0], &veu.at(5)[0], &veu.at(1)[0]);
      gc->SetTitle(Form("MB%d%s Core calibration",clu,(char*)abc[cry].c_str()));
      gc->Draw("AP*");
      TF1 *fl = new TF1("fl",flinear,0,5000,2);
      float slope=(veu.at(1)[0]-veu.at(1)[(int)(veu.at(1).size()-1)])/(veu.at(0)[0]-veu.at(0)[(int)(veu.at(0).size()-1)]);
      fl->SetParameter(1,slope);
      gc->Fit(fl,"R");

      /*  ca2->cd(cry*12+1+clu);
      TGraphErrors* ge = new TGraphErrors(veu.at(0).size(), &veu.at(5)[0], &veu.at(6)[0],0, &veu.at(7)[0]);
      ge->SetTitle(Form("MB%d%s Core efficiency",clu,(char*)abc[cry].c_str()));
      ge->Draw("AP*");
      //ge->Fit("pol1");
      ca3->cd(cry*12+1+clu);
      TGraphErrors* gr = new TGraphErrors(veu.at(0).size(), &veu.at(5)[0], &veu.at(2)[0],0, 0);
      gr->SetTitle(Form("MB%d%s Core resolution",clu,(char*)abc[cry].c_str()));
      gr->Draw("AP*");//*/
      
      //pritn calibration parameters
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),fl->GetParameter(0));
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),fl->GetParameter(1));
      gain.push_back(fl->GetParameter(0));
      offs.push_back(fl->GetParameter(1));
      cout << "done clu=" << clu << " cry=" << cry << endl;
      //delete fl;
    }// crystal
  }//clu
  cf->SaveLevel(kEnvLocal);
  f_fitEu->Close();
  TCanvas *ca4 = new TCanvas("ca4","ca4",600,300);
  ca4->cd();
  TGraph* g = new TGraph(gain.size(),&gain[0],&offs[0]);
  g->Draw("AP*");
}
void CalibGeCo(){
  frange = 15;
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
      frange = 20;
      SetRange(1800,4500);
      if(clu==3 && cry==2){
	frange = 20; 
	SetRange(2000,3000);
	cout << "done" << endl;
      }
      /*if(clu==1 && cry==1){
       	frange = 20;
	}//*/
      if(clu==11 && cry==1){
       	frange = 100;
       	SetRange(9000,12000);
      }
      if(clu==11 && cry==2){
       	frange = 100;
       	SetRange(2500,3000);
      }
      if(clu==10){
       	//frange = 15;
       	if(cry>0)
       	  frange = 30;
	//	if(cry==3)
	// frange = 80;
       	SetRange(3200,4600);
       }
      vector<double> r = fitCo(h,1);
      //cout << " fitted " << r[0]<< endl;
      if(r[2]<1) cout << "///////////////////////!!!! clu=" << clu << " cry=" << cry << endl;
      if(r[0]<0)
	continue;
      gain.push_back(r[0]);
      offs.push_back(r[1]);
      res0.push_back(r[2]);
      res1.push_back(r[3]);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),r[0]);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),r[1]);//*/
    }
  }
  cout << "finished" << endl;
  cf->SaveLevel(kEnvLocal);
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
  /*  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      //cout << clu << "\t" << cry << endl;
      // cout << "clu=" << clu << " cry=" << cry << endl;
      TH2F* h2 = (TH2F*)f->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",clu,cry));
      if(h2==NULL)
	continue;
      for(int s=0;s<40;s++){
	resolution = 2;
	//cout << s << endl;
	TH1F* h = (TH1F*)h2->ProjectionY(Form("%s_%02d",h2->GetName(),s),s+1,s+1);
	if(h==NULL || h->Integral()<10)
	  continue;
	frange = 20;
	SetRange(1800,3100);
	if(clu==2 && cry==1 && s==1)
	  SetRange(2000,2500);
	/*if(clu==3 && cry==2 && s==4)
	  SetRange(2800,3800);
	if(clu==4 && cry==1 && s==2)
	  SetRange(2800,3800);
	if(clu==1 && cry==1 && s==2)
       	  frange = 20;
	if(clu==3 && cry==0 && s==3)
	frange = 20;// /
	if(clu>=6 && clu<=9) SetRange(2000,3000);
	if(clu==11)
	  {
	    SetRange(3000,4200);
	    frange=25;
	    if(cry==0&&s==33) SetRange(3400,4200);
	  }
	if(clu==11&&cry==2&&s ==33)
	  SetRange(4500,5500);
	if(clu==11&&cry==1&&s ==25)
	  frange = 20;
	if(clu==10){
	  SetRange(3200,4800);
	  frange = 50;
	  if(cry==0){
	    frange = 30;
	    if(s==33)
	      SetRange(4800,6200);
	    if(s==34)
	      SetRange(4400,5800);
	  }
	  if(cry==1){
	    frange = 25;
	    //if(s==36||s==33)
	    // SetRange(5000,6000);
	    if(s==20){
	      SetRange(6000,8000);
	      resolution = 10;
	    }
	    if(s==6)
	      SetRange(5000,6000);
	  }
	  if(cry==2){
	    frange = 20;
	    if(s==12) 
	      SetRange(3000,4000);
	    if(s==25)
	      SetRange(4000,5000);
	    if(s==28)
	      SetRange(3000,4000);
	    if(s==35){
	      frange = 15;
	      SetRange(3000,4000);
	    }
	  }
	  if(cry==3)
	    frange = 80;
	  if(s%10==9)
	    continue;
	}
	vector<double> r = fitCo(h,0);
	cout << " fitted " << r[0]<< endl;
	if(r[0]<0)
	  continue;
	gain.push_back(r[0]);
	offs.push_back(r[1]);
	res0.push_back(r[2]);
	res1.push_back(r[3]);
	//cout << "r[2]=" << r[2] << " r[3]=" << r[3] << endl;
	if(r[2]<0.7) cout << "///////////////////////!!!! clu=" << clu << " cry=" << cry << " s=" << s << endl;
	if(r[2]>3 && r[3] <2)
	  return;
	//if(r[2]>3)
	//  return;
	//	if(r[2]<0.8)
	//	  return;
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Gain",clu,cry,s),r[0]);
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Offset",clu,cry,s),r[1]);// /
      }//segments
    }//crystals
  }//clusters
//*/
  cout << "finished segments" << endl;
  ca->cd(3);
  cout << gain.size() << endl;
  TGraph* gs = new TGraph(gain.size(),&gain[0],&offs[0]);
  gs->Draw("AP*");
  
  ca->cd(4);
  TGraph* grs = new TGraph(res0.size(),&res0[0],&res1[0]);
  grs->Draw("AP*");
  
  cf->SaveLevel(kEnvLocal);
}
void CalibGeCoBa(){
  frange = 15;
  TEnv *cf = new TEnv(Form("%s.cal",fileBa));
  TFile *fco = new TFile(fileCo);
  TFile *fba = new TFile(fileBa);
  vector <double> gain;
  vector <double> offs;
  vector <double> res0;
  vector <double> res1;
  
  ca = new TCanvas("ca","ca",900,400);
  ca->Divide(2,1);

  TH2F* h21;
  TH2F* h22;
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      //cout << clu << "\t" << cry << endl;
      // cout << "clu=" << clu << " cry=" << cry << endl;
      h21 = (TH2F*)fco->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",clu,cry));
      h22 = (TH2F*)fba->Get(Form("hraw_segen_vs_nr_clus%02d_crys%02d",clu,cry));
      
      /*ca->cd(1);
      h21->Draw("colz");
      ca->cd(2);
      h22->Draw("colz");
      ca->WaitPrimitive();//*/

      if(h21==NULL || h22==NULL)
	continue;
      for(int s=0;s<40;s++){
	resolution = 2;
	cout << "start s=" << s << endl;
	TH1F* h11 = (TH1F*)h21->ProjectionY(Form("%s_%02d",h21->GetName(),s),s+1,s+1);
	TH1F* h12 = (TH1F*)h22->ProjectionY(Form("%s_%02d",Form("%s_2",h22->GetName()),s),s+1,s+1);
	if(h11==NULL || h11->Integral()<10 || h12==NULL || h12->Integral()<10)
	  continue;
	frange = 20;
	SetRange(1800,3100);
	//	if(clu==0 && cry==0 && s==0) SetRange(2200,2800);
	if(clu==2 && cry==1 && s==1)
	  SetRange(2000,2500);
	/*if(clu==3 && cry==2 && s==4)
	  SetRange(2800,3800);
	if(clu==4 && cry==1 && s==2)
	  SetRange(2800,3800);
	if(clu==1 && cry==1 && s==2)
       	  frange = 20;
	if(clu==3 && cry==0 && s==3)
	frange = 20;//*/
	if(clu==8 && cry==2 && s==2) continue;
	if(clu>=6 && clu<=9) SetRange(2000,3000);
	if(clu==11)
	  {
	    SetRange(3000,4200);
	    frange=25;
	    if(cry==0&&s==33) SetRange(3400,4200);
	  }
	if(clu==11&&cry==2&&s ==33)
	  SetRange(4500,5500);
	if(clu==11&&cry==1&&s ==25)
	  frange = 20;
	if(clu==10){
	  SetRange(3200,4800);
	  frange = 50;
	  if(cry==0){
	    frange = 30;
	    if(s==33)
	      SetRange(4800,6200);
	    if(s==34)
	      SetRange(4400,5800);
	  }
	  if(cry==1){
	    frange = 25;
	    //if(s==36||s==33)
	    // SetRange(5000,6000);
	    if(s==20){
	      SetRange(6000,8000);
	      resolution = 10;
	    }
	    if(s==6)
	      SetRange(5000,6000);
	  }
	  if(cry==2){
	    frange = 20;
	    if(s==12) 
	      SetRange(3000,4000);
	    if(s==25)
	      SetRange(4000,5000);
	    if(s==28)
	      SetRange(3000,4000);
	    if(s==35){
	      frange = 15;
	      SetRange(3000,4000);
	    }
	  }
	  if(cry==3)
	    frange = 80;
	  if(s%10==9)
	    continue;
	}
	vector<vector<double> > CoPar;
	cout << "Let's start some fits!" << endl;
	vector<double> vco = fitCo2(h11,CoPar,0,0);
	cout << "out" << endl;
	frange=40;
	vector<double> r=fitCoBa(h12,vco[0],CoPar,1,0);
	cout << " fitted " << r[0]<< endl;
	if(r[0]<0)
	  continue;
	cout << "r.size()=" << r.size() << endl;
	gain.push_back(r[0]);
	offs.push_back(r[1]);
	res0.push_back(r[2]);
	res1.push_back(r[3]);
	//cout << "r[2]=" << r[2] << " r[3]=" << r[3] << endl;
	if(r[2]<0.7) cout << "///////////////////////!!!! clu=" << clu << " cry=" << cry << " s=" << s << endl;
	if(r[2]>3 && r[3] <2)
	  return;
	//if(r[2]>3)
	//  return;
	//	if(r[2]<0.8)
	//	  return;
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Gain",clu,cry,s),r[0]);
	cf->SetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Offset",clu,cry,s),r[1]);//*/
	cout << "grind through" << endl;
      }//segments
    }//crystals
  }//clusters
  cout << "finished segments" << endl;
  ca->cd(1);
  cout << gain.size() << endl;
  TGraph* gs = new TGraph(gain.size(),&gain[0],&offs[0]);
  gs->Draw("AP*");
  
  ca->cd(2);
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
void CalibCore(){
  string abc[4] = {"A","B","C","D"};
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  // fileCo = (char*)"./hist/hraw0475.root";
  // fileEu = (char*)"./hist/hraw0448.root";
  // fileBa = (char*)"./hist/hraw0468.root";
  // fileYy = (char*)"./hist/hraw0469.root";
  // TEnv *cf = new TEnv("settings/all_cal0410.dat");
  fileCo = (char*)"./rootfiles/mode3_0668_0670.root";
  fileEu = (char*)"./rootfiles/mode3_0667.root";
  fileBa = (char*)"./rootfiles/mode3_0673_0674.root";
  fileYy = (char*)"./hist/hrawXXXX.root";
  TEnv *cf = new TEnv("settings/cores_20210310.dat");
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      ca->cd(cry*12+1+clu);
      TGraph *g = core(clu,cry);
      if(!g)
	continue;
      g->SetTitle(Form("DET %d%s Core calibration",clu,(char*)abc[cry].c_str()));

      g->Draw("AP*");
      TF1 *fl = new TF1("fl",flinear,0,500000,2);
      fl->SetParameters(0.004,1);
      g->Fit(fl,"q");
      double gain = fl->GetParameter(0);
      double offs = fl->GetParameter(1);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),gain);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),offs);
      
    }
  }
  cf->SaveLevel(kEnvLocal);
}
TGraph* core(int m, int c, bool draw){
  vector<double> ene;
  vector<double> chn;
  TFile *f = new TFile(fileCo);
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  if(h==NULL || h->Integral(h->FindBin(100000),h->FindBin(1000000))<100)
    return NULL;
  
  //h->Draw();
  vector<double> v = fitCo(h,1,0);
  const int nEu = 12;
  double enEu[nEu] = {121.77, 244.66, 344.28, 411.11, 443.96, 778.90, 867.3, 963.38, 1085.84, 1112.08, 1299.14, 1408.01};
  const int nCo = 3;
  double enCo[nCo] = {1173.23, 1332.49, 2505.69};
  const int nBa = 5;
  double enBa[nBa] = {81.00, 276.40, 302.85, 356.01, 383.85};
  const int nYy = 2;
  double enYy[nYy] = {898.04, 1836.06};

  frange = 30;

  f = new TFile(fileEu);
  h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nEu;i++){
    double fr = fitonepeak(h,(enEu[i]-20)/v[0] ,(enEu[i]+20)/v[0],0);
    if(fr>0){
      ene.push_back(enEu[i]);
      chn.push_back(fr);
    }
  }
  f = new TFile(fileCo);
  h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nCo;i++){
    double fr = fitonepeak(h,(enCo[i]-20)/v[0] ,(enCo[i]+20)/v[0],0);
    if(fr>0){
      ene.push_back(enCo[i]);
      chn.push_back(fr);
    }
  }
  f = new TFile(fileBa);
  h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  for(int i=0;i<nBa;i++){
    double fr = fitonepeak(h,(enBa[i]-20)/v[0] ,(enBa[i]+20)/v[0],0);
    if(fr>0){
      ene.push_back(enBa[i]);
      chn.push_back(fr);
    }
  }
  f = new TFile(fileYy);
  if(f->IsOpen()){
    h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
    for(int i=0;i<nYy;i++){
      double fr = fitonepeak(h,(enYy[i]-20)/v[0] ,(enYy[i]+20)/v[0],0);
      if(fr>0){
	ene.push_back(enYy[i]);
	chn.push_back(fr);
      }
    }
  }
  TGraph* g = new TGraph(ene.size(),&chn[0],&ene[0]);
  TF1 *fl = new TF1("fl",flinear,0,500000,2);
  fl->SetParameters(0.004,1);
  if(draw){
    ca = new TCanvas("ca","ca",1200,800);
    ca->Divide(1,2);
    ca->cd(1);
    g->Draw("AP*");
    g->Fit(fl);
    double gain = fl->GetParameter(0);
    double offs = fl->GetParameter(1);
    ca->cd(2);
    
    vector<double> res;
    fout.open("python/data/calcurve.dat");
    for(UShort_t i=0; i<ene.size(); i++){
      res.push_back(chn[i]*gain+offs - ene[i]);
      fout << ene[i] << "\t" << chn[i] << "\t" << res.back() << endl;
    }
    fout.close();
    TGraph *gr = new TGraph(ene.size(),&ene[0],&res[0]);
    gr->Draw("AP*");
  }
  else
    g->Fit(fl,"qn");

  f->Close();
  

  return g;

}
void CalibCore(int run){
  TFile *f = new TFile(Form("hist/hraw%04d.root",run));
  if(!f->IsOpen()){
    return;
  }
  ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(12,4);
  TEnv *cf = new TEnv(Form("settings/runbyrun/all_run%04d_cal0411.dat",run));
  for(int clu=0;clu<12;clu++){
    for(int cry=0;cry<4;cry++){
      ca->cd(cry*12+1+clu);
      TGraph *g = corerun(clu,cry,run);
      if(!g)
	continue;
      g->Draw("AP*");
      TF1 *fl = new TF1("fl",flinear,0,500000,2);
      fl->SetParameters(0.004,1);
      g->Fit(fl,"q");
      double gain = fl->GetParameter(0);
      double offs = fl->GetParameter(1);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",clu,cry),gain);
      cf->SetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",clu,cry),offs);
      
    }
  }
  cf->SaveLevel(kEnvLocal);
}
void allruns(){
  for(int r=447; r<483;r++){
    if(r==470)
      continue;
    CalibCore(r);
  }
}
TGraph* corerun(int m, int c, int run, bool draw){
  vector<double> ene;
  vector<double> chn;
  TFile *f = new TFile(Form("hist/hraw%04d.root",run));
  TH1F* h = (TH1F*)f->Get(Form("hraw_en_clus%02d_crys%02d",m,c));
  if(h==NULL || h->Integral(h->FindBin(100000),h->FindBin(1000000))<100)
    return NULL;

  TEnv *prev;
  if(run<476)
    prev = new TEnv("./settings/all_cal0410.dat");
  else
    prev = new TEnv("./settings/all2_cal0410.dat");
  double rough = prev->GetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",m,c),0.0);

  const int nEu = 12;
  double enEu[nEu] = {121.77, 244.66, 344.28, 411.11, 443.96, 778.90, 867.3, 963.38, 1085.84, 1112.08, 1299.14, 1408.01};
  const int nCo = 3;
  double enCo[nCo] = {1173.23, 1332.49, 2505.69};
  const int nBa = 5;
  double enBa[nBa] = {81.00, 276.40, 302.85, 356.01, 383.85};
  const int nYy = 2;
  double enYy[nYy] = {898.04, 1836.06};

  frange = 30;
  if((run>447 && run<458) || run==474 || run==479 || run ==482){
    for(int i=0;i<nEu;i++){
      double g = fitonepeak(h,(enEu[i]-20)/rough ,(enEu[i]+20)/rough,0);
      if(g>0){
	ene.push_back(enEu[i]);
	chn.push_back(c);
      }
    }
  }
  if(run==475 || run==478){
    for(int i=0;i<nCo;i++){
      double g = fitonepeak(h,(enCo[i]-20)/rough ,(enCo[i]+20)/rough,0);
      if(g>0){
	ene.push_back(enCo[i]);
	chn.push_back(g);
      }
    }
  }
  if((run>457 && run<469) || run==480 || run==481){
    for(int i=0;i<nBa;i++){
      double g = fitonepeak(h,(enBa[i]-20)/rough ,(enBa[i]+20)/rough,0);
      if(g>0){
	ene.push_back(enBa[i]);
	chn.push_back(g);
      }
    }
  }
  if(run==469){
    for(int i=0;i<nYy;i++){
      double g = fitonepeak(h,(enYy[i]-20)/rough ,(enYy[i]+20)/rough,0);
      if(g>0){
	ene.push_back(enYy[i]);
	chn.push_back(g);
      }
    }
  }
  TGraph* g = new TGraph(ene.size(),&chn[0],&ene[0]);
  TF1 *fl = new TF1("fl",flinear,0,500000,2);
  fl->SetParameters(0.004,1);
  if(draw){
    ca = new TCanvas("ca","ca",1200,800);
    ca->Divide(1,2);
    ca->cd(1);
    g->Draw("AP*");
    g->Fit(fl);
    double gain = fl->GetParameter(0);
    double offs = fl->GetParameter(1);
    ca->cd(2);
    
    vector<double> res;
    for(UShort_t i=0; i<ene.size(); i++){
      res.push_back(chn[i]*gain+offs - ene[i]);
    }
    TGraph *gr = new TGraph(ene.size(),&ene[0],&res[0]);
    gr->Draw("AP*");
  }
  else
    g->Fit(fl,"qn");

  f->Close();
  return g;

}

double fitonepeak(TH1F* h, double low, double hig, bool draw){
  h->GetXaxis()->SetRangeUser(low,hig);
  if(h->Integral() < 10)
    return -1;
  TSpectrum *sp = new TSpectrum(3,resolution);
  sp->SetResolution(resolution);
  Int_t nfound = 0;
  if(draw)
    nfound = sp->Search(h,resolution,"nobackground",0.5);
  else
    nfound = sp->Search(h,resolution,"nobackgroundgoff",0.5);

  if(nfound!=1){
    cout << "Found " << nfound << " peaks in spectrum, not one, try again" << endl;
    TH1F* hc = (TH1F*)h->Clone("testclone");
    hc->Rebin(2);
    hc->GetXaxis()->SetRangeUser(low,hig);
    nfound = 0;
    if(draw)
      nfound = sp->Search(hc,resolution,"nobackground",0.4);
    else
      nfound = sp->Search(hc,resolution,"nobackgroundgoff",0.4);
      
    if(nfound!=1){
      cout << "Found " << nfound << " peaks in spectrum, not one, aborting" << endl;
      if(draw)
	hc->DrawCopy();
      return -1;
    }      
  }
  Double_t xpeak = sp->GetPositionX()[0];
  Double_t ypeak = sp->GetPositionY()[0];
  if(draw){
    cout << xpeak << "\t" << ypeak << endl;
    ca = new TCanvas("ca","ca",600,600);
    ca->cd();
    h->DrawCopy();
  }
  TF1 *fu;
  TF1 *fus[3];
  h->GetXaxis()->SetRangeUser(xpeak-frange,xpeak+frange);
  fu = new TF1(Form("f%s",h->GetName()),fgammagaussbg,xpeak-frange,xpeak+frange,6);
  fu->SetLineColor(3);
  fu->SetLineWidth(1);
  fu->SetParameter(0,0);//bg const
  fu->SetParameter(1,0);//bg slope
  fu->SetParameter(2,h->Integral(xpeak-frange,xpeak+frange));//norm
  fu->SetParameter(3,xpeak);//mean
  fu->SetParLimits(3,xpeak-500,xpeak+500);//mean
  fu->SetParameter(4,200);//sigma
  fu->SetParLimits(4,100,frange);//sigma
  fu->SetParameter(5,h->GetBinContent(h->FindBin(xpeak-frange)));//step
  if(draw)
    h->Fit(fu,"Rn");
  else
    h->Fit(fu,"Rqn");

  //draw results.
  if(draw){
    fu->Draw("same");
    fus[0] = new TF1(Form("f%s_bg",h->GetName()),fgammabg,xpeak-frange,xpeak+frange,6);
    fus[1] = new TF1(Form("f%s_st",h->GetName()),fgammastep,xpeak-frange,xpeak+frange,6);
    fus[2] = new TF1(Form("f%s_ga",h->GetName()),fgammagaus,xpeak-frange,xpeak+frange,6);
	

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
  return fu->GetParameter(3);
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

Double_t fgammagausnostep(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t result = par[0] + par[1]*x[0]; 

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  return result;

}

Double_t flinear(Double_t *x, Double_t *par){
  return x[0]*par[0] + par[1];

}
