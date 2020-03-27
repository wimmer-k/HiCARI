#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "TMath.h"
#include "TEnv.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;
double deg2rad = TMath::Pi()/180.;
double rad2deg = 180./TMath::Pi();
void Positions(){
  ifstream infile;
  infile.open("/home/gamma20/HiCARI/settings/MBpositions_jiseok_mar26.dat");
  infile.ignore(1000,'\n');
  int clu, cry, seg;
  float phi, the, rho;
  TH2F* c_thetaphi = new TH2F("c_thetaphi", "c_thetaphi",180,0,180,180,-180,180);
  TH2F* s_thetaphi = new TH2F("s_thetaphi", "s_thetaphi",180,0,180,180,-180,180);
  TH2F* c_xy = new TH2F("c_xy", "c_xy",200,-200,200,200,-200,200);
  TH2F* s_xy = new TH2F("s_xy", "s_xy",200,-200,200,200,-200,200);
  TEnv* fout = new TEnv("/home/gamma20/HiCARI/settings/HiCARIpos_0327.dat");
  while(!infile.eof()){
  //for(int i=0;i<10;i++){
    infile >> clu >> cry >> seg >> phi >> the >> rho;
    
    infile.ignore(1000,'\n');
    // if(clu!=1)
    //   continue;
    cout << clu << "\t" << cry << "\t" << seg << "\t" <<phi << "\t" << the << "\t" << rho << endl;
    seg -= 1;
    TVector3 v(0,0,1);
    v.SetTheta(the*deg2rad);
    v.SetPhi(-phi*deg2rad);
    v.SetMag(rho);
    
    if(seg==-1){
      c_thetaphi->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
      c_xy->Fill(v.X(),v.Y());
    }
    else{
      s_thetaphi->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
      s_xy->Fill(v.X(),v.Y());
    }
    if(seg>-1){
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.X",clu,cry,seg),v.X());
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Y",clu,cry,seg),v.Y());
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Z",clu,cry,seg),v.Z());
    }
    if(infile.eof())
      break;
  }
  TCanvas* c = new TCanvas("c","c",800,400);
  c->Divide(2,1);
  c->cd(1);
  c_thetaphi->SetMarkerStyle(20);
  c_thetaphi->SetMarkerSize(0.5);
  c_thetaphi->Draw();
  s_thetaphi->SetMarkerColor(3);
  s_thetaphi->SetMarkerStyle(20);
  s_thetaphi->SetMarkerSize(0.5);
  s_thetaphi->Draw("same");
  
  c->cd(2);
  c_xy->SetMarkerStyle(20);
  c_xy->SetMarkerSize(0.5);
  c_xy->Draw();
  s_xy->SetMarkerColor(3);
  s_xy->SetMarkerStyle(20);
  s_xy->SetMarkerSize(0.5);
  s_xy->Draw("same");
  
  fout->SaveLevel(kEnvLocal);

}
