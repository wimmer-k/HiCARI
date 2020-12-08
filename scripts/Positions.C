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
void MBPositions(int show=-1){
  ifstream infile;
  //infile.open("/home/wimmer/programs/HiCARI/settings/MBpositions_jiseok_mar26.dat");
  infile.open("/home/gamma20/packages/HiCARI/settings/MB_geometery_fall.txt");
  infile.ignore(1000,'\n');
  infile.ignore(1000,'\n');
  infile.ignore(1000,'\n');
  int clu, cry, seg;
  float phi, the, rho;
  TH2F* c_thetaphi = new TH2F("c_thetaphi", "c_thetaphi",180,0,180,180,-180,180);
  TH2F* s_thetaphi = new TH2F("s_thetaphi", "s_thetaphi",180,0,180,180,-180,180);
  TH2F* c_xy = new TH2F("c_xy", "c_xy",200,-400,400,200,-400,400);
  TH2F* s_xy = new TH2F("s_xy", "s_xy",200,-400,400,200,-400,400);
  TEnv* fout = new TEnv("/home/gamma20/packages/HiCARI/settings/HiCARIpos_1022.dat");
  while(!infile.eof()){
  //for(int i=0;i<10;i++){
    infile >> clu >> cry >> seg >> phi >> the >> rho;
    
    infile.ignore(1000,'\n');
    if(show > -1 && clu!=show)
      continue;
    cout << clu << "\t" << cry << "\t" << seg << "\t" <<phi << "\t" << the << "\t" << rho << endl;
    seg -= 1;
    TVector3 v(0,0,1);
    v.SetTheta(the*deg2rad);
    v.SetPhi(phi*deg2rad);
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
void SCPositions(int show=-1){
  ifstream infile;
  infile.open("/home/gamma20/packages/HiCARI/settings/SC_seg_fall_xyz.txt");
  infile.ignore(1000,'\n');
  infile.ignore(1000,'\n');
  int clu, cry, seg;
  float x, y, z;
  //TH2F* c_thetaphi = new TH2F("c_thetaphi", "c_thetaphi",180,0,180,180,-180,180);
  TH2F* s_thetaphi = new TH2F("s_thetaphi", "s_thetaphi",180,0,180,180,-180,180);
  //TH2F* c_xy = new TH2F("c_xy", "c_xy",200,-400,400,200,-400,400);
  TH2F* s_xy = new TH2F("s_xy", "s_xy",200,-400,400,200,-400,400);
  TEnv* fout = new TEnv("/home/gamma20/packages/HiCARI/settings/HiCARIpos_1022.dat");
  while(!infile.eof()){
  //for(int i=0;i<10;i++){
    if(infile.eof())
      break;
    infile >> clu >> cry >> seg >> x >> y >> z;
    
    infile.ignore(1000,'\n');
    if(infile.eof())
      break;
    // if(clu!=1)
    //   continue;
    if(show > -1 && clu!=show)
      continue;
    cout << clu << "\t" << cry << "\t" << seg << "\t" << x << "\t" << y << "\t" << z << endl;
    seg -= 1;
    TVector3 v(x,y,z);
    
    // if(seg==-1){
    //   c_thetaphi->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
    //   c_xy->Fill(v.X(),v.Y());
    // }
    // else{
    //   s_thetaphi->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
    //   s_xy->Fill(v.X(),v.Y());
    // }
    s_thetaphi->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
    s_xy->Fill(v.X(),v.Y());
    if(seg>-1){
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.X",clu+6,cry,seg),v.X());
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Y",clu+6,cry,seg),v.Y());
      fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Z",clu+6,cry,seg),v.Z());
    }
    if(infile.eof())
      break;
  }
  TCanvas* c = new TCanvas("c","c",800,400);
  c->Divide(2,1);
  c->cd(1);
  //c_thetaphi->SetMarkerStyle(20);
  //c_thetaphi->SetMarkerSize(0.5);
  //c_thetaphi->Draw();
  s_thetaphi->SetMarkerColor(3);
  s_thetaphi->SetMarkerStyle(20);
  s_thetaphi->SetMarkerSize(0.5);
  //s_thetaphi->Draw("same");
  s_thetaphi->Draw();
  
  c->cd(2);
  //c_xy->SetMarkerStyle(20);
  //c_xy->SetMarkerSize(0.5);
  //c_xy->Draw();
  s_xy->SetMarkerColor(3);
  s_xy->SetMarkerStyle(20);
  s_xy->SetMarkerSize(0.5);
  //s_xy->Draw("same");
  s_xy->Draw();
  
  fout->SaveLevel(kEnvLocal);

}
void comparewithsim(char* dat = "settings/HiCARIpos_0328.dat", char* sim = "simulation/settings/coordinates_sim.dat"){
  //TEnv* env = new TEnv("settings/139Sn_MBcoordinates.dat");
  TEnv* env[2];
  env[0] = new TEnv(dat);
  env[1] = new TEnv(sim);
  TVector3 pdat[10][4][8];
  TVector3 psim[10][4][8];
  TVector3 pdatc[10][4];
  TVector3 psimc[10][4];
  TVector3 pdatf[10];
  TVector3 psimf[10];
  TH2F* h[2][4];
  h[0][0] = new TH2F("hdat_tp","hdat_tp",100,0,100,180,-180,180);
  h[1][0] = new TH2F("hdat_tp","hdat_tp",100,0,100,180,-180,180);
  h[0][1] = new TH2F("hdat_xy","hdat_xy",200,-300,300,200,-300,300);
  h[1][1] = new TH2F("hdat_xy","hdat_xy",200,-300,300,200,-300,300);
  h[0][2] = new TH2F("hdat_xz","hdat_xz",200,-300,300,200,0,300);
  h[1][2] = new TH2F("hdat_xz","hdat_xz",200,-300,300,200,0,300);
  h[0][3] = new TH2F("hdat_yz","hdat_yz",200,-300,300,200,0,300);
  h[1][3] = new TH2F("hdat_yz","hdat_yz",200,-300,300,200,0,300);
  ofstream fout;
  fout.open("python/data/coordinates0505.dat");
  for(int cl=0; cl<10; cl++){
    int fdat =0;
    int fsim =0;
    pdatf[cl].SetXYZ(0,0,0);
    psimf[cl].SetXYZ(0,0,0);
    for(int cr=0; cr<4; cr++){
      int cdat =0;
      int csim =0;
      pdatc[cl][cr].SetXYZ(0,0,0);
      psimc[cl][cr].SetXYZ(0,0,0);
      for(int se=0; se<8; se++){
	double x,y,z;
	x = env[0]->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.X",cl,cr,se),-0.0);
	y = env[0]->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Y",cl,cr,se),-0.0);
	z = env[0]->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Z",cl,cr,se),-0.0);
	pdat[cl][cr][se].SetXYZ(x,y,z);
	if(pdat[cl][cr][se].Theta()>0){
	  h[0][0]->Fill(pdat[cl][cr][se].Theta()*180./TMath::Pi(),pdat[cl][cr][se].Phi()*180./TMath::Pi());
	  h[0][1]->Fill(pdat[cl][cr][se].X(),pdat[cl][cr][se].Y());
	  h[0][2]->Fill(pdat[cl][cr][se].X(),pdat[cl][cr][se].Z());
	  h[0][3]->Fill(pdat[cl][cr][se].Y(),pdat[cl][cr][se].Z());
	  fout << cl << "\t" << cr << "\t" << se << "\t" << pdat[cl][cr][se].Theta()*180./TMath::Pi() << "\t" << pdat[cl][cr][se].Phi()*180./TMath::Pi() << "\t" << x << "\t" << y;
	  pdatc[cl][cr] += pdat[cl][cr][se];
	  cdat++;
	  pdatf[cl] += pdat[cl][cr][se];
	  fdat++;
	}
	x = env[1]->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.X",cl,cr,se),-0.0);
	y = env[1]->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Y",cl,cr,se),-0.0);
	z = env[1]->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Z",cl,cr,se),-0.0);
	//cout << x << "\t" << y << "\t" << z << endl;
	psim[cl][cr][se].SetXYZ(x,y,z);
	if(psim[cl][cr][se].Theta()>0){
	  h[1][0]->Fill(psim[cl][cr][se].Theta()*180./TMath::Pi(),psim[cl][cr][se].Phi()*180./TMath::Pi());
	  h[1][1]->Fill(psim[cl][cr][se].X(),psim[cl][cr][se].Y());
	  h[1][2]->Fill(psim[cl][cr][se].X(),psim[cl][cr][se].Z());
	  h[1][3]->Fill(psim[cl][cr][se].Y(),psim[cl][cr][se].Z());
	  fout << "\t" << psim[cl][cr][se].Theta()*180./TMath::Pi() << "\t" << psim[cl][cr][se].Phi()*180./TMath::Pi() << "\t" << x << "\t" << y << endl;
	  psimc[cl][cr] += psim[cl][cr][se];
	  csim++;
	  psimf[cl] += psim[cl][cr][se];
	  fsim++;
	}
	// if(pdat[cl][cr][se].Theta()>0 && psim[cl][cr][se].Theta()>0){
	//   cout << cl << "\t" << cr << "\t" << se << "\t" << pdat[cl][cr][se].Mag() << "\t" << psim[cl][cr][se].Mag() << "\t" << pdat[cl][cr][se].Mag()-psim[cl][cr][se].Mag() << endl;
	// }
      }//segments
      // if(cdat>0 && csim>0){
      // 	pdatc[cl][cr] *= 1./cdat;
      // 	psimc[cl][cr] *= 1./csim;
      // 	cout << cl << "\t" << cr << "\t" << pdatc[cl][cr].Mag() << "\t" << psimc[cl][cr].Mag() << "\t" << pdatc[cl][cr].Mag()-psimc[cl][cr].Mag() << endl;
      // }
    }//crystals
    
    if(fdat>0 && fsim>0){
      pdatf[cl] *= 1./fdat;
      psimf[cl] *= 1./fsim;
      cout << cl << "\t" << pdatf[cl].Mag() << "\t" << psimf[cl].Mag() << "\t" << pdatf[cl].Mag()-psimf[cl].Mag() << endl;
    }
    
  }
  TCanvas *c = new TCanvas("c","c",350*4,300*4);
  c->Divide(2,2);
  for(int i=0;i<4;i++){
    c->cd(1+i);    
    h[0][i]->Draw();
    h[0][i]->SetMarkerColor(3);
    h[0][i]->SetMarkerSize(0.5);
    h[0][i]->SetMarkerStyle(20);
    h[1][i]->Draw("same");
    h[1][i]->SetMarkerColor(2);
    h[1][i]->SetMarkerSize(0.5);
    h[1][i]->SetMarkerStyle(20);
  }
  h[0][0]->GetYaxis()->SetTitle("#phi_{lab} (deg)");
  h[0][0]->GetXaxis()->SetTitle("#theta_{lab} (deg)");

}





void ReadMatrix(const char* filename);
TVector3 TransformCoordinates(int hole, int cry, TVector3 local);
float fcrmat[20][4][4][4];


void QuadPositions(){

  ReadMatrix("/home/gamma20/packages/HiCARI/settings/quadmatrix.dat");
  
  ifstream infile;
  int clu=15, cry=1, seg=0;
  float x, y, z;

  //TEnv* fout = new TEnv("/home/gamma20/packages/HiCARI/settings/HiCARIpos_1025.dat");
  TH2F* c_thetaphi[4];
  TH2F* c_xy[4];
  TH2F* c_xz[4];
  TH2F* c_yz[4];
  TH2F* s_thetaphi[4];
  TH2F* s_xy[4];
  TH2F* s_xz[4];
  TH2F* s_yz[4];
  TCanvas* c = new TCanvas("c","c",1200,400);
  c->Divide(2,2);
  for(cry=0;cry<4;cry++){
    seg = 0;
    s_thetaphi[cry] = new TH2F(Form("s_thetaphi_%d",cry), Form("s_thetaphi_%d",cry),180,0,180,180,-180,180);
    s_xy[cry] = new TH2F(Form("s_xy_%d",cry), Form("s_xy_%d",cry),200,-400,0,200,-200,200);
    s_xz[cry] = new TH2F(Form("s_xz_%d",cry), Form("s_xz_%d",cry),200,-400,0,200,-200,200);
    s_yz[cry] = new TH2F(Form("s_yz_%d",cry), Form("s_yz_%d",cry),200,-100,300,200,-200,200);
    c_thetaphi[cry] = new TH2F(Form("s_thetaphi_%d",cry), Form("s_thetaphi_%d",cry),180,0,180,180,-180,180);
    c_xy[cry] = new TH2F(Form("s_xy_%d",cry), Form("s_xy_%d",cry),200,-400,0,200,-200,200);
    c_xz[cry] = new TH2F(Form("s_xz_%d",cry), Form("s_xz_%d",cry),200,-400,0,200,-200,200);
    c_xz[cry] = new TH2F(Form("s_xz_%d",cry), Form("s_xz_%d",cry),200,-400,0,200,-200,200);
    c_yz[cry] = new TH2F(Form("s_yz_%d",cry), Form("s_yz_%d",cry),200,-100,300,200,-200,200);

    infile.open(Form("/home/gamma20/packages/HiCARI/settings/cogP%d.dat",cry+1));
    if(!infile.is_open()){
      cout << " not found file" << endl;
      cout << Form("/home/gamma20/packages/HiCARI/settings/cogP%d.dat",cry+1) << endl;
      continue;
    }
    
    infile.ignore(1000,'\n');
    while(!infile.eof()){
      infile >> x >> y >> z;
      infile.ignore(1000,'\n');
      cout << clu << "\t" << cry << "\t" << seg << "\t" << x << "\t" << y << "\t" << z << endl;
      TVector3 v(x,y,z);
      v = TransformCoordinates(15,cry,v);
      c_thetaphi[cry]->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
      c_xy[cry]->Fill(v.X(),v.Y());
      c_xz[cry]->Fill(v.X(),v.Z());
      c_yz[cry]->Fill(v.Y(),v.Z());
      if(seg==0){
	s_thetaphi[cry]->Fill(v.Theta()*rad2deg,v.Phi()*rad2deg);
	s_xy[cry]->Fill(v.X(),v.Y());
	s_xz[cry]->Fill(v.X(),v.Z());
 	s_yz[cry]->Fill(v.Y(),v.Z());
      }
	
      // fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.X",10,cry,seg),v.X());
      // fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Y",10,cry,seg),v.Y());
      // fout->SetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Z",10,cry,seg),v.Z());

      seg++;
      if(seg%10==9)
	seg++;
      if(infile.eof())
	break;
    }
    c->cd(1);
    c_thetaphi[cry]->SetMarkerColor(2+cry);
    c_thetaphi[cry]->SetMarkerStyle(20);
    c_thetaphi[cry]->SetMarkerSize(0.5);
    if(cry==0)
      c_thetaphi[cry]->Draw();
    else
      c_thetaphi[cry]->Draw("same");
    s_thetaphi[cry]->SetMarkerColor(2+cry);
    s_thetaphi[cry]->SetMarkerStyle(20);
    s_thetaphi[cry]->SetMarkerSize(2);
    s_thetaphi[cry]->Draw("same");

    c->cd(2);
    c_xy[cry]->SetMarkerColor(2+cry);
    c_xy[cry]->SetMarkerStyle(20);
    c_xy[cry]->SetMarkerSize(0.5);
    if(cry==0)
      c_xy[cry]->Draw();
    else
      c_xy[cry]->Draw("same");
    s_xy[cry]->SetMarkerColor(2+cry);
    s_xy[cry]->SetMarkerStyle(20);
    s_xy[cry]->SetMarkerSize(2);
    s_xy[cry]->Draw("same");

    c->cd(3);
    c_xz[cry]->SetMarkerColor(2+cry);
    c_xz[cry]->SetMarkerStyle(20);
    c_xz[cry]->SetMarkerSize(0.5);
    if(cry==0)
      c_xz[cry]->Draw();
    else
      c_xz[cry]->Draw("same");
    s_xz[cry]->SetMarkerColor(2+cry);
    s_xz[cry]->SetMarkerStyle(20);
    s_xz[cry]->SetMarkerSize(2);
    s_xz[cry]->Draw("same");
 
    c->cd(4);
    c_yz[cry]->SetMarkerColor(2+cry);
    c_yz[cry]->SetMarkerStyle(20);
    c_yz[cry]->SetMarkerSize(0.5);
    if(cry==0)
      c_yz[cry]->Draw();
    else
      c_yz[cry]->Draw("same");
    s_yz[cry]->SetMarkerColor(2+cry);
    s_yz[cry]->SetMarkerStyle(20);
    s_yz[cry]->SetMarkerSize(2);
    s_yz[cry]->Draw("same");
 
    infile.close();
  }
  
  //fout->SaveLevel(kEnvLocal);
 
}
void ReadMatrix(const char* filename){
  ifstream infile;
  infile.open(filename);
  if(!infile.is_open()){
    cout << "no matrix found" << endl;
    return;
  }
  int hole,cry;
  while(!infile.eof()){
    //int hole,cry;
    infile >> hole >> cry;
    infile.ignore(100,'\n');
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++){
	infile >> fcrmat[hole][cry][i][j];
      }
      infile.ignore(100,'\n');
    }
    if(infile.eof())
      break;
  }
  for(int hole=0;hole<20;hole++){
    for(int cry=0;cry<4;cry++){
      if(fcrmat[hole][cry][3][3] > 0){
	for(int i=0;i<4;i++){
	  for(int j=0;j<4;j++){
	    cout << fcrmat[hole][cry][i][j] << "\t";
	  }
	  cout << endl;
	}
      }
    }//crystals
  }//holes
}
TVector3 TransformCoordinates(int hole, int cry, TVector3 local){
  double x = local.X();
  double y = local.Y();
  double z = local.Z();
   double xt = fcrmat[hole][cry][0][0] * x + fcrmat[hole][cry][0][1] * y + fcrmat[hole][cry][0][2] * z + fcrmat[hole][cry][0][3];
  double yt = fcrmat[hole][cry][1][0] * x + fcrmat[hole][cry][1][1] * y + fcrmat[hole][cry][1][2] * z + fcrmat[hole][cry][1][3];
  double zt = fcrmat[hole][cry][2][0] * x + fcrmat[hole][cry][2][1] * y + fcrmat[hole][cry][2][2] * z + fcrmat[hole][cry][2][3];
  return TVector3(xt,yt,zt);
}
