#include "MergeHistograms.hh"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include <sstream>



using namespace TMath;
using namespace std;

void MergeHistograms::Write(int runbr, int runhi){
  fhlist->Sort();
  fhlist->Write();
  PrintHistos(runbr,runhi);
}

double MergeHistograms::GetCorrRate(){
  try{
    fhmap.at("hEcorr")->Integral();
  } catch(const out_of_range& e) {
    return sqrt(-1);
  }
  double r = fhmap.at("hEcorr")->Integral()/fhmap.at("hEcorr")->GetEntries() *100;
  Fill("hcorrate_entry",10000,0,1e7,fentry,100,0,100,r);
  
  return r;
}

void MergeHistograms::FillHistograms(int checkADC, HiCARICalc* hi, GretinaCalc* gr, unsigned long long int brTS, unsigned long long int hiTS, unsigned long long int m2TS){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hashicari = hi->GetMult()!=0;
  bool hasbigrips = true;//checkADC>0;

  
  FillCorrelationHistograms(checkADC, hi, brTS,hiTS,m2TS);
}
void MergeHistograms::FillCorrelationHistograms(int checkADC, HiCARICalc* hi, unsigned long long int brTS, unsigned long long int hiTS, unsigned long long int m2TS){
  Fill("hTSdiff_BR_HI",2000,-1000,1000,brTS-hiTS);
  Fill("hTSdiff_BR_M2",2000,-1000,1000,brTS-m2TS);
  Fill("hTSdiff_HI_M2",2000,-1000,1000,hiTS-m2TS);

  Fill("hentry_BR",1000,0,1e6,fentry,1000,0,5e12,brTS);
  Fill("hentry_HI",1000,0,1e6,fentry,1000,0,5e12,hiTS);
  Fill("hentry_M2",1000,0,1e6,fentry,1000,0,5e12,m2TS);
  
  for(UShort_t i=0; i<hi->GetMult();i++){
    HiCARIHitCalc* hit = hi->GetHit(i);
    if(hit->GetCluster()==fSett->CorrelationCluster() && hit->GetCrystal()==fSett->CorrelationCrystal())
      Fill("hEcorr",1024,0,4096,checkADC, 1000,0,4000,hit->GetEnergy());
  }
      
}

void MergeHistograms::PrintHistos(int runbr, int runhi){
  cout << "printing" << endl;
  TCanvas *c = new TCanvas("c","c",1200,1200);
  c->Divide(2,2);
  const int n2d = 2;
  char* names2d[n2d] = {"hcorrate_entry","hEcorr"};
  char* title2d[n2d] = {"Event Correlation Rate (%)", "Energy Correlation"};
  for(int i=0;i<n2d;i++){
    c->cd(1+i);
    try{
      fhmap.at(names2d[i])->Draw("colz");
      fhmap.at(names2d[i])->SetTitle(title2d[i]);
    } catch(const out_of_range& e) {
      cout << "histogram " << names2d[i]<< " not found " << endl;
      continue;
    }
  }
  c->cd(n2d+1);
  gPad->SetLogy();
  const int n1d = 3;
  char* names1d[n1d] = {"hTSdiff_BR_HI","hTSdiff_BR_M2","hTSdiff_HI_M2"};
  char* title1d[n1d] = {"TS difference BR - HI", "TS difference BR - mode2", "TS difference HI - mode2"};
  TLegend * l = new TLegend(0.1,0.5,0.45,0.9);
  for(int i=0;i<n1d;i++){
    try{
      if(i==0)
	fhmap.at(names1d[i])->Draw();
      fhmap.at(names1d[i])->Draw("same");	
      fhmap.at(names1d[i])->SetLineColor(1+i);
      fhmap.at(names1d[i])->SetTitle(title1d[i]);
      l->AddEntry(fhmap.at(names1d[i]),title1d[i],"L");
    } catch(const out_of_range& e) {
      cout << "histogram " << names1d[i]<< " not found " << endl;
      continue;
    }
  }
  l->Draw();
  c->cd(4);
  if(runbr>0 || runhi>0){
    cout << runbr << "\t "<< runhi << endl;
    TLatex *la = new TLatex(0.1,0.5,Form("run BR%04d HI%04d",runbr,runhi));
    la->Draw();
  }
  c->SaveAs("/home/gamma20/fall2020/plots/current_merge.png");
    
}
