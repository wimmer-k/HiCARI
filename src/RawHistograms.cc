#include "RawHistograms.hh"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include <sstream>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

using namespace TMath;
using namespace std;

void RawHistograms::Write(){
  fhlist->Sort();
  fhlist->Write();

}

void RawHistograms::FillHistograms(Mode3Event* m3e, HiCARI* hi, Gretina* gr){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hasmode3 = m3e->GetMult()!=0;
  bool hashicari = hi->GetMult()!=0;
  bool hasgretina = gr != NULL && gr->GetMult()!=0;


  if(fSett->Mode3Histos() && hasmode3){
    FillMode3Histograms(m3e);
  }
  if(hashicari){
    FillHiCARIHistograms(hi);
  }
  if(hasgretina){
    FillGretinaHistograms(gr);
  }
}
void RawHistograms::FillMode3Histograms(Mode3Event* m3e){
  //Mode 2 histograms
  //static bool firstEvent = true;

  Fill("hmode3_emult",30,0,30,m3e->GetMult());
  for(int i=0; i<m3e->GetMult(); i++){
    Mode3Hit* hit = m3e->GetHit(i);
    Fill("hmode3_hmult",30,0,30,hit->GetMult());
    for(int j=0; j<hit->GetMult(); j++){
      Trace * trace = hit->GetTrace(j);
      if(fSett->VLevel()>1){
	cout << "Trace " << j << " Length " << trace->GetLength() << 
	  "\tEnergy " << trace->GetEnergy() <<
	  "\tBoard " << trace->GetBoard() <<
	  "\tSlot " << trace->GetSlot() <<
	  "\tChannel " << trace->GetChn() <<
	  "\tHole " << trace->GetHole() <<
	  "\tCrystal " << trace->GetCrystal();
	if(trace->GetChn()==9)
	  cout << " <- CC " << endl;
	else if(trace->GetEnergy()>1000)
	  cout << " <- with net energy " << endl;
	else
	  cout << endl;
      }
      Fill("hmode3_hole",20,0,20,trace->GetHole());
      Fill(Form("hmode3_crys_hole%02d",trace->GetHole()),10,0,10,trace->GetCrystal());
      Fill(Form("hmode3_slot_hole%02d_crys%02d",trace->GetHole(),trace->GetCrystal()),10,0,10,trace->GetSlot());
      Fill(Form("hmode3_chan_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),10,0,10,trace->GetChn());
      Fill(Form("hmode3_en_hole%02d_crys%02d_slot%02d_chan%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot(),trace->GetChn()),20000,0,1e6,trace->GetEnergy());
      Fill(Form("hmode3_en_vs_chn_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),10,0,10,trace->GetChn(),1000,0,1e6,trace->GetEnergy());
      //for MB only, tracking have several ch9 and CL have the cores on 0 and 5....
      if(trace->GetChn()==9){
	Fill(Form("hmode3_core_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),2000,0,1e6,trace->GetEnergy());
	if(fSett->TracePlots()){
	  double baseline =0;
	  for(int j=0;j<trace->GetLength();j++){
	    Fill(Form("hmode3_coretraces_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),200,0,200,j,1200,-10000,2000,(int)trace->GetTrace()[j]);
	    if(j<fSett->BaselineLength())
	      baseline+=trace->GetTrace()[j];
	  }
	  baseline/=fSett->BaselineLength();
	  Fill(Form("hmode3_smoke_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),1000,0,1e6,trace->GetEnergy(),1000,0,1000,baseline);
	}
      }
    }
  }
}

void RawHistograms::FillHiCARIHistograms(HiCARI* hi){
  Fill("hraw_mult",30,0,30,hi->GetMult());
  for(int i=0; i<hi->GetMult(); i++){
    int segs = 6;
    HiCARIHit* hit = hi->GetHit(i);
    if(hit->IsTracking())
      segs = 40;
    if(hit->IsSuperClo())
      segs = 4;
    Fill("hraw_cluster",12,0,12,hit->GetCluster());
    Fill("hraw_crystal",4,0,4,hit->GetCrystal());
    Fill("hraw_crystal_vs_cluster",12,0,12,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("hraw_en",5000,0,3e4,hit->GetEnergy());
    Fill("hraw_en_LR",5000,0,3e5,hit->GetEnergy());
    Fill("hraw_en_summary",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),5000,0,3e4,hit->GetEnergy());
    Fill("hraw_en_summary_LR",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),5000,0,3e5,hit->GetEnergy());
    //temp increase spectrum range gain seems larger for P3 pos 2
    if(hit->GetCluster()==11 && hit->GetCrystal()==1){
      Fill(Form("hraw_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1e4,0,3e4,hit->GetEnergy());
    }
    else{
      Fill(Form("hraw_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1e4,0,1e4,hit->GetEnergy());
    }
    //Fill(Form("hraw_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),5000,0,1e6,hit->GetEnergy());
    Fill(Form("hraw_segsum_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1000,0,1e4,hit->GetEnergy(),1000,0,1e4,hit->GetSegmentSum());
    Fill(Form("hraw_segmult_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetMult());
    
    for(int j=0; j<hit->GetMult(); j++){
      Fill(Form("hraw_segen_vs_nr_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr(j),5000,0,1e4,hit->GetSegmentEn(j));
      Fill(Form("hraw_segen_vs_nr_LR_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr(j),5000,0,1e5,hit->GetSegmentEn(j));
      Fill("hraw_segen",5000,0,1e4,hit->GetSegmentEn(j));
      Fill("hraw_segen_LR",5000,0,1e5,hit->GetSegmentEn(j));
    }
  }
}
void RawHistograms::FillGretinaHistograms(Gretina* gr){
  Fill("graw_mult",30,0,30,gr->GetMult());
}
