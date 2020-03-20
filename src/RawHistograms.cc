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

#ifdef SIMULATION
#include "Gretina.hh"
#endif
using namespace TMath;
using namespace std;

void RawHistograms::Write(){
  fhlist->Sort();
  fhlist->Write();

}

#ifdef SIMULATION
void RawHistograms::FillHistograms(Gretina* gr, Miniball* mb, ZeroDeg* zd, MINOS* mi){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hasmode2 = gr->GetMult()!=0;
  if(hasmode2){
    FillMode2Histograms(gr);
  }
}
void RawHistograms::FillHistograms(Mode3Event* m3e, Miniball* mb, Gretina* gr){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hasmode2 = gr->GetMult()!=0;
  bool hasmode3 = m3e->GetMult()!=0;
  bool hasminib = mb->GetMult()!=0;


  if(hasmode3){
    FillMode3Histograms(m3e);
  }
  if(hasmode2){
    FillMode2Histograms(gr);
  }
  if(hasminib){
    FillMiniballHistograms(mb);
  }
}
#else
void RawHistograms::FillHistograms(Mode3Event* m3e, Germanium* ge){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hasmode3 = m3e->GetMult()!=0;
  bool hasgerma = ge->GetMult()!=0;


  if(hasmode3){
    FillMode3Histograms(m3e);
  }
  if(hasgerma){
    FillGermaniumHistograms(ge);
  }
}
#endif
void RawHistograms::FillMode3Histograms(Mode3Event* m3e){
  //Mode 2 histograms
  //static bool firstEvent = true;

  Fill("hraw_emult",30,0,30,m3e->GetMult());
  for(int i=0; i<m3e->GetMult(); i++){
    Mode3Hit* hit = m3e->GetHit(i);
    Fill("hraw_hmult",30,0,30,hit->GetMult());
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
      Fill("hraw_hole",20,0,20,trace->GetHole());
      Fill(Form("hraw_crys_hole%02d",trace->GetHole()),10,0,10,trace->GetCrystal());
      Fill(Form("hraw_slot_hole%02d_crys%02d",trace->GetHole(),trace->GetCrystal()),10,0,10,trace->GetSlot());
      Fill(Form("hraw_chan_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),10,0,10,trace->GetChn());
      Fill(Form("hraw_en_hole%02d_crys%02d_slot%02d_chan%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot(),trace->GetChn()),20000,0,1e6,trace->GetEnergy());
      Fill(Form("hraw_en_vs_chn_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),10,0,10,trace->GetChn(),1000,0,1e6,trace->GetEnergy());
      //for MB, SC only
      if(trace->GetChn()==9){
	Fill(Form("hraw_core_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),2000,0,1e6,trace->GetEnergy());
	if(fSett->TracePlots()){
	  double baseline =0;
	  for(int j=0;j<trace->GetLength();j++){
	    Fill(Form("hraw_coretraces_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),200,0,200,j,1200,-10000,2000,(int)trace->GetTrace()[j]);
	    if(j<fSett->BaselineLength())
	      baseline+=trace->GetTrace()[j];
	  }
	  baseline/=fSett->BaselineLength();
	  Fill(Form("hraw_smoke_hole%02d_crys%02d_slot%02d",trace->GetHole(),trace->GetCrystal(),trace->GetSlot()),1000,0,1e6,trace->GetEnergy(),1000,0,1000,baseline);
	}
      }
    }
  }
}

void RawHistograms::FillGermaniumHistograms(Germanium* ge){
  Fill("h_mult",30,0,30,ge->GetMult());
  for(int i=0; i<ge->GetMult(); i++){
    int segs = 6;
    GeCrystal* cr = ge->GetHit(i);
    if(cr->IsTracking())
      segs = 40;
    Fill("h_cluster",12,0,12,cr->GetCluster());
    Fill("h_crystal",4,0,4,cr->GetCrystal());
    Fill("h_crystal_vs_cluster",12,0,12,cr->GetCluster(),4,0,4,cr->GetCrystal());
    Fill(Form("h_en_clus%02d_crys%02d",cr->GetCluster(),cr->GetCrystal()),5000,0,1e6,cr->GetEnergy());
    Fill(Form("h_segsum_vs_en_clus%02d_crys%02d",cr->GetCluster(),cr->GetCrystal()),1000,0,1e6,cr->GetEnergy(),1000,0,1e6,cr->GetSegmentSum());
    Fill(Form("h_segmult_clus%02d_crys%02d",cr->GetCluster(),cr->GetCrystal()),segs,0,segs,cr->GetMult());
    
    for(int j=0; j<cr->GetMult(); j++){
      Fill(Form("h_segen_vs_nr_clus%02d_crys%02d",cr->GetCluster(),cr->GetCrystal()),segs,0,segs,cr->GetSegmentNr(j),5000,0,1e6,cr->GetSegmentEn(j));
    }
  }
}

#ifdef SIMULATION
void RawHistograms::FillMode2Histograms(Gretina* gr){
  //Mode 2 histograms

  Fill("hmult",30,0,30,gr->GetMult());

  for(int i=0; i<gr->GetMult(); i++){
    Crystal* cr = gr->GetHit(i);
    Fill("hgamma",
	 10000,0,10000,cr->GetEnergy());
    Fill("hsegsum",
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    Fill("hgamma_IPsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetIPSum());
    Fill("hgamma_IPsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetIPSum());
    if (cr->GetEnergy()>300 && cr->GetEnergy()<10000){
      Fill("hgamma_segsum_diff_gated",
	   20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    }
    Fill("hm2_segmult",
	 41,-0.5,40.5,cr->GetMult());
    Fill("herror",
	 10,-0.5,9.5,cr->GetError());
    Fill("herror_ID",
	 150,-0.5,149.5,cr->GetID(),
	 10,-0.5,9.5,cr->GetError());
    Fill("hgamma_error",
	 10,-0.5,9.5,cr->GetError(),
	 500,0,10000,cr->GetEnergy());
    Fill(Form("hgamma_d%d_c%d",fSett->Clu2Det(cr->GetCluster()),cr->GetCrystal()),
	 3500,0,3500,cr->GetEnergy());
    Fill("hhitpattern",
	 150,-0.5,149.5,cr->GetID());
    Fill("hhitpattern_en",
	 150,-0.5,149.5,cr->GetID(),
	 2000,0,2000,cr->GetEnergy());

    for(int i=0; i<4; i++){
      Fill(Form("hcore_e_%d",i),
	   10000,0,3e4,cr->GetCoreE(i));
    }
    Fill("hprestep",
	 1000,0,1000,cr->GetPreStep());
    Fill("hpoststep",
	 1000,0,1000,cr->GetPostStep());

    for(int j=0; j<cr->GetMult(); j++){
      IPoint* IP = cr->GetIPoint(j);
      int segID = cr->GetID()*36 + IP->GetSeg();
      Fill("hseghitpattern",
	   30*4*36,-0.5,30*4*36-0.5,segID);
      Fill("hseghitpattern_en",
	   30*4*36,-0.5,30*4*36-0.5,segID,
	   2000,0,2000,IP->GetSegEnergy());
    }
  }
}
#endif
