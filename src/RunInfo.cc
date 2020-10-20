#include "RunInfo.hh"
using namespace std;
void RunInfo::Clear(){
  fHIrunnr = -1;
  fBRrunnr = -1;
  fM2runnr = -1;
  fHIbytes = -1;
  fHIrawevents = -1;
  fHIcalevents = -1;
  fMode2events = -1;
  fBRevents = -1;
}

void RunInfo::AppendInfo(RunInfo* sec){
  if(fHIrunnr<0 && sec->GetHIRunNumber()>0){
    fHIrunnr = sec->GetHIRunNumber();
    fHIbytes = sec->GetHIBytes();
    fHIrawevents = sec->GetHIRawEvents();
  }
  else
    cout <<  "HiCARI info already there, run number set to " << fHIrunnr << endl;
  if(fHIcalevents<0 && sec->GetHICalEvents()>0){
    fHIcalevents = sec->GetHICalEvents();
    fBigRIPSctr = sec->GetBigRIPSCtr(); 
    fBigRIPSHitctr = sec->GetBigRIPSHitCtr();
    fHiCARIctr = sec->GetHiCARICtr();
    fHiCARIHitctr = sec->GetHiCARIHitCtr();
  }
  if(fM2runnr<0 && sec->GetM2RunNumber()>0)
    fM2runnr = sec->GetM2RunNumber();
  else
    cout <<  "Mode2 info already there, run number set to " << fM2runnr << endl;
  if(fMode2events<0 && sec->GetMode2Events()>0){
    fMode2events = sec->GetMode2Events();
    fMode2Hitctr = sec->GetMode2HitCtr();
  }
    

  if(fBRrunnr<0 && sec->GetBRRunNumber()>0)
    fBRrunnr = sec->GetBRRunNumber();
  else
    cout <<  "BigRIPS info already there, run number set to " << fBRrunnr << endl;


}


void RunInfo::PrintRunInfo(){
  if(fHIrunnr>0){
    cout << "HiCARI Run Number as read from inputfile:\t" << fHIrunnr << endl; 
    cout << "Number of Raw Bytes read from HiCARI file:\t" << fHIbytes << endl;
    cout << "Number of Raw Events built for HiCARI:\t" << fHIrawevents << endl;
    cout << "Number of Cal Events built for HiCARI:\t" << fHIcalevents << endl;
    if(fHIcalevents){
      cout << "HiCARI events \t" << fHiCARIctr  << endl;
      cout << "HiCARI hits   \t" << fHiCARIHitctr  << endl;
      cout << "BigRIPS events\t" << fBigRIPSctr  << endl;
      cout << "BigRIPS hits  \t" << fBigRIPSHitctr  << endl;
    }
  }
  if(fM2runnr>0){
    cout << "Mode2 Run Number as read from inputfile:\t" << fM2runnr << endl; 
    cout << "Number of Mode2 Events built for HiCARI:\t" << fMode2events << endl;
    if(fMode2events){
      cout << "Mode2 hits   \t" << fMode2Hitctr  << endl;
    }
  }
  
  if(fBRrunnr>0){
    cout << "BigRIPS Run Number as read from inputfile:\t" << fBRrunnr << endl; 
    cout << "Number of BigRIPS Events read from RIDF file:\t" << fBRevents << endl;
  }

}
