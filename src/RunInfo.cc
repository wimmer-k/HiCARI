#include "RunInfo.hh"
using namespace std;
void RunInfo::Clear(){
  fHIrunnr = -1;
  fBRrunnr = -1;
  fM2runnr = -1;
  fHIbytes = -1;
  fHIrawevents = -1;
  fHIcalevents = -1;
  fM2bytes = -1;
  fM2rawevents = -1;
  fM2calevents = -1;
  fBRevents = -1;
}

void RunInfo::AppendInfo(RunInfo* sec){
  // cout << fHIrunnr << "\t" << fM2runnr << "\t" << fBRrunnr  << endl;
  // cout << " adding " <<  sec->GetHIRunNumber() << "\t" <<  sec->GetM2RunNumber() << "\t" << sec->GetBRRunNumber()<< endl;
  if(fHIrunnr<0 && sec->GetHIRunNumber()>0){
    fHIrunnr = sec->GetHIRunNumber();
    fHIbytes = sec->GetHIBytes();
    fHIrawevents = sec->GetHIRawEvents();
  }
  else if(fHIrunnr>-1 && sec->GetHIRunNumber()>0)
    cout <<  "HiCARI info already there, run number set to " << fHIrunnr << endl;
  if(fHIcalevents<0 && sec->GetHICalEvents()>0){
    fHIcalevents = sec->GetHICalEvents();
    fBigRIPSctr = sec->GetBigRIPSCtr(); 
    fBigRIPSHitctr = sec->GetBigRIPSHitCtr();
    fHiCARIctr = sec->GetHiCARICtr();
    fHiCARIHitctr = sec->GetHiCARIHitCtr();
  }
  if(fM2runnr<0 && sec->GetM2RunNumber()>0){
    fM2runnr = sec->GetM2RunNumber();
    fM2bytes = sec->GetM2Bytes();
  }
  else if (fM2runnr>-1 && sec->GetM2RunNumber()>0)
    cout <<  "Mode2 info already there, run number set to " << fM2runnr << endl;
  if(fM2calevents<0 && sec->GetM2CalEvents()>0){
    fM2calevents = sec->GetM2CalEvents();
    fGretinactr = sec->GetGretinaCtr();
    fGretinaHitctr = sec->GetGretinaHitCtr();
    fGretinaHitABctr = sec->GetGretinaHitABCtr();
  }
    

  if(fBRrunnr<0 && sec->GetBRRunNumber()>0)
    fBRrunnr = sec->GetBRRunNumber();
  else if (fBRrunnr>-1 && sec->GetBRRunNumber()>0)
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
    //cout << "Number of Mode2 Raw Events built for HiCARI:\t" << fM2rawevents << endl;
    cout << "Number of Gretina Events built for HiCARI:\t" << fM2calevents << endl;
    if(fM2calevents){
      cout << "Gretina events   \t" << fGretinactr  << endl;       
      cout << "Gretina hits     \t" << fGretinaHitctr  << endl;
      cout << "Gretina hits AB  \t" << fGretinaHitABctr  << endl;
    }
  }
  
  if(fBRrunnr>0){
    cout << "BigRIPS Run Number as read from inputfile:\t" << fBRrunnr << endl; 
    cout << "Number of BigRIPS Events read from RIDF file:\t" << fBRevents << endl;
  }

}
