#include "UnpackedEvent.hh"

using namespace std;
using namespace TMath;
double rad2deg = 180./TMath::Pi();
double deg2rad = TMath::Pi()/180.;

void swapbytes(char* a, char *b){
  char tmp=*a;
  *a=*b;
  *b=tmp;
}
// Mode 3 data is high endian format
void  HEtoLE(char* cBuf, int bytes){
  for(int i=0; i<bytes; i+=2)
    swapbytes((cBuf+i), (cBuf+i+1));
}


void bits_uint(unsigned int value){
  unsigned int bit;
  for(bit = (~0U >> 1) + 1; bit > 0; bit >>= 1){
    putchar(value & bit ? '1' : '0');
  }
  putchar('\n');
}

UnpackedEvent::UnpackedEvent(Settings* settings){
  if(settings!=NULL){
    fvl = settings->VLevel();
  }
  else{
    fvl = 0;
  }
  fSett = settings;
}

void UnpackedEvent::Init(){
  fmakemode2 = false;

  frhist = new RawHistograms(fSett);
  fchist = new CalHistograms(fSett);

  fEventTimeDiff = fSett->EventTimeDiff();
  fRand = new TRandom();
  fGammaSim = new GammaSim;

  fMode3Event = new Mode3Event;
  fGretina = new Gretina;
  fGretinaCalc = new GretinaCalc;
  fMiniball = new Miniball;
  fMiniballCalc = new MiniballCalc;
  fZeroDeg = new ZeroDeg;
#ifdef USEMINOS
  fMINOS = new MINOS;
#endif
  if(fwtree){
    //setting up tree
    cout << "setting up raw tree " << endl;
    ftr = new TTree("build","built events");
    ftr->Branch("mode3Event",&fMode3Event, 320000);
    ftr->Branch("gretina",&fGretina, 320000);
    ftr->Branch("miniball",&fMiniball, 320000);
    ftr->BranchRef();
    if(fvl>1)
      cout << "done setting up raw tree" << endl;
  }
  if(fwcaltree){
    //setting up tree
    cout << "setting up calibrated tree " << endl;
    fcaltr = new TTree("caltr","calibrated and built events");
    fcaltr->Branch("zerodeg",&fZeroDeg, 320000);
#ifdef USEMINOS
    fcaltr->Branch("minos",&fMINOS, 320000);
#endif
    fcaltr->Branch("gretinacalc",&fGretinaCalc, 320000);
    fcaltr->Branch("miniballcalc",&fMiniballCalc, 320000);
    fcaltr->BranchRef();
    if(fvl>1)
      cout << "done setting up calibrated tree" << endl;
  }
  if(fwsimtree){
    cout << "setting up simulation tree " << endl;
    fsimtr = new TTree("simtr","Geant4 emitted gamma rays");
    fsimtr->Branch("GammaSim",&fGammaSim, 320000);
    fsimtr->BranchRef();
    if(fvl>1)
      cout << "done setting up simulation tree" << endl;
  }
  
  fnentries = 0;
  fGRhits = 0;
  fMBhits = 0;
  fstrangehits = 0;
  fMode3Event->Clear();
  fGretina->Clear();
  fMiniball->Clear();
  fZeroDeg->Clear();
#ifdef USEMINOS
  fMINOS->Clear();
#endif
  fhasdata = false;
  fcurrent_ts = 0;
  ffirst_ts = 0;
  fctr = 0;

  fncalentries = 0;
  fGretinaCalc->Clear();
  fMiniballCalc->Clear();

  ReadSimResolution(fSett->SimResolutionFile());
  ReadSimThresholds(fSett->SimThresholdFile());
  fGammaSim->Clear();

}

int UnpackedEvent::DecodeMode3(char* cBuf, int len, long long int gts){
  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
    fcurrent_ts = gts;
  }
  fctr++;

  unsigned short *wBuf = (unsigned short*)cBuf;
  HEtoLE(cBuf, len);
  // length now in units of 16bit words
  len/=2;
  if(fvl>0)
    cout << "UnpackedEvent: " << " length " << len <<endl;

  Trace curTrace;
  curTrace.Clear();

  int tracesFound = 0;

  //As long as we still have data in the buffer.
  while(len>0){
    tracesFound++;
    if(fvl>2){
      cout << "UnpackedEvent: " << *wBuf << " " << (hex) << *wBuf << (dec) << endl;
      cout << "UnpackedEvent: " << *(wBuf+1) << " " << (hex) << *(wBuf+1) << (dec) << endl;
    }

    // 1st & 2nd word are 0xaaaa
    if( (*wBuf != 0xaaaa) && (*(wBuf+1) != 0xaaaa) ) {
      cerr << "0xAAAA header missing" << endl;
      return 9;
    } else if(fvl>1){
      cout << "UnpackedEvent: " << "-----------------------------next 0xaaaa aaaa" << endl;
    }

    wBuf+=2;

    int length  = (*wBuf & 0x07ff) * 2 + 2;
    len -= length;

    //Read out from buffer and build current trace.
    curTrace = DecodeTrace(&wBuf,length,gts);

    //Trace is now built from input.
    //Time to add it into an event.
 
    static bool errorShown = false;
    if(!errorShown && (curTrace.GetLED() < fcurrent_ts)){
      cout << "UnpackedEvent: " << "Found events that are not in order" << endl;
      cout << "UnpackedEvent: " << "Have you run this through GEB_HFC first?" << endl;
      errorShown = true;
    }

    //Enough global time has passed, and so I am now in a new event.
    //Close and write the event, then clear it out.
    if(curTrace.GetLED() - fcurrent_ts > fEventTimeDiff || curTrace.GetLED() < fcurrent_ts){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in Mode3 Entry." << endl;
      if(fvl>1&&fMode3Event->GetMult()<1)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__ << " entry " << fnentries << " mode3 is empty " << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__ << " entry " << fnentries << " s800 is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    //The event is now written to the tree and is cleared if we decided to make a new event.
    //The event, whether newly made or made previously, gets passed a new trace.
    fMode3Event->AddTrace(curTrace);

    if(fvl>1)
      cout << "UnpackedEvent: " << "setting current ts " << curTrace.GetLED() << endl;

    fcurrent_ts = curTrace.GetLED();

    if(fvl>2)
      cout << "UnpackedEvent: " << "remaining " << len << endl;

  }
  if(fvl>2){
    cout << "UnpackedEvent: " << "Mode3 event found with timestamp " << gts << endl;
    cout << "UnpackedEvent: " << tracesFound << " traces decoded from the mode3 entry" << endl;
  }
  return 0;
}

Trace UnpackedEvent::DecodeTrace(unsigned short** wBuf_p, int length, long long int gts){
  unsigned short* wBuf = *wBuf_p;

  Trace curTrace;
  curTrace.Clear();


  curTrace.SetTS(gts);
  curTrace.SetLength(length-16);
  // this length includes aaaa aaaa until next aaaa aaaa.
  // trace length is this minus 16

  if(fvl>2)
    cout << "UnpackedEvent: " << "trace length " << curTrace.GetLength() << endl;

  curTrace.SetBoard(*wBuf >> 11);//called GA in mario's doc 13.dez
  wBuf++;
  if(fvl>4){
    cout << "UnpackedEvent: " << "---------------------------------------------------------------------------------"<<endl;
    cout << "UnpackedEvent: " << *wBuf << "\t" << (hex) << *wBuf << (dec) <<"\t";
    bits_uint(*wBuf);
    cout << "UnpackedEvent: " << (*wBuf & 0x000f) << "\t" << (hex) << (*wBuf & 0x000f) << (dec) <<"\t";
    bits_uint(*wBuf & 0x000f);
    cout << "UnpackedEvent: " << ((*wBuf & 0x0030)>>4) << "\t" << (hex) << ((*wBuf & 0x0030)>>4) << (dec) <<"\t";
    bits_uint((*wBuf & 0x0030)>>4);
    cout << "UnpackedEvent: " << ((*wBuf & 0x00c0)>>6) << "\t" << (hex) << ((*wBuf & 0x00c0)>>6) << (dec) <<"\t";
    bits_uint((*wBuf & 0x00c0)>>6);
    cout << "UnpackedEvent: " << ((*wBuf & 0x1f00)>>8) << "\t" << (hex) << ((*wBuf & 0x1f00)>>8) << (dec) <<"\t";
    bits_uint((*wBuf & 0x1f00)>>8);
  }

  curTrace.SetChn(*wBuf & 0x000f);
  curTrace.SetSlot((*wBuf & 0x0030)>>4);
  curTrace.SetCrystal((*wBuf & 0x00c0)>>6);
  curTrace.SetHole((*wBuf & 0x1f00)>>8);
  wBuf++;


  if(fvl>2)
    cout << "UnpackedEvent: " << "board " << curTrace.GetBoard() << " chn " << curTrace.GetChn() << " slot " << curTrace.GetSlot() << " cry " << curTrace.GetCrystal() << " det " << curTrace.GetHole() << endl;

  int id = curTrace.GetHole()*4 + curTrace.GetCrystal();

  // point 5th
  long long int ts;
  ts = (long long int) *(wBuf+1);
  ts += ((long long int) *(wBuf+0)) << 16;
  ts += ((long long int) *(wBuf+3)) << 32;
  curTrace.SetLED(ts);

  if(fvl>2)
    cout << "UnpackedEvent: " << "led " << curTrace.GetLED() << endl;

  wBuf+=2; //point 7th

  int en = (int) *(wBuf+3) & 0x00ff;
  en = en << 16;
  en += (int) *(wBuf);
  int sign = *(wBuf+3) & 0x0100;

  if(sign){
    if(curTrace.GetChn()==9) //core
      en = (int)(en - (int)0x01000000); // (2^24)
    else{
      en = (int)(en - (int)0x01000000); // (2^24)
      en = - en;
    }
  }
  else{
    if(curTrace.GetChn()!=9) // not core
      en = - en;
  }


  curTrace.SetEnergy(en);
  curTrace.SetEnSign(sign);
  curTrace.SetPileUp((*(wBuf+3) & 0x8000) >> 15);
  wBuf+=2; //point 9th

  if(fvl>2)
    cout << "UnpackedEvent: " << "energy " << curTrace.GetEnergy() << " sign " << curTrace.GetEnSign() << " pileup " << curTrace.GetPileUp() << endl;

  ts = (long long int) *(wBuf+0);
  ts += ((long long int) *(wBuf+3)) << 16;
  ts += ((long long int) *(wBuf+2)) << 32;
  curTrace.SetCFD(ts);

  if(fvl>2)
    cout << "UnpackedEvent: " << "cfd " << curTrace.GetCFD() << endl;

  wBuf+=4; //point 13th

  int cfd = (int) *(wBuf+1) << 16;
  cfd += (int) *wBuf;
  curTrace.SetCFD(0,cfd);
  wBuf+=2; //point 15th

  cfd = (int) *(wBuf+1) << 16;
  cfd += (int) *wBuf;
  curTrace.SetCFD(1,cfd);
  wBuf+=2; //point 17th

  if(fvl>2)
    cout << "UnpackedEvent: " << "cfd points " << curTrace.GetCFD(0) << " and " << curTrace.GetCFD(1) << endl;

  //cout << "UnpackedEvent: " << curTrace.GetTrace().size() << " size" << endl;
  //  int bg =0;
  //  int sig =0;
  for(int i=0; i<(length-16); i+=2){
    //cout << "UnpackedEvent: " << i <<"\t"<< *(wBuf) << endl;
    curTrace.SetTrace(i ,-(short) *(wBuf+1)+512 );
    curTrace.SetTrace(i+1 ,-(short) *(wBuf)+512 );
    // if(id==124){
    //   curTrace.SetLaBrTrace(i ,-(short) *(wBuf+1)+512 );
    //   curTrace.SetLaBrTrace(i+1 ,-(short) *(wBuf)+512 );
    //   if(i>19&&i<71){

    // 	bg+=-(short) *(wBuf)+512;
    // 	bg+=-(short) *(wBuf+1)+512;
    // 	if(fvl>2)
    // 	  cout << "UnpackedEvent: " << "labr bg " << bg << endl;
    //   }
    //   if(i>69&&i<121){

    // 	sig+=-(short) *(wBuf)+512;
    // 	sig+=-(short) *(wBuf+1)+512;
    // 	if(fvl>2)
    // 	  cout << "UnpackedEvent: " << "labr sig " << sig << endl;
    //   }
    // }
    wBuf++;
    wBuf++;


    if(fvl>5){
      cout << "UnpackedEvent: " << i <<"\t"<< curTrace.GetTrace()[i] << endl;
      cout << "UnpackedEvent: " << i+1 <<"\t"<< curTrace.GetTrace()[i+1] << endl;
    }
  }
  // if(id==124){
  //   if(fvl>1){
  //     cout << "UnpackedEvent: " << "end fof trace bg " << bg << " sig " << sig << endl;
  //     cout << "UnpackedEvent: " << " energy " << sig-bg << endl;
  //   }
  //   curTrace.SetLaBr(sig-bg);
  // }
  if(fvl>2){
    cout << " these next ones should be aaaa and aaaa again " << endl;
    cout << "UnpackedEvent: " << *wBuf << " " << (hex) << *wBuf << (dec) << endl;
    cout << "UnpackedEvent: " << *(wBuf+1) << " " << (hex) << *(wBuf+1) << (dec) << endl;
  }


  //needs to be modified for MB and SC !!
  if(curTrace.GetBoard() == 6 && curTrace.GetChn() == 9){ //board 3 = 10 , 6 = 5 MeV range
    if(fvl>1){
      cout << "UnpackedEvent: " << fctr << " core hit id: "<<id<<"\tts: " << curTrace.GetLED() << endl;
    }
    curTrace.SetCore(true);
  }
  else{
    if(fvl>1){
      cout << "UnpackedEvent: " << fctr << " segm hit id: "<<id<<"\tts: " << curTrace.GetLED() << endl;
    }
  }
  //end needs modification

  *wBuf_p = wBuf;

  return curTrace;
}

int UnpackedEvent::DecodeGretina(Crystal* cryst, long long int gts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << gts << endl;
  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
  }
  if(fvl>1){
    cout << "UnpackedEvent: " << "-----------------------------next hit: "<< fnentries<< endl;
    cryst->PrintEvent();
  }
  int det = cryst->GetCluster();
  int cry = cryst->GetCrystal();

  // if(det<0 || det > MAXDETPOS){ // i don't know whether det runs from 0 or 1
  //   cout << "UnpackedEvent: " << "invalid detector number " << det << endl;
  //   return 11;
  // }
  if(cry<0 || cry > MAXCRYSTALNO-1){
    cout << "UnpackedEvent: " << "invalid crystal number " << cry << endl;
    return 12;
  }



  if(fvl>1){
    cout << "UnpackedEvent: " << "mult " << cryst->GetMult() <<"\ten " <<  cryst->GetEnergy() <<"\tts " <<  cryst->GetTS() <<"\ten max " <<  cryst->GetMaxEn()<<"\tip max " <<  cryst->GetMaxIPNr()<< endl;
    cout << "UnpackedEvent: " <<"-----------------------------"<< endl;
  }

  fGRhits++;

  //the events which have no good interaction points
  if(cryst->GetMaxIPNr()<0){
    fstrangehits++;
  }
  if(fvl>0 && cryst->GetMaxIPNr()<0){
    cryst->PrintEvent();;
  }
  // now check time stamps
  long long int deltaEvent = gts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (gretina) " << gts << " difference " << deltaEvent<< endl;

  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< "this ts " << gts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
  }
  else{
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " gretina single time difference " << deltaEvent << endl;
    if(fcurrent_ts>-1){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in Gretina." << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }




  fGretina->AddHit(cryst);
  //set the current timestamp
  fcurrent_ts = gts;

  if(fvl>2){
    cout << "UnpackedEvent: " << "Gretina event found with timestamp " << gts << endl;
  }
  return 0;
}

int UnpackedEvent::DecodeMiniball(MBCrystal* cryst, long long int gts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << gts << endl;
  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
  }
  if(fvl>1){
    cout << "UnpackedEvent: " << "-----------------------------next hit: "<< fnentries<< endl;
  }
  int det = cryst->GetCluster();
  int cry = cryst->GetCrystal();

  if(cry<0 || cry > MBCRYST){
    cout << "UnpackedEvent: " << "invalid MB crystal number " << cry << endl;
    return 12;
  }



  if(fvl>1){
    cout << "UnpackedEvent: mult " << cryst->GetMult() <<"\ten " <<  cryst->GetEnergy() <<"\tts " <<  cryst->GetTS() <<"\tmax seg en " <<  cryst->GetMaxSegEn()<<"\tmax seg nr" <<  cryst->GetMaxSegNr()<< endl;
    cout << "UnpackedEvent: " <<"-----------------------------"<< endl;
  }

  fMBhits++;

  //the events which have no good interaction points
  if(cryst->GetMaxSegNr()<0){
    fstrangehits++;
  }
  if(fvl>0 && cryst->GetMaxSegNr()<0){
    cout << "strange hit " << endl;
    cryst->PrintEvent();
  }
  // now check time stamps
  long long int deltaEvent = gts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (gretina) " << gts << " difference " << deltaEvent<< endl;

  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< "this ts " << gts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
  }
  else{
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " miniball single time difference " << deltaEvent << endl;
    if(fcurrent_ts>-1){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in Gretina." << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }

  fMiniball->AddHit(cryst);
  if(fvl>1)
    fMiniball->PrintEvent();
  //set the current timestamp
  fcurrent_ts = gts;

  if(fvl>2){
    cout << "UnpackedEvent: " << "Miniball event found with timestamp " << gts << endl;
  }
  return 0;
}

int UnpackedEvent::DecodeGammaG4Sim(G4SIM_EGS* g4Sim, long long int ts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  fGammaSim->SetTimeStamp(ts);
  for(Int_t i = 0; i < g4Sim->num; i++){
    fGammaSim->AddEmittedGamma(&g4Sim->gammas[i]);
  }

  if(fvl>2){
    cout << "Simulation: GammaSim: " << endl;
    for(Int_t i = 0; i < fGammaSim->GetMult(); i++)
      cout << "e = " << fGammaSim->GetEmittedGamma(i)->GetEnergy()
	   << " (x, y, z) = (" 
	   << fGammaSim->GetEmittedGamma(i)->GetPos().X() << ", "
	   << fGammaSim->GetEmittedGamma(i)->GetPos().Y() << ", "
	   << fGammaSim->GetEmittedGamma(i)->GetPos().Z() 
	   << ")  (phi, theta) = ("
	   << fGammaSim->GetEmittedGamma(i)->GetPhi() << ", "
	   << fGammaSim->GetEmittedGamma(i)->GetTheta() << ")" << endl;

    cout << "Simulation: GammaSim event found with timestamp " << ts << endl;
  }

  // These events go into a separate tree as singles. No need to pay 
  // attention to the time window.
  
  fsimtr->Fill();
  fnsimentries++;
  fGammaSim->Clear();

  return 0;
}

int UnpackedEvent::DecodeZeroDegPhysicsData(ZD_PHYSICSDATA* zeropd, long long int ts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  // now check time stamps
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (zerodeg) " << ts << " difference " << deltaEvent<< endl;
  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
    if(fhasdata==true){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to two ZeroDeg entries." << endl;
      cout << "UnpackedEvent: " << " coincidence with another ZeroDeg! deltaT = " << deltaEvent << " writing and clearing last event" << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " zerodeg is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
  } else {
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " zerodeg single event difference " << deltaEvent << endl;
    if(fcurrent_ts>-1&&fhasdata==true){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in ZeroDeg." << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " zerodeg is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }
  //set the current timestamp
  fcurrent_ts = ts;

  if(fvl>0)
    cout << "UnpackedEvent: ZD_PHYSICSDATA: ata = " << zeropd->ata 
	 << " bta = "  << zeropd->bta 
	 << " xta = "  << zeropd->xta 
	 << " yta = "  << zeropd->yta 
	 << " betata = "  << zeropd->betata
	 << " ts = "  << ts << endl;

  fZeroDeg->SetTimeStamp(ts);
  fZeroDeg->SetATA(fRand->Gaus(zeropd->ata, fSett->TargetAngleResolution()));
  fZeroDeg->SetBTA(fRand->Gaus(zeropd->bta, fSett->TargetAngleResolution()));
  fZeroDeg->SetXTA(fRand->Gaus(zeropd->xta, fSett->TargetPosResolution()));
  fZeroDeg->SetYTA(fRand->Gaus(zeropd->yta, fSett->TargetPosResolution()));
  fZeroDeg->SetBetaTA(fRand->Gaus(zeropd->betata, fSett->TargetBetaResolution()));
  fhasdata = true;
  if(fvl>2){
    cout << "UnpackedEvent: ZeroDegPhysicsData event found with timestamp " << ts << endl;
  }
  return 0;
}

#ifdef USEMINOS
int UnpackedEvent::DecodeMINOSPhysicsData(MINOS_DATA* minospd, long long int ts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  // now check time stamps
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (minos) " << ts << " difference " << deltaEvent<< endl;
  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }
  
  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
  } else {
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " zerodeg single event difference " << deltaEvent << endl;
    if(fcurrent_ts>-1){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in MINOS." << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }
  //set the current timestamp
  fcurrent_ts = ts;

  if(fvl>0)
    cout << "UnpackedEvent: MINOS_DATA: x = " << minospd->x 
	 << " y = "  << minospd->y 
	 << " z = "  << minospd->z 
	 << " betare = "  << minospd->betare
	 << " ts = "  << ts << endl;

  fMINOS->SetTimeStamp(ts);
  double x = fRand->Gaus(minospd->x, fSett->MINOSXYResolution()); 
  double y = fRand->Gaus(minospd->y, fSett->MINOSXYResolution()); 
  double z = fRand->Gaus(minospd->z, fSett->MINOSZResolution()); 
  fMINOS->SetVertex(x,y,z);
  fMINOS->SetBetaRE(minospd->betare);
  if(fvl>2){
    cout << "UnpackedEvent: MINOSPhysicsData event found with timestamp " << ts << endl;
  }
  return 0;
}
#endif

void UnpackedEvent::ReadSimResolution(const char* filename){
  TEnv* env = new TEnv(filename);
   if(fvl>1){
     env->Print();
   }
  double globA = env->GetValue("Detector.All.Crystal.All.A",0.0);
  double globB = env->GetValue("Detector.All.Crystal.All.B",0.0);
  double globC = env->GetValue("Detector.All.Crystal.All.C",0.0);
  for(int det=0;det<MAXDETPOS;det++){
    for(int cr=0;cr<4;cr++){
      int toRes = env->GetValue(Form("Detector.%d.Crystal.%d.SimResolution",det,cr),0);
      if(toRes){
	simresolution newSimres;
	newSimres.Detector = det;
	newSimres.Crystal = cr;
	newSimres.A = env->GetValue(Form("Detector.%d.Crystal.%d.A",det,cr),0.0);
	newSimres.B = env->GetValue(Form("Detector.%d.Crystal.%d.B",det,cr),0.0);
	newSimres.C = env->GetValue(Form("Detector.%d.Crystal.%d.C",det,cr),0.0);
	fSimResolutions.push_back(newSimres);
      } else {
	simresolution newSimres;
	newSimres.Detector = det;
	newSimres.Crystal = cr;
	newSimres.A = globA;
	newSimres.B = globB;
	newSimres.C = globC;
	fSimResolutions.push_back(newSimres);
      }
    }
  }
}

void UnpackedEvent::ReadSimThresholds(const char* filename){
  TEnv* env = new TEnv(filename);
  double globE = env->GetValue("Detector.All.Crystal.All.E",0.0);
  double globdE = env->GetValue("Detector.All.Crystal.All.dE",0.001);
  cout << "reading file " << filename << endl;
  
  for(int det=0;det<MAXDETPOS;det++){
    for(int cr=0;cr<4;cr++){
      int toThresh = env->GetValue(Form("Detector.%d.Crystal.%d.SimThreshold",det,cr),0);
      if(toThresh){
	simthreshold newSimthresh;
	newSimthresh.Detector = det;
	newSimthresh.Crystal = cr;
	newSimthresh.E = env->GetValue(Form("Detector.%d.Crystal.%d.E",det,cr),0.0);
	newSimthresh.dE = env->GetValue(Form("Detector.%d.Crystal.%d.dE",det,cr),0.001);
	fSimThresholds.push_back(newSimthresh);
      }
      else{
	simthreshold newSimthresh;
	newSimthresh.Detector = det;
	newSimthresh.Crystal = cr;
	newSimthresh.E = globE;
	newSimthresh.dE = globdE;
	fSimThresholds.push_back(newSimthresh);
      }
    }
  }
}
//! apply the resolution
bool UnpackedEvent::SimResolution(Gretina* gr){
  //cout << __PRETTY_FUNCTION__ << " with resolution " << fSett->SimGretinaPositionResolution() << endl;
  for(int i=0; i<gr->GetMult(); i++){
    Crystal* crys = gr->GetHit(i);
    //Position resolution
    if(fSett->SimGretinaPositionResolution()>0){
      for(int j=0; j<crys->GetMult(); j++){
	crys->GetIPoint(j)->SetPosition(fRand->Gaus(crys->GetIPoint(j)->GetPosition().X(),
						    fSett->SimGretinaPositionResolution()),
					fRand->Gaus(crys->GetIPoint(j)->GetPosition().Y(),
						    fSett->SimGretinaPositionResolution()),
					fRand->Gaus(crys->GetIPoint(j)->GetPosition().Z(),
						    fSett->SimGretinaPositionResolution()) );
      }
    }
    //Energy resolutions
    int det = fSett->Clu2Det(crys->GetCluster());
    int crysnum = crys->GetCrystal();
    //cout << crys->GetCluster() << "\t" << det << "\t" << crysnum << endl;
    //cout << crys->GetEnergy() << " ----> ";
    for(vector<simresolution>::iterator it = fSimResolutions.begin(); it!=fSimResolutions.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	crys->SetEnergy(fRand->Gaus(crys->GetEnergy(),
				    it->A*sqrt(1.0+crys->GetEnergy()*it->B) + it->C*crys->GetEnergy()));
	for(int j=0; j<crys->GetMult(); j++){
	  IPoint* ipoint = crys->GetIPoint(j);
	  ipoint->SetEnergy(fRand->Gaus(ipoint->GetEnergy(),
					it->A*sqrt(1.0+ipoint->GetEnergy()*it->B) + it->C*ipoint->GetEnergy()));
	}
	break;
      }
    }
    //cout << crys->GetEnergy() << endl;
  }

  return true;
}

//! apply the resolution
bool UnpackedEvent::SimResolution(Miniball* mb){
  //cout << __PRETTY_FUNCTION__ << " with resolution " << fSett->SimGretinaPositionResolution() << endl;
  for(int i=0; i<mb->GetMult(); i++){
    MBCrystal* crys = mb->GetHit(i);
    //Energy resolutions
    int det = fSett->Clu2Det(crys->GetCluster());
    int crysnum = crys->GetCrystal();
    //cout << crys->GetCluster() << "\t" << det << "\t" << crysnum << endl;
    for(vector<simresolution>::iterator it = fSimResolutions.begin(); it!=fSimResolutions.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	//cout << it->A << "\t" << it->B << "\t" << it->C << "\t" <<crys->GetEnergy();
	crys->SetEnergy(fRand->Gaus(crys->GetEnergy(),
				    it->A*sqrt(1.0+crys->GetEnergy()*it->B) + it->C*crys->GetEnergy()));
	//segments
	for(int j=0; j<crys->GetMult(); j++){
	  crys->SetSegmentEn(j,fRand->Gaus(crys->GetSegmentEn(j),
					it->A*sqrt(1.0+crys->GetSegmentEn(j)*it->B) + it->C*crys->GetSegmentEn(j)));
	}
	//cout << "\t" <<crys->GetEnergy()<< endl;
	break;
      }
    }
  }

  return true;
}

bool UnpackedEvent::SimThresholds(Gretina* gr){
  //cout << __PRETTY_FUNCTION__ << endl;
  for(int i=0; i<gr->GetMult(); i++){
    Crystal* crys = gr->GetHit(i);
    //cout << crys->GetEnergy() << " -> ";
    int det = fSett->Clu2Det(crys->GetCluster());
    int crysnum = crys->GetCrystal();
    //cout << "det " << det << ", clu " << crys->GetCluster() << ", crys " << crysnum << endl;
    for(vector<simthreshold>::iterator it = fSimThresholds.begin(); it!=fSimThresholds.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	if( fRand->Uniform(0,1) >
	    0.5*(1.0 + tanh( (crys->GetEnergy() - it->E) / it->dE ) ) ){
	  crys->SetEnergy(0.);
	}
	break;
      }
    }
    //cout << crys->GetEnergy() <<endl;
  }

  return true;
}
bool UnpackedEvent::SimThresholds(Miniball* mb){
  //cout << __PRETTY_FUNCTION__ << endl;
  for(int i=0; i<mb->GetMult(); i++){
    MBCrystal* crys = mb->GetHit(i);
    int det = fSett->Clu2Det(crys->GetCluster());
    int crysnum = crys->GetCrystal();
    for(vector<simthreshold>::iterator it = fSimThresholds.begin(); it!=fSimThresholds.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	if( fRand->Uniform(0,1) >
	    0.5*(1.0 + tanh( (crys->GetEnergy() - it->E) / it->dE ) ) ){
	  crys->SetEnergy(0.);
	}
	break;
      }
    }
  }

  return true;
}
void UnpackedEvent::ClearEvent(){
  fhasdata = false;
  fMode3Event->Clear();
  fGretina->Clear();
  fGretinaCalc->Clear();
  fMiniball->Clear();
  fMiniballCalc->Clear();
  fZeroDeg->Clear();
#ifdef USEMINOS
  fMINOS->Clear();
#endif
  return;
}

/*!
  Performs end of event actions.
  The actions performed depend on the flags given.
  - rt, write out the raw tree.
  - rh, fill raw histograms, to be written later.
  - ct, write out the calibrated tree.
  - ch, fill calibrated histograms, to be written later.

  Note that calibrations are only performed in GrROOT if either "ct" or "ch" are given as flags.
*/
void UnpackedEvent::CloseEvent(){
  //cout << __PRETTY_FUNCTION__ << endl;
  SimResolution(fGretina);
  SimThresholds(fGretina);
  SimResolution(fMiniball);
  SimThresholds(fMiniball);
  if(fwtree || fwhist){
    if(fmakemode2){
      //cout << "fMode3Event->GetMult() " << fMode3Event->GetMult() <<"\tfMiniball->GetMult() " << fMiniball->GetMult() << endl;   
      MakeMode2();
      //cout << "fMode3Event->GetMult() " << fMode3Event->GetMult() <<"\tfMiniball->GetMult() " << fMiniball->GetMult() << "-----------after " << endl;   
    }


    if (fwhist){
#ifdef USEMINOS
      frhist->FillHistograms(fGretina,fMiniball,fZeroDeg,fMINOS);
#else 
      frhist->FillHistograms(fGretina,fMiniball,fZeroDeg,NULL);
#endif
    }
    //Write the raw tree.
    if (fwtree){
      ftr->Fill();
    }
    fnentries++;
  }
  if(fwcaltree||fwcalhist){
    

    //Build all of the calibrated objects, using the calibration in cal.
    //Use the data from the first three parameters, output into the last three parameters.
    //cout << "calculation BuildAllCalc called" << endl;
    //    if(trackMe)
    //      fcal->GammaTrack(fgretinaCalc,fgretinaEvent);

#ifdef USEMINOS
    fcal->BuildAllCalc(fGretina,fGretinaCalc,fMiniball, fMiniballCalc,fZeroDeg,fMINOS);
#else
    fcal->BuildAllCalc(fGretina,fGretinaCalc,fMiniball, fMiniballCalc,fZeroDeg,NULL);
#endif
    if(fwcaltree){
      fcaltr->Fill();
      fncalentries++;
    }
    if(fwcalhist){
#ifdef USEMINOS
      fchist->FillHistograms(fGretinaCalc, fMiniballCalc,fZeroDeg,fMINOS);
#else
      fchist->FillHistograms(fGretinaCalc, fMiniballCalc,fZeroDeg,NULL);
#endif
    }
  }
  
  this->ClearEvent();
}
/*!
  Write the last event to file.
  This event would not be written otherwise, as there are no following buffers.
  Usually, events are written whenever the following entry is past some time window.
*/
void UnpackedEvent::WriteLastEvent(){
  if(fvl>2)
    cout << "UnpackedEvent: " << "last event " << endl;
  fMode3Event->SetCounter(fctr);
  fctr = 0;
  this->CloseEvent();

  if(fwhist){
    frhist->Write();
  }
  if(fwcalhist){
    fchist->Write();
  }
  if(fvl>1)
    fcal->PrintCtrs();
}
/*! 
  Create Miniball object from mode3 data
*/
void UnpackedEvent::MakeMode2(){
  fMiniball->Clear();
  fGretina->Clear();
  //cout << __PRETTY_FUNCTION__ << endl;
  //cout << " mult " << fMode3Event->GetMult() << endl;
  for(int i=0; i<fMode3Event->GetMult(); i++){
    Mode3Hit* hit = fMode3Event->GetHit(i);
    //cout << "hit " << i << " mult " << hit->GetMult() << endl;
    for(int j=0; j<hit->GetMult(); j++){
      Trace * trace = hit->GetTrace(j);
      if(fSett->VLevel()>1){
	cout << "Trace " << j << " Length " << trace->GetLength() << 
	  "\tEnergy " << trace->GetEnergy() <<
	  "\tBoard " << trace->GetBoard() <<
	  "\tSlot " << trace->GetSlot() <<
	  "\tChannel " << trace->GetChn() <<
	  "\tHole " << trace->GetHole() <<
	  "\tCrystal " << trace->GetCrystal() << 
	  "\tTimeStamp " << trace->GetTS();
	if(trace->GetChn()==9)
	  cout << " <- CC " << endl;
	else if(trace->GetEnergy()>1000)
	  cout << " <- with net energy " << endl;
	else
	  cout << endl;
	cout << fSett->MiniballModule(trace->GetHole(),trace->GetCrystal(),trace->GetSlot()) << "\t" << fSett->MiniballCrystal(trace->GetHole(),trace->GetCrystal(),trace->GetSlot()) << endl;
	
      }
      bool tracking = false;
      int mod = fSett->MiniballModule(trace->GetHole(),trace->GetCrystal(),trace->GetSlot());
      int cry = fSett->MiniballCrystal(trace->GetHole(),trace->GetCrystal(),trace->GetSlot());
      // for tracking detectors, channel will be 0-39 = chn+slot*10
      if(mod>9){
	mod-=10;
	tracking =true;
	continue;
      }

      int chn = trace->GetChn();
      if(mod<0 || cry<0){
	cout << "invalid MB module or crystal " << endl;
	continue;
      }
      MBCrystal* mbhit = fMiniball->GetHit(mod,cry);
      if(mbhit){
	if(chn==9)
	  mbhit->InsertCore(mod, cry, trace->GetEnergy(), trace->GetTS());
	else
	  mbhit->InsertSegment(mod, cry, chn, trace->GetEnergy());	  
      }      
      else{
	fMiniball->AddHit(new MBCrystal(mod, cry, chn, trace->GetEnergy(), trace->GetTS()));
      }      
    }//traces
  }//hits
  if(fSett->VLevel()>1)
    fMiniball->PrintEvent();
}
