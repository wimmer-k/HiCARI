#include "BuildEvents.hh"
using namespace std;

  //! Default constructor
BuildEvents::BuildEvents(Settings* set){
  fSett = set;
  SetWindow(set->EventTimeDiff());
  SetCorrelationMode(set->CorrelationMode()); // for p-gamma coinc choose 1  
};

/*!
  Initialyze the event building
  \param brtr tree with input bigrips data
  \param hitr tree with input hicari data
*/
void BuildEvents::Init(TTree* brtr, TTree* hitr, TTree* m2tr){

  flastevent = -1;
  fhasBR = false;
  fhasHI = false;
  fhasM2 = false;
  fBRentries = 0;
  fHIentries = 0;
  fM2entries = 0;
  
  if(brtr!=NULL)
    fhasBR = true;
  if(hitr!=NULL)
    fhasHI = true;
  if(m2tr!=NULL)
    fhasM2 = true;
  
  fBRtr = brtr;
  fHItr = hitr;
  fM2tr = m2tr;

  fBRentry = 0;
  fHIentry = 0;
  fM2entry = 0;
  fnbytes = 0;

  fBRts = 0;
  fcheckADC = -1;
  ftrigbit = -1;
  fbeam = new Beam;
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      ffp[f] = new FocalPlane;
    }
  }
  fHIts = 0;
  fhicari = new HiCARICalc;
  fM2ts = 0;
  fmode2 = new GretinaCalc;

  //local copy for intermediate storing
  flocalBRts = 0;
  flocalcheckADC = -1;
  flocaltrigbit = -1;
  flocalbeam = new Beam;
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      flocalfp[f] = new FocalPlane;
    }
  }
  flocalHIts = 0;
  flocalhicari = new HiCARICalc;
  flocalM2ts = 0;
  flocalmode2 = new GretinaCalc;

  if(fhasBR){
    fBRtr->SetBranchAddress("timestamp",&flocalBRts);
    fBRtr->SetBranchAddress("checkADC",&flocalcheckADC);
    fBRtr->SetBranchAddress("trigbit",&flocaltrigbit);
    fBRtr->SetBranchAddress("beam",&flocalbeam);
    if(fSett->BigRIPSDetail()>0){
      for(unsigned short f=0;f<NFPLANES;f++){
	fBRtr->SetBranchAddress(Form("fp%d",fpID[f]),&flocalfp[f]);
      }
    }
    fBRentries = fBRtr->GetEntries();
    cout << fBRentries << " entries in BigRIPS tree" << endl;
  }
  if(fhasHI){
    fHItr->SetBranchAddress("hicaricalc",&flocalhicari);
    fHIentries = fHItr->GetEntries();
    cout << fHIentries << " entries in HiCARI tree" << endl;
  }
  if(fhasM2){
    fM2tr->SetBranchAddress("gretinacalc",&flocalmode2);
    fM2entries = fM2tr->GetEntries();
    cout << fM2entries << " entries in Mode2 tree" << endl;
  }
    
  fmtr = new TTree("tr","merged tree");
  fmtr->Branch("beam",&fbeam,320000);
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      fmtr->Branch(Form("fp%d",fpID[f]),&ffp[f],320000);
    }
  }
  fmtr->Branch("checkADC",&fcheckADC,320000);
  fmtr->Branch("trigbit",&ftrigbit,320000);
  fmtr->Branch("brTS",&fBRts,320000);
  fmtr->Branch("brentry",&fBRentry,320000);
  fmtr->Branch("hicari",&fhicari,320000);
  fmtr->Branch("hiTS",&fHIts,320000);
  fmtr->Branch("hientry",&fHIentry,320000);
  fmtr->Branch("mode2",&fmode2,320000);
  fmtr->Branch("m2TS",&fM2ts,320000);
  fmtr->Branch("m2entry",&fM2entry,320000);
  fmtr->BranchRef();


  fcurrentts = 0;
  Int_t status = 0;
  if(fhasBR){
    status = fBRtr->GetEvent(0);
    if(status<0)
      cout << "first BigRIPS entry faulty!" << endl;
    fcurrentts = flocalBRts;
    cout << "first BigRIPS timestamp: " << flocalBRts << endl;
    status = fBRtr->GetEvent(fBRentries-1);
    if(status<0)
      cout << "last BigRIPS entry faulty!" << endl;
    cout << "last BigRIPS timestamp: " << flocalBRts << endl;
  }
  if(fhasHI){
    status = fHItr->GetEvent(0);
    if(status<0)
      cout << "first HiCARI entry faulty!" << endl;
    if(fcurrentts <1 || flocalhicari->GetTS() <  fcurrentts)
      fcurrentts = flocalhicari->GetTS();
    cout << "first HiCARI timestamp: " << flocalhicari->GetTS() << endl;
    status = fHItr->GetEvent(fHIentries-1);
    if(status<0)
      cout << "last HiCARI entry faulty!" << endl;
    cout << "last HiCARI timestamp: " << flocalhicari->GetTS() << endl;;
  }
  if(fhasM2){
    status = fM2tr->GetEvent(0);
    if(status<0)
      cout << "first Mode2 entry faulty!" << endl;
    if(fcurrentts <1 || flocalmode2->GetTS() <  fcurrentts)
      fcurrentts = flocalmode2->GetTS();
    cout << "first Mode2 timestamp: " << flocalmode2->GetTS() << endl;
    status = fM2tr->GetEvent(fM2entries-1);
    if(status<0)
      cout << "last Mode2 entry faulty!" << endl;
    cout << "last Mode2 timestamp: " << flocalmode2->GetTS() << endl;;
  }
  if(fcurrentts<1){
    cout << RED << "cannot find first timestamp, Aborting." << DEFCOLOR << endl;
    return;
  }
  
  flocalBRts = 0;
  flocalHIts = 0;
  flocalM2ts = 0;
  flocalcheckADC = -1;
  flocaltrigbit = -1;
  flocalbeam->Clear();
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      flocalfp[f]->Clear();
    }
  }
  flocalhicari->Clear();
  flocalmode2->Clear();

  flastBRts = 0;
  flastHIts = 0;
  flastM2ts = 0;
  fBRtsjump = false;
  fHItsjump = false;
  fM2tsjump = false;


  fmhist = new MergeHistograms(fSett);
}

bool BuildEvents::ReadBigRIPS(){
  if(fverbose>1)
    cout << __PRETTY_FUNCTION__ << endl;
  
  flocalcheckADC = -1;
  flocaltrigbit = -1;
  flocalbeam->Clear();
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      flocalfp[f]->Clear();
    }
  }
  flocalBRts = 0;
  if(fBRentry==fBRentries){
    return false;
  }
  Int_t status = fBRtr->GetEvent(fBRentry);
  if(fverbose>2)
    cout << "status " << status << endl;
  if(status == -1){
    cerr<<"Error occured, couldn't read entry "<<fBRentry<<" from tree "<<fBRtr->GetName()<<endl;
    return false;
  }
  else if(status == 0){
    cerr<<"Error occured, entry "<<fBRentry<<" in tree "<<fBRtr->GetName()<<" in file doesn't exist"<<endl;
    return false;
  }
  fnbytes += status;
  if(flocalBRts<flastBRts){
    cout << endl << "BigRIPS timestamp jump detected. this = " << flocalBRts << ", last = " << flastBRts << endl;
    fBRtsjump = true;
    return false;
  }

  if(fverbose>0)
    cout << "read new bigrips with TS = " << flocalBRts << " tof = "<< flocalbeam->GetTOF(0)+48<< endl;
  fBRentry++;

  detector* det = new detector;
  det->TS = flocalBRts;
  det->ID = 0;
  fdetectors.push_back(det);

  flastBRts = flocalBRts;
  return true;
}

bool BuildEvents::ReadHiCARI(){
  if(fverbose>1)
    cout << __PRETTY_FUNCTION__ << endl;
  flocalhicari->Clear();
  flocalHIts = 0;
  if(fHIentry==fHIentries){
    return false;
  }
  Int_t status = fHItr->GetEvent(fHIentry);
  if(fverbose>2)
    cout << "status " << status << endl;
  if(status == -1){
    cerr<<"Error occured, couldn't read entry "<<fHIentry<<" from tree "<<fHItr->GetName()<<endl;
    return false;
  }
  else if(status == 0){
    cerr<<"Error occured, entry "<<fHIentry<<" in tree "<<fHItr->GetName()<<" in file doesn't exist"<<endl;
    return false;
  }
  fnbytes += status;
  
  flocalHIts = flocalhicari->GetTS();

  if(flocalHIts<flastHIts){
    cout <<"HiCARI timestamp jump detected. this = " << flocalHIts << ", last = " << flastHIts << endl;
    fHItsjump = true;
    return false;
  }
  if(fverbose>0)
    cout << "read new hicari with TS = " << flocalHIts << endl;
  fHIentry++;

  detector* det = new detector;
  det->TS = flocalHIts;
  det->ID = 1;
  fdetectors.push_back(det);

  flastHIts = flocalHIts;
  return true;
}

bool BuildEvents::ReadMode2(){
  if(fverbose>1)
    cout << __PRETTY_FUNCTION__ << endl;
  flocalmode2->Clear();
  flocalM2ts = 0;
  if(fM2entry==fM2entries){
    return false;
  }
  Int_t status = fM2tr->GetEvent(fM2entry);
  if(fverbose>2)
    cout << "status " << status << endl;
  if(status == -1){
    cerr<<"Error occured, couldn't read entry "<<fM2entry<<" from tree "<<fM2tr->GetName()<<endl;
    return false;
  }
  else if(status == 0){
    cerr<<"Error occured, entry "<<fM2entry<<" in tree "<<fM2tr->GetName()<<" in file doesn't exist"<<endl;
    return false;
  }
  fnbytes += status;
  flocalM2ts = flocalmode2->GetTS();

  if(flocalM2ts<flastM2ts){
    cout <<"Mode2 timestamp jump detected. this = " << flocalM2ts << ", last = " << flastM2ts << endl;
    fM2tsjump = true;
    return false;
  }
  if(fverbose>0)
    cout << "read new mode2 with TS = " << flocalM2ts << endl;
  fM2entry++;

  detector* det = new detector;
  det->TS = flocalM2ts;
  det->ID = 2;
  fdetectors.push_back(det);

  flastM2ts = flocalM2ts;
  return true;
}

bool BuildEvents::ReadEach(){
  flastBRts = 0;
  flastHIts = 0;
  flastM2ts = 0;
  bool success = false;
  if(fhasHI)
    success += ReadHiCARI();
  if(fhasM2)
    success += ReadMode2();
  if(fhasBR)
    success += ReadBigRIPS();
  if(success)
    return true;
  else{
    cout << endl << "failed to read from HiCARI, Mode2, and BigRIPS." << endl;
    return false;
  }
}

void BuildEvents::CloseEvent(){
  if(fverbose>0)
    cout << __PRETTY_FUNCTION__ << endl;
  // // bool printme = false;
  // if(fBRts>0 && fHIts >0){
  //   //printme = true;
  //   cout << "fBRentry = " << fBRentry << "\tfBRts = " << fBRts << ", fHIentry = " << fHIentry << "\tfHIts = " << fHIts << "\tmult = " << fhicari->GetMult()<< endl;
  //   cout << "BR AoQ = " << fbeam->GetAQ(1) << " Z = " << fbeam->GetZ(1) << endl;
  //   fhicari->Print();
  // }
  if(fverbose>0 && fBRts>0){
    cout << "closing event with local TS = " << flocalBRts << " tof = "<< flocalbeam->GetTOF(0)<< endl;
    cout << "closing event with set TS = " << fBRts << " tof = "<< fbeam->GetTOF(0)<< endl;
  }
  if(fhicari->HadBigRIPS()){
    //cout << fhicari->GetBigRIPSHit()->GetTS() << endl;
    long long int br_TS = fhicari->GetBigRIPSHit()->GetTS();
    fhicari->CorrectTime(br_TS);
    fmode2->CorrectTime(br_TS);
  }

  
  switch(fmode){
  default:
  case 0: //write all events
    fmtr->Fill();
    fmhist->FillHistograms(fcheckADC, fhicari, fmode2, fBRts, fHIts, fM2ts);
    break;
  case 1://BR and HI coincidence
    if(fBRts>0 && fHIts >0){
      fmtr->Fill();
      fmhist->FillHistograms(fcheckADC, fhicari, fmode2, fBRts, fHIts, fM2ts);
    }
    break;
  }

  fBRts = 0;
  fcheckADC = -1;
  ftrigbit = -1;
  fbeam->Clear();
  if(fSett->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      ffp[f]->Clear();
    }
  }
  fHIts = 0;
  fhicari->Clear();
  fM2ts = 0;
  fmode2->Clear();
  // if(printme){
  //   cout << "after clearing " << endl;
  //   fhicari->Print();
  // }

  if(fverbose>0)
    cout << "end "<< __PRETTY_FUNCTION__ << endl;
}

bool BuildEvents::Merge(){
  if(fverbose>1)
    cout << __PRETTY_FUNCTION__ << endl;

  if(fverbose>1){
    cout << "before sorting" << endl;
    for(vector<detector*>::iterator det=fdetectors.begin(); det!=fdetectors.end(); det++){
      cout << "ID = " << (*det)->ID << ", TS = " << (*det)->TS << endl;
    }
  }
  sort(fdetectors.begin(), fdetectors.end(),TSComparer());
  if(fverbose>1){
    cout << "after sorting" << endl;
    for(vector<detector*>::iterator det=fdetectors.begin(); det!=fdetectors.end(); det++){
      cout << "ID = " << (*det)->ID << ", TS = " << (*det)->TS << endl;
    }
  }

  switch(fdetectors.at(0)->ID){
  case 0: //BigRIPS
    if(fBRts>0){
      if(fverbose>1)
	cout << "has already BigRIPS" << endl;
      CloseEvent();
    }
    else if(flocalBRts - fcurrentts > fwindow){
      if(fverbose>0)
	cout << "BR larger than window"  << flocalBRts <<" - "<<fcurrentts<< " = " << flocalBRts - fcurrentts << " > "<<fwindow << endl;
      CloseEvent();
    }
    fBRts = flocalBRts;
    fcheckADC = flocalcheckADC;
    ftrigbit = flocaltrigbit;
    fbeam = (Beam*)flocalbeam->Clone();
    if(fSett->BigRIPSDetail()>0){
      for(unsigned short f=0;f<NFPLANES;f++){
	ffp[f] = (FocalPlane*)flocalfp[f]->Clone(Form("fp_%d",f));
      }
    }
    fcurrentts = fBRts;
    fdetectors.erase(fdetectors.begin());
    if(!ReadBigRIPS()&&fBRtsjump==false)
      cout << endl << MAGENTA << "failed to read BigRIPS, end of file" << DEFCOLOR << endl;
    break;
  case 1: //HiCARI
    if(flocalHIts - fcurrentts > fwindow){
      if(fverbose>0)
	cout << "HI larger than window " << flocalHIts <<" - "<<fcurrentts<< " = " << flocalHIts - fcurrentts << " > "<<fwindow << endl;
      CloseEvent();
    }
    fHIts = flocalHIts;
    fhicari = (HiCARICalc*)flocalhicari->Clone();
    flocalhicari->Clear();
    fcurrentts = fHIts;
    fdetectors.erase(fdetectors.begin());
    if(!ReadHiCARI()&&fHItsjump==false)
      cout << endl << MAGENTA << "failed to read HiCARI, end of file" << DEFCOLOR << endl;
    break;
  case 2: //Mode2
    if(flocalM2ts - fcurrentts > fwindow){
      if(fverbose>0)
	cout << "M2 larger than window" << endl;
      CloseEvent();
    }
    fM2ts = flocalM2ts;
    fmode2 = (GretinaCalc*)flocalmode2->Clone();
    flocalmode2->Clear();
    fcurrentts = fM2ts;
    fdetectors.erase(fdetectors.begin());
    if(!ReadMode2()&&fM2tsjump==false)
      cout << endl << MAGENTA << "failed to read Mode2, end of file" << DEFCOLOR << endl;
    break;
  default:
    break;
  }
  if(flastevent>0 && flastevent == (int)(fBRentry + fHIentry + fM2entry)){
    cout << endl << BLUE << "last event reached " << DEFCOLOR << endl;
    return false;
  }

  if(fdetectors.size()==0){
    if(fhasBR && fhasHI && fhasM2 && fBRtsjump==true && fHItsjump==true && fM2tsjump==true){
      cout << RED << "all timestamps jumped" << DEFCOLOR << endl;
      fBRtsjump = false;
      fHItsjump = false;
      fM2tsjump = false;
      return ReadEach();
    }
    cout << endl << BLUE << "all files finished " << DEFCOLOR << endl;
    return false;
  }
  return true;
}
