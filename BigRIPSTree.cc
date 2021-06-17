#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtEventInfo.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibIC.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtIC.hh"
#include "TArtPPAC.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtFocalPlane.hh"
#include "TArtTOF.hh"
#include "TArtRIPS.hh"
#include "TArtBeam.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "CommandLineInterface.hh"
#include "Settings.hh"
#include "RunInfo.hh"

#include "FocalPlane.hh"
#include "PPAC.hh"
#include "Beam.hh"
#include "Globaldefs.h"
using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig);
double get_time();
int main(int argc, char* argv[]){
  double time_start = get_time();  
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);
  char* InputFile = NULL;
  char* OutputFile = NULL;
  char* SetFile = NULL;
  int nmax = 0;
  int vl = 0;
  int detail = -1;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &OutputFile);
  interface->Add("-s", "settings file", &SetFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->Add("-d", "detail in the tree", &detail);
  interface->CheckFlags(argc, argv);

  if(InputFile == NULL || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  cout<<"input file:"<<InputFile<<endl;
  //get the run number from the filename

  string infilestr = string(InputFile);
  char extension[20];
  strcpy(extension, infilestr.substr(infilestr.find_last_of(".")+1).c_str());
  int run;
  TString ifname(InputFile);
  if(strcmp(extension,"gz")==0)
    ifname.Remove(0,ifname.Length()-12); // Last 9 characters: nameXXXX.ridf.gz
  else
    ifname.Remove(0,ifname.Length()-9); // Last 9 characters: nameXXXX.ridf
  //cout << "ifname.Data() " << ifname.Data() << endl;
  sscanf(ifname.Data(),"%04d.ridf",&run);
  //cout << run << endl;

  cout<<"output file: "<<OutputFile<< endl;
  if(SetFile == NULL)
    cerr<<"No settings file! Using standard values"<<endl;
  else
    cout<<"settings file: "<<SetFile<<endl;

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
    
  if(outfile->IsZombie()){
    return 4;
  }
  Settings* set = new Settings(SetFile);
  set->SetVLevel(vl);
  if(detail>-1)
    set->SetBigRIPSDetail(detail);
  cout << "BigRIPSTree Detail Level: " << set->BigRIPSDetail() << endl;
  
  RunInfo* info = new RunInfo();
  info->SetBRRunNumber(run);
  cout << "run number: "<< run << endl;
  set->Write("settings",TObject::kOverwrite);


  TEnv *evtnumbers = new TEnv(set->EvtNrFile());
  int startevent = evtnumbers->GetValue(Form("Start.Event.Number.%d",run),0);
  cout << "run number " << run << " starts with event nr " << startevent << endl;

  /// load time dependent corrections
  TFile *cor = new TFile(set->TimeCorFile());
  TH1F* f7_cor[2];
  TH1F* f11_cor[2];
  TH1F* br_cor[2];
  TH1F* zd_cor[2];
  if(cor->IsZombie()){
    cerr << "ignore previous warning!!" << endl;
    cerr << "File " << set->TimeCorFile() << " not existing, cannot perform time dependent corrections!" << endl;
    f7_cor[0] = NULL;
    f7_cor[1] = NULL;
    f11_cor[0] = NULL;
    f11_cor[1] = NULL;
    br_cor[0] = NULL;
    br_cor[1] = NULL;
    zd_cor[0] = NULL;
    zd_cor[1] = NULL;
  }
  else{
    f7_cor[0] = (TH1F*)cor->Get("hoffsF7");
    f7_cor[1] = (TH1F*)cor->Get("hgainF7");
    f11_cor[0] = (TH1F*)cor->Get("hoffsF11");
    f11_cor[1] = (TH1F*)cor->Get("hgainF11");
    br_cor[0] = (TH1F*)cor->Get("hoffsBR");
    br_cor[1] = (TH1F*)cor->Get("hgainBR");
    zd_cor[0] = (TH1F*)cor->Get("hoffsZD");
    zd_cor[1] = (TH1F*)cor->Get("hgainZD");

    outfile->cd();
    for(int i=0;i<2;i++){
      if(f7_cor[i]!=NULL)
	f7_cor[i]->Write(); 
      if(f11_cor[i]!=NULL)
	f11_cor[i]->Write(); 
      if(br_cor[i]!=NULL)
	br_cor[i]->Write(); 
      if(zd_cor[i]!=NULL)
	zd_cor[i]->Write(); 
    }
  }
  outfile->cd();

  
  TArtStoreManager* sman = TArtStoreManager::Instance();

  TArtEventStore* estore = new TArtEventStore();
  estore->SetInterrupt(&signal_received);
  estore->Open(InputFile);
  //std::cout<<"estore ->"<< InputFile <<std::endl;

  TArtBigRIPSParameters* para = TArtBigRIPSParameters::Instance();
  para->LoadParameter(set->PPACFile());
  para->LoadParameter(set->PlasticFile());
  para->LoadParameter(set->ICFile());
  para->LoadParameter(set->FocalFile());
  TArtCalibPID* brcalib          = new TArtCalibPID();
  TArtCalibPPAC* ppaccalib       = brcalib->GetCalibPPAC();
  TArtCalibPlastic* plasticcalib = brcalib->GetCalibPlastic();
  TArtCalibIC* iccalib           = brcalib->GetCalibIC(); 
  TArtCalibFocalPlane* cfpl      = brcalib->GetCalibFocalPlane();

  TArtRecoPID* recopid = new TArtRecoPID();
 
  // Definition of observables we want to reconstruct
  cout << "Defining bigrips parameters" << endl;

  TArtRIPS* recorips[4];
  recorips[0] = recopid->DefineNewRIPS(3,5,set->MatrixFile(0),"D3"); // F3 - F5
  recorips[1] = recopid->DefineNewRIPS(5,7,set->MatrixFile(1),"D5"); // F5 - F7
  recorips[2] = recopid->DefineNewRIPS(8,9,set->MatrixFile(2),"D7"); // F8 - F9
  recorips[3] = recopid->DefineNewRIPS(9,11,set->MatrixFile(3),"D8"); // F9 - F11  

  // Reconstruction of TOF DefineNewTOF(first plane, second plane, time offset)
  TArtTOF* recotof[6];
  recotof[0] = recopid->DefineNewTOF("F3pl","F7pl",set->TimeOffset(0),5); // F3 - F5
  recotof[1] = recopid->DefineNewTOF("F3pl","F7pl",set->TimeOffset(1),5); // F5 - F7
  recotof[2] = recopid->DefineNewTOF("F3pl","F7pl",set->TimeOffset(2),5); // F3 - F7

  recotof[3] = recopid->DefineNewTOF("F8pl","F11pl-1",set->TimeOffset(3),9); // F8 - F9
  recotof[4] = recopid->DefineNewTOF("F8pl","F11pl-1",set->TimeOffset(4),9); // F9 - F11
  recotof[5] = recopid->DefineNewTOF("F8pl","F11pl-1",set->TimeOffset(5),9); // F8 - F11
  

  TArtTOF* tof7to8  = recopid->DefineNewTOF("F7pl","F8pl"); // F7 - F8

  // Reconstruction of IC observables for ID
  TArtBeam* recobeam[6];
  recobeam[0] = recopid->DefineNewBeam(recorips[0],recotof[0],"F7IC");
  recobeam[1] = recopid->DefineNewBeam(recorips[1],recotof[1],"F7IC");   
  recobeam[2] = recopid->DefineNewBeam(recorips[0],recorips[1],recotof[2],"F7IC");   

  recobeam[3] = recopid->DefineNewBeam(recorips[2],recotof[3],"F11IC");
  recobeam[4] = recopid->DefineNewBeam(recorips[3],recotof[4],"F11IC");
  recobeam[5] = recopid->DefineNewBeam(recorips[2],recorips[3],recotof[5],"F11IC");


  // output tree
  TTree* tr = new TTree("BRZDtr","BigRIPS / ZeroDeg Data Tree");
  //branch for trig bit
  int trigbit = 0;
  tr->Branch("trigbit",&trigbit,"trigbit/I");
  int checkADC = 0;
  tr->Branch("checkADC",&checkADC,"checkADC/I");
  //branch for original event number
  int eventnumber = 0;
  tr->Branch("eventnumber",&eventnumber,"eventnumber/I");
  int toteventnumber = 0;
  tr->Branch("toteventnumber",&toteventnumber,"toteventnumber/I");
  //branch for timestamp
  unsigned long long int timestamp = 0;
  tr->Branch("timestamp",&timestamp,"timestamp/l");
  //branches for each focal plane
  FocalPlane* fp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    fp[f] = new FocalPlane;
    if(set->BigRIPSDetail()>0)
      tr->Branch(Form("fp%d",fpID[f]),&fp[f],320000);
  }
  //branch for the beam, beta, a/q, z
  Beam* beam = new Beam;
  tr->Branch("beam",&beam,320000);
  //branch for the PPACs
  PPAC* ppacs = new PPAC;
  if(set->BigRIPSDetail()>1)
    tr->Branch("ppacs",&ppacs,320000);


  int event = startevent;
  unsigned long long int last_timestamp = 0;
  unsigned long long int rawevent_timestamp = 0;
  int ctr =0;
  while(estore->GetNextEvent() && !signal_received){
    //clearing
    trigbit = 0;
    checkADC = -1;
    timestamp = 0;
    eventnumber++;
    toteventnumber = event;
    for(int f=0;f<NFPLANES;f++){
      fp[f]->Clear();
    }
    ppacs->Clear();
    beam->Clear();

    //Making the BigRIPS tree calibration
    brcalib->ClearData();
    brcalib->ReconstructData();
    

    //trigger bit information
    TArtRawEventObject* rawevent = (TArtRawEventObject*)sman->FindDataContainer("RawEvent");
    int rawevent_number = rawevent->GetEventNumber();
    rawevent_timestamp = rawevent->GetTimeStamp();
    for(int i=0; i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject* seg = rawevent->GetSegment(i);
      Int_t fpl = seg->GetFP();
      Int_t detector = seg->GetDetector();
      //cout << fpl << "\t " << detector << endl;
      if(fpl==8 && detector ==21){
	if(seg->GetNumData()<1)
	  break;
	TArtRawDataObject* d = seg->GetData(set->CorrelationChannel());
	checkADC = d->GetVal();
	if(set->VLevel()>1){
	  for(int j=0; j < seg->GetNumData(); j++){
	    TArtRawDataObject* d = seg->GetData(j);
	    cout << j << "\t" << d->GetVal() << endl;
	  }
	}//vlevel
      }
      if(fpl==63 && detector==10){
        for(int j=0; j < seg->GetNumData(); j++){
	  TArtRawDataObject* d = seg->GetData(j);
	  Short_t ch = d->GetCh();
	  if(ch==0)
	    trigbit = d->GetVal();
        }
      }
      
    }
    //timestamp information
    TClonesArray* info_a = (TClonesArray*)sman->FindDataContainer("EventInfo");
    TArtEventInfo* info = (TArtEventInfo*)info_a->At(0);
    unsigned int bit = info->GetTriggerBit();
    timestamp = info->GetTimeStamp();
    if(timestamp<last_timestamp){
      cout << "timestamp was reset, this TS = " << timestamp << ", last one was " << last_timestamp << " difference " << timestamp-last_timestamp << endl;
    }
    //
    if(vl>1)
      cout << "tb = "<< trigbit << "\t bit = " << bit << "\t TS(event info) = " << timestamp << "\t diff to last" << timestamp-last_timestamp << "\t eventnumber = "  << eventnumber << "\t rawevent eventnumber = " << rawevent_number << "\t rawevent timestamp = " << rawevent_timestamp << endl;
    last_timestamp = timestamp;


    TArtPPAC* tppac;
    for(unsigned short p=0;p<NPPACS;p++){
      SinglePPAC* dppac = new SinglePPAC;
      tppac = ppaccalib->GetPPAC(p);
      if(tppac){
	dppac->Set(tppac->GetID(),tppac->GetX(),tppac->GetY(),tppac->GetXZPos(),tppac->GetYZPos(),tppac->GetTSumX(),tppac->GetTSumY());
	if(tppac->IsFiredX()||tppac->IsFiredY())
	  ppacs->AddPPAC(dppac);
      }
    }


    //focal plane detector information
    TArtFocalPlane* tfpl;
    TArtPlastic* tpla;
    TArtIC* tic;
    TVectorD* vec;
    Track track;
    Plastic plastic;
    MUSIC music;

    for(unsigned short f=0;f<NFPLANES;f++){
      
      fp[f]->Clear();
      track.Clear();
      tfpl = cfpl->FindFocalPlane(fpID[f]);

      TMatrixD xvec(2,1); xvec.Zero();
      TMatrixD yvec(2,1); yvec.Zero();
      TMatrixD xmat(2,2); xmat.Zero();
      TMatrixD ymat(2,2); ymat.Zero();
      int first = firstPPAC(fpID[f]);
      if(first>-1){
   
	double zpos = tfpl->GetZoffset();
	int nfiredx[3] = {0, 0, 0}; //total, upstream, downstream
	int nfiredy[3] = {0, 0, 0};
	for(unsigned short p=0;p<4;p++){
	  SinglePPAC* dppac = ppacs->GetPPACID(first+p);
	  double x = dppac->GetX();
	  double y = dppac->GetY();
	  double zx = dppac->GetXZ() - zpos;
	  double zy = dppac->GetYZ() - zpos;
	  if(dppac->FiredX()){
	    xvec(0,0) += zx*x;
	    xvec(1,0) += x;
	    xmat(0,1) += zx;
	    xmat(1,0) += zx;
	    xmat(0,0) += zx*zx;
	    xmat(1,1) ++;
	    nfiredx[0]++;
	    if(p<2)
	      nfiredx[1]++;
	    else
	      nfiredx[2]++;
	  }
	  if(dppac->FiredY()){
	    yvec(0,0) += zy*y;
	    yvec(1,0) += y;
	    ymat(0,1) += zy;
	    ymat(1,0) += zy;
	    ymat(0,0) += zy*zy;
	    ymat(1,1) ++;
	    nfiredy[0]++;
	    if(p<2)
	      nfiredy[1]++;
	    else
	      nfiredy[2]++;
	  }
	}
	if(nfiredx[1]>0 && nfiredx[2]>0){
	  TMatrixD rxvec = xmat.Invert()*xvec;
	  tfpl->SetOptVector(0,rxvec(1,0));
	  tfpl->SetOptVector(1,TMath::ATan(rxvec(0,0))*1000);
	}
	else{
	  tfpl->SetOptVector(0,-99999);
	  tfpl->SetOptVector(1,-99999);
	}
	if(nfiredy[1]>0 && nfiredy[2]>0){
	  TMatrixD ryvec = ymat.Invert()*yvec;
	  tfpl->SetOptVector(2,ryvec(1,0));
	  tfpl->SetOptVector(3,TMath::ATan(ryvec(0,0))*1000);
	}
	else{
	  tfpl->SetOptVector(2,-99999);
	  tfpl->SetOptVector(3,-99999);
	}
	tfpl->SetNumFiredPPACX(nfiredx[0]);
	tfpl->SetNumFiredPPACY(nfiredy[0]);
	tfpl->CopyPos();
 
	if(vl>2)
	  cout << "FP " << fpID[f] ;
	if(tfpl){
	  if(vl>2)
	    cout << "\tnfiredx  = " << tfpl->GetNumFiredPPACX()<< ", nfiredy  = " << tfpl->GetNumFiredPPACY();
	  vec=tfpl->GetOptVector(); 
	  if(vl>2)
	    cout << "\tx = " << (*vec)(0) <<", y = "<< (*vec)(2)<<", a = "<< (*vec)(1) <<",b = " << (*vec)(3);
	  track.Set((*vec)(0), (*vec)(2), (*vec)(1), (*vec)(3));
	}
	if(vl>2)
	  cout << endl;
      }//ppac exists


      // PLASTICS
      plastic.Clear();
      tpla = plasticcalib->FindPlastic(Form("F%dpl",fpID[f]));
      if(fpID[f]==11)
	tpla = plasticcalib->FindPlastic(Form("F%dpl-1",fpID[f]));
      if(tpla){
	//cout << f << "\t" << fpID[f] << "\t";
	//cout << tpla->GetTimeL() <<", "<< tpla->GetTimeR() <<", "<< tpla->GetQLRaw() <<", "<< tpla->GetQRRaw() << endl;
	plastic.SetTime(tpla->GetTimeL(), tpla->GetTimeR());
	plastic.SetCharge(tpla->GetQLRaw(), tpla->GetQRRaw());
      }


      // MUSICS
      music.Clear();
      tic = iccalib->FindIC(Form("F%dIC",fpID[f]));
      if(tic){
	music.SetNHits(tic->GetNumHit());
	music.SetEnergy(tic->GetEnergyAvSum(),tic->GetEnergySqSum());
	if(set->BigRIPSDetail()>1){
	  //raw ADC values
	  for(int c=0;c<NUM_IC_CHANNEL;c++){
	    double adc = tic->GetRawADC(c);
	    if(adc>0){
	      music.SetRawADC(adc, c);
	      music.SetGainMatchADC(tic->GetGainMatchADC(c));
	    }
	  }// ic channels
	}// BR detail
      }// IC present

      fp[f]->SetTrack(track);
      fp[f]->SetPlastic(plastic);
      fp[f]->SetMUSIC(music);
    }
    
    //Reconstructiong the PID
    recopid->ClearData();
    recopid->ReconstructData();

    //beam parameters and PID
    beam->SetTOFBeta(0,recotof[2]->GetTOF(),recotof[2]->GetBeta());
    beam->SetTOFBeta(1,recotof[5]->GetTOF(),recotof[5]->GetBeta());
    beam->SetTOFBeta(2,tof7to8->GetTOF(),tof7to8->GetBeta());
    
    for(unsigned short b=0;b<6;b++){
      double z = recobeam[b]->GetZet();
      double gz = 1;
      double oz = 0;
      if(b<3 && f7_cor[0]!=NULL && f7_cor[1]!=NULL){
	gz = f7_cor[1]->GetBinContent(event/10000+1);
	oz = f7_cor[0]->GetBinContent(event/10000+1);
      }
      if(b>2 && f11_cor[0]!=NULL && f11_cor[1]!=NULL){
	gz = f11_cor[1]->GetBinContent(event/10000+1);
	oz = f11_cor[0]->GetBinContent(event/10000+1);
      }
      double a = recobeam[b]->GetAoQ();
      double ga = 1;
      double oa = 0;
      if(b<3 && br_cor[0]!=NULL && br_cor[1]!=NULL){
	ga = br_cor[1]->GetBinContent(event/10000+1);
	oa = br_cor[0]->GetBinContent(event/10000+1);
      }
      if(b>2 && zd_cor[0]!=NULL && zd_cor[1]!=NULL){
	ga = zd_cor[1]->GetBinContent(event/10000+1);
	oa = zd_cor[0]->GetBinContent(event/10000+1);
      }
      z = z*gz+oz;
      a = a*ga+oa;

      //cout << z << "\t" << event << "\t" << o << "\t" << g << "\t" << z*g+o << endl;
      
      beam->SetAQZ(b,a,z);
      beam->SetRIPSBeta(b,recobeam[b]->GetBeta());
    }
 
    for(unsigned short b=0;b<3;b++){
      double corr = 0;
      corr = set->GetBRAoQCorrection_F3X()*fp[fpNr(3)]->GetTrack()->GetX() +
	     set->GetBRAoQCorrection_F3A()*fp[fpNr(3)]->GetTrack()->GetA() +
	     set->GetBRAoQCorrection_F3Q()*sqrt(fp[fpNr(3)]->GetPlastic()->GetChargeL() * fp[fpNr(3)]->GetPlastic()->GetChargeR()) +
	     set->GetBRAoQCorrection_F5X()*fp[fpNr(5)]->GetTrack()->GetX() +
	     set->GetBRAoQCorrection_F5A()*fp[fpNr(5)]->GetTrack()->GetA() +
	     set->GetBRAoQCorrection_F7X()*fp[fpNr(7)]->GetTrack()->GetX() +
	     set->GetBRAoQCorrection_F7A()*fp[fpNr(7)]->GetTrack()->GetA() +
	     set->GetBRAoQCorrection_F7Q()*sqrt(fp[fpNr(7)]->GetPlastic()->GetChargeL() * fp[fpNr(7)]->GetPlastic()->GetChargeR());      
      beam->CorrectAQ(b,corr);
      beam->ScaleAQ(b,set->GetBRAoQCorrection_gain(),set->GetBRAoQCorrection_offs());
      
      corr = set->GetZDAoQCorrection_F8X()*fp[fpNr(8)]->GetTrack()->GetX() +
	     set->GetZDAoQCorrection_F8A()*fp[fpNr(8)]->GetTrack()->GetA() +
	     set->GetZDAoQCorrection_F8Q()*sqrt(fp[fpNr(8)]->GetPlastic()->GetChargeL() * fp[fpNr(8)]->GetPlastic()->GetChargeR()) +
	     set->GetZDAoQCorrection_F9X()*fp[fpNr(9)]->GetTrack()->GetX() +
	     set->GetZDAoQCorrection_F9A()*fp[fpNr(9)]->GetTrack()->GetA() +
	     set->GetZDAoQCorrection_F11X()*fp[fpNr(11)]->GetTrack()->GetX() +
	     set->GetZDAoQCorrection_F11A()*fp[fpNr(11)]->GetTrack()->GetA() +
	     set->GetZDAoQCorrection_F11Q()*sqrt(fp[fpNr(11)]->GetPlastic()->GetChargeL() * fp[fpNr(11)]->GetPlastic()->GetChargeR()); 
      beam->CorrectAQ(b+3,corr);
      beam->ScaleAQ(b+3,set->GetZDAoQCorrection_gain(),set->GetZDAoQCorrection_offs());
    }
    
    for(unsigned short b=0;b<4;b++){
      beam->SetDelta(b ,recorips[b]->GetDelta());
      beam->SetBrho(b ,recorips[b]->GetBrho());
    }

    //F8 ppacs
    beam->SetF8PPAC1A((*ppacs->GetPPACID(19)));
    beam->SetF8PPAC1B((*ppacs->GetPPACID(20)));
    beam->SetF8PPAC2A((*ppacs->GetPPACID(21)));
    beam->SetF8PPAC2B((*ppacs->GetPPACID(22)));
    beam->SetF8PPAC3A((*ppacs->GetPPACID(35)));
    beam->SetF8PPAC3B((*ppacs->GetPPACID(36)));


    // // calculate Z directly:
    // TArtIC* F7ic = iccalib->FindIC("F7IC");
    // double adc[6];
    // for(int i=0;i<6;i++){
    //   adc[i] = F7ic->GetRawADC(i);
    // }
    
    
    //fill the tree
    tr->Fill();

    //output
    if(ctr%1000 == 0){
      double time_end = get_time();
      cout << setw(5) << ctr << " events done " << setiosflags(ios::fixed) << setprecision(1) << (Float_t)ctr/(time_end - time_start) << " events/s \r" << flush;
      if(ctr%100000 == 0)
	tr->AutoSave();
    }
    ctr++;
    event++;
    
    if(nmax>0 && ctr>nmax-1)
      break;
  }

  double time_end = get_time();
  cout << "Program Run time: " << time_end - time_start << " s." << endl;
  cout << "Total of " << ctr << " events processed, " << ctr/(time_end - time_start) << " events/s." << endl;
  cout << BLUE << tr->GetEntries() << DEFCOLOR << " entries written to tree ("<<BLUE<<tr->GetZipBytes()/(1024*1024)<< DEFCOLOR<<" MB)"<< endl;
  info->SetBREvents(tr->GetEntries());
  tr->Write("",TObject::kOverwrite);
  info->Write("info",TObject::kOverwrite);
  outfile->Close();
  timer.Stop();
  cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

  return 0;
}
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){  
    struct timeval t;  
    gettimeofday(&t, NULL);  
    double d = t.tv_sec + (double) t.tv_usec/1000000;  
    return d;  
}  
