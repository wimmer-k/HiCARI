#include "SimHistograms.hh"
using namespace TMath;
using namespace std;

#ifdef USELISA
void SimHistograms::FillHistograms(GretinaCalc* gr, MiniballCalc* mb, ZeroDeg* zd, GammaSim* gs, LISA* li){
#else
void SimHistograms::FillHistograms(GretinaCalc* gr, MiniballCalc* mb, ZeroDeg* zd, GammaSim* gs){
#endif

  //-------------------------------------------------------------------------
  //*************************************************************************
  //Fill the histograms here.
  //*************************************************************************
  //-------------------------------------------------------------------------
  vector<double> emitted_energies;
  double highestemitted =0;
  for(UShort_t g = 0; g < gs->GetMult(); g++){
    if(gs->GetEmittedGamma(g)->GetEnergy()>highestemitted)
      highestemitted = gs->GetEmittedGamma(g)->GetEnergy();
    Fill("egam_emitted",8000,0,8000,gs->GetEmittedGamma(g)->GetEnergy());
    Fill("posx_emitted",100,-5,5,gs->GetEmittedGamma(g)->GetX());
    Fill("posy_emitted",100,-5,5,gs->GetEmittedGamma(g)->GetY());
    Fill("posxy_emitted",100,-5,5,gs->GetEmittedGamma(g)->GetX(),100,-5,5,gs->GetEmittedGamma(g)->GetY());
    Fill("posz_emitted",1000,-50,50,gs->GetEmittedGamma(g)->GetZ());
    Fill("phi_emitted",360,-3.14159,3.14159,gs->GetEmittedGamma(g)->GetPhi());
    Fill("theta_emitted",180,0,3.14159,gs->GetEmittedGamma(g)->GetTheta());
    Fill("phi_theta_emitted",360,-3.14159,3.14159,gs->GetEmittedGamma(g)->GetPhi(),180,0,3.14159,gs->GetEmittedGamma(g)->GetTheta());
    emitted_energies.push_back(gs->GetEmittedGamma(g)->GetEnergy());
  }
  Fill("mult_emitted",10,0,10,gs->GetMult());
  
  double betan;
#ifdef USELISA
  Fill("reactiontarget",MAXTARGETS,0,MAXTARGETS,li->GetReaction());
  for(int t=0;t<li->GetNTargets();t++){
    Fill(Form("deltaE_%d",t),fSett->DeltaEBins(),fSett->DeltaERange(0),fSett->DeltaERange(1),li->GetDeltaE(t));
    Fill(Form("deltaE_%d_hit_%d",t,li->GetReaction()),fSett->DeltaEBins(),fSett->DeltaERange(0),fSett->DeltaERange(1),li->GetDeltaE(t));
    if(li->GetReaction()>-1)
      Fill(Form("deltaE_%d_vs_hit_%d",t,li->GetReaction()),fSett->DeltaEBins(),fSett->DeltaERange(0),fSett->DeltaERange(1),li->GetDeltaE(t),MAXTARGETS,0,MAXTARGETS,li->GetReaction());
  }
  if(li->GetReaction()>-1)
    betan = fSett->TargetBeta(li->GetReaction())*(1 + (zd->GetBetaTA() - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
  else betan =0;
#else
  betan = fSett->TargetBeta()*(1 + (zd->GetBetaTA() - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
#endif
  
  double beta_real = gs->GetEmittedGamma(0)->GetBeta();
  // gamma detected in gretina
  if( gr->GetHit(0) && gr->GetHit(0)->GetTS() == gs->GetTimeStamp() ){
    //cout << "coinc " << endl;

    Fill("mult_emitted_detected",100,0,100, gs->GetMult(),100,0,100, gr->GetMult());
    Fill("GRmult_emitted_detected",100,0,100, gs->GetMult(),100,0,100, gr->GetMult());

    Fill("ata",200,-100,100,zd->GetATA());
    Fill("yta",100,-50,50,zd->GetYTA());
    Fill("bta",200,-100,100,zd->GetBTA());
    Fill("ata_vs_bta",200,-100,100,zd->GetATA(),200,-100,100,zd->GetBTA());
    Fill("betata",1000,0,1,zd->GetBetaTA());
    Fill("azita",360,0,2*3.1415926,zd->GetPhi());
    Fill("scatter",500,0,0.5,zd->GetTheta());
    Fill("ptot",500,10,20,zd->GetPtot());
    Fill("ppar",500,10,20,zd->GetPpar());
    Fill("ptra",250,0,5,zd->GetPtra());
    Fill("etot",1000,2000,4000,zd->GetEtot());
    
    Fill("betan",1000,0.4,0.7,betan);

    Fill("betareal",1000,0.4,0.7,beta_real);
    Fill("betadiff",1000,-0.01,0.01, (zd->GetBetaTA() - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
    Fill("betan_betareal",1000,0.4,0.7,betan, 1000,0.4,0.7,beta_real);
    
    for(UShort_t g=0;g<gr->GetMult();g++){ // looping over gamma events
      HitCalc* hit = gr->GetHit(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      if(energy<1 || energy>6000)//useful events
	continue;


      Fill("egam",8000,0,8000,energy);
      Fill("egamdc",8000,0,8000,energy_dc);
      Fill("egamdc_fine",8000,0,4000,energy_dc);
#ifdef USELISA
      if(li->GetReaction()>-1){
	Fill(Form("egamdc_mu%d",li->GetReaction()),8000,0,8000,energy_dc);
	Fill(Form("egamdc_fine_mu%d",li->GetReaction()),8000,0,4000,energy_dc);
	Fill(Form("egamdc_dettheta_mu%d",li->GetReaction()),4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      }
#endif
      Fill("egam_mult",40,0,40,gr->GetMult(),8000,0,8000,energy);
      Fill("egamdc_mult",40,0,40,gr->GetMult(),8000,0,8000,energy_dc);
      Fill("egam_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("egam_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("egam_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("egam_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egam_summary",60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy);
      Fill("egam_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy);
      Fill("egamdc_summary",60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy_dc);
      Fill("egamdc_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy_dc);

      Fill("GRegam",8000,0,8000,energy);
      Fill("GRegamdc",8000,0,8000,energy_dc);
      Fill("GRegam_mult",40,0,40,gr->GetMult(),8000,0,8000,energy);
      Fill("GRegamdc_mult",40,0,40,gr->GetMult(),8000,0,8000,energy_dc);
      Fill("GRegam_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("GRegam_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("GRegam_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("GRegam_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("GRegamdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);

      //positions
      Fill("posxy",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Y());
      Fill("posxz",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Z());
      Fill("posyz",600,-300,300,hit->GetPosition().Y(),600,-300,300,hit->GetPosition().Z());
      Fill("pos_phi",360,-3.14159,3.14159,hit->GetPosition().Phi());
      Fill("pos_theta",180,0,3.14159,hit->GetPosition().Theta());
      Fill("pos_phi_theta",360,-3.14159,3.14159,hit->GetPosition().Phi(),180,0,3.14159,hit->GetPosition().Theta());
      TVector3 hitd = hit->GetPosition();
      hitd.SetX(hitd.X() - zd->GetXTA() - fSett->TargetX());
      hitd.SetY(hitd.Y() - zd->GetYTA() - fSett->TargetY());
#ifdef USELISA
      hitd.SetZ(hitd.Z() - fSett->TargetZ(0));
#else
      hitd.SetZ(hitd.Z() - fSett->TargetZ());
#endif
      Fill("posshift_phi_theta",360,-3.14159,3.14159,hitd.Phi(),180,0,3.14159,hitd.Theta());
      
      Fill("GRposxy",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Y());
      Fill("GRposxz",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Z());
      Fill("GRposyz",600,-300,300,hit->GetPosition().Y(),600,-300,300,hit->GetPosition().Z());
      Fill("GRpos_phi",360,-3.14159,3.14159,hit->GetPosition().Phi());
      Fill("GRpos_theta",180,0,3.14159,hit->GetPosition().Theta());
      Fill("GRpos_phi_theta",360,-3.14159,3.14159,hit->GetPosition().Phi(),180,0,3.14159,hit->GetPosition().Theta());
      Fill("GRposshift_phi_theta",360,-3.14159,3.14159,hitd.Phi(),180,0,3.14159,hitd.Theta());
    }//gamma events



    for(UShort_t g=0;g<gr->GetMultAB();g++){ // looping over gamma events
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)//useful events
	continue;

      HitCalc* hit = gr->GetHitAB(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      Fill("egamAB",8000,0,8000,energy);
      Fill("egamABdc",8000,0,8000,energy_dc);
      Fill("egamABdc_fine",8000,0,4000,energy_dc);
#ifdef USELISA
      if(li->GetReaction()>-1){
	Fill(Form("egamABdc_mu%d",li->GetReaction()),8000,0,8000,energy_dc);
	Fill(Form("egamABdc_fine_mu%d",li->GetReaction()),8000,0,4000,energy_dc);
	Fill(Form("egamABdc_dettheta_mu%d",li->GetReaction()),4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      }
#endif
      Fill("egamAB_mult",40,0,40,gr->GetMult(),8000,0,8000,energy);
      Fill("egamABdc_mult",40,0,40,gr->GetMult(),8000,0,8000,energy_dc);
      Fill("egamAB_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("egamAB_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("egamAB_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("egamAB_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamABdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamAB_summary",60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy);
      Fill("egamAB_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy);
      Fill("egamABdc_summary",60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy_dc);
      Fill("egamABdc_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy_dc);
     
      Fill("GRegamAB",8000,0,8000,energy);
      Fill("GRegamABdc",8000,0,8000,energy_dc);
      Fill("GRegamAB_mult",40,0,40,gr->GetMult(),8000,0,8000,energy);
      Fill("GRegamABdc_mult",40,0,40,gr->GetMult(),8000,0,8000,energy_dc);
      Fill("GRegamAB_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("GRegamAB_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("GRegamAB_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("GRegamAB_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("GRegamABdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);

    }//gamma events


    Fill("mult",30,-0.5,29.5,gr->GetMult());
    Fill("mult_ab",30,-0.5,29.5,gr->GetMultAB());

    Fill("GRmult",30,-0.5,29.5,gr->GetMult());
    Fill("GRmult_ab",30,-0.5,29.5,gr->GetMultAB());

    //coincidences GR with GR
    double enhigh = 0;
    double enlow = 0;
    for(int i=0; i<gr->GetMult(); i++){
      for(int j=i+1; j<gr->GetMult(); j++){
	enhigh = max(gr->GetHit(i)->GetEnergy(),gr->GetHit(j)->GetEnergy());
	enlow = min(gr->GetHit(i)->GetEnergy(),gr->GetHit(j)->GetEnergy());
	if(enlow<1 || enhigh>6000)//useful events
	  continue;
	Fill("egamegam",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("GRegamGRegam",1000,0,4000,enhigh,1000,0,4000,enlow);
	enhigh = max(gr->GetHit(i)->GetDCEnergy(),gr->GetHit(j)->GetDCEnergy());
	enlow = min(gr->GetHit(i)->GetDCEnergy(),gr->GetHit(j)->GetDCEnergy());
	Fill("egamegamdc",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("GRegamGRegamdc",1000,0,4000,enhigh,1000,0,4000,enlow);
      }
    }

    for(int i=0; i<gr->GetMultAB(); i++){
      for(int j=i+1; j<gr->GetMultAB(); j++){
	enhigh = max(gr->GetHitAB(i)->GetEnergy(),gr->GetHitAB(j)->GetEnergy());
	enlow = min(gr->GetHitAB(i)->GetEnergy(),gr->GetHitAB(j)->GetEnergy());
	if(enlow<1 || enhigh>6000)//useful events
	  continue;
	Fill("egamegamAB",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("GRegamGRegamAB",1000,0,4000,enhigh,1000,0,4000,enlow);
	enhigh = max(gr->GetHitAB(i)->GetDCEnergy(),gr->GetHitAB(j)->GetDCEnergy());
	enlow = min(gr->GetHitAB(i)->GetDCEnergy(),gr->GetHitAB(j)->GetDCEnergy());
	Fill("egamegamABdc",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("GRegamGRegamABdc",1000,0,4000,enhigh,1000,0,4000,enlow);
      }
    }
  }//gretina detected

  // gamma detected in miniball
  if( mb->GetHit(0) && mb->GetHit(0)->GetTS() == gs->GetTimeStamp() ){
    Fill("mult_emitted_detected",100,0,100, gs->GetMult(),100,0,100, mb->GetMult());
    Fill("MBmult_emitted_detected",100,0,100, gs->GetMult(),100,0,100, mb->GetMult());

    Fill("ata",200,-100,100,zd->GetATA());
    Fill("yta",100,-50,50,zd->GetYTA());
    Fill("bta",200,-100,100,zd->GetBTA());
    Fill("ata_vs_bta",200,-100,100,zd->GetATA(),200,-100,100,zd->GetBTA());
    Fill("betata",1000,0,1,zd->GetBetaTA());
    Fill("azita",360,0,2*3.1415926,zd->GetPhi());
    Fill("scatter",500,0,0.5,zd->GetTheta());
    Fill("ptot",500,10,20,zd->GetPtot());
    Fill("ppar",500,10,20,zd->GetPpar());
    Fill("ptra",250,0,5,zd->GetPtra());
    Fill("etot",1000,2000,4000,zd->GetEtot());


    Fill("betan",1000,0.4,0.7,betan);
    Fill("betareal",1000,0.4,0.7,beta_real);
    Fill("betadiff",1000,-0.01,0.01, (zd->GetBetaTA() - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
    Fill("betan_betareal",1000,0.4,0.7,betan, 1000,0.4,0.7,beta_real);


    for(UShort_t g=0;g<mb->GetMult();g++){ // looping over gamma events
      MBHitCalc* hit = mb->GetHit(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      if(energy<1 || energy>6000)//useful events
	continue;


      Fill("egam",8000,0,8000,energy);
      Fill("egamdc",8000,0,8000,energy_dc);
      Fill("egam_mult",40,0,40,mb->GetMult(),8000,0,8000,energy);
      Fill("egamdc_mult",40,0,40,mb->GetMult(),8000,0,8000,energy_dc);
      Fill("egam_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("egam_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("egam_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("egam_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egam_summary", 60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy);
      Fill("egam_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy);
      Fill("egamdc_summary", 60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy_dc);
      Fill("egamdc_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy_dc);
      
      Fill("MBegam",8000,0,8000,energy);
      Fill("MBegamdc",8000,0,8000,energy_dc);
      Fill("MBegam_mult",40,0,40,mb->GetMult(),8000,0,8000,energy);
      Fill("MBegamdc_mult",40,0,40,mb->GetMult(),8000,0,8000,energy_dc);
      Fill("MBegam_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("MBegam_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("MBegam_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("MBegam_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("MBegamdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);

      //positions
      Fill("posxy",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Y());
      Fill("posxz",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Z());
      Fill("posyz",600,-300,300,hit->GetPosition().Y(),600,-300,300,hit->GetPosition().Z());
      Fill("pos_phi",360,-3.14159,3.14159,hit->GetPosition().Phi());
      Fill("pos_theta",180,0,3.14159,hit->GetPosition().Theta());
      Fill("pos_phi_theta",360,-3.14159,3.14159,hit->GetPosition().Phi(),180,0,3.14159,hit->GetPosition().Theta());
      TVector3 hitd = hit->GetPosition();
      hitd.SetX(hitd.X() - zd->GetXTA() - fSett->TargetX());
      hitd.SetY(hitd.Y() - zd->GetYTA() - fSett->TargetY());
#ifdef USELISA
      hitd.SetZ(hitd.Z() - fSett->TargetZ(0));
#else
      hitd.SetZ(hitd.Z() - fSett->TargetZ());
#endif
      Fill("posshift_phi_theta",360,-3.14159,3.14159,hitd.Phi(),180,0,3.14159,hitd.Theta());

      Fill("MBposxy",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Y());
      Fill("MBposxz",600,-300,300,hit->GetPosition().X(),600,-300,300,hit->GetPosition().Z());
      Fill("MBposyz",600,-300,300,hit->GetPosition().Y(),600,-300,300,hit->GetPosition().Z());
      Fill("MBpos_phi",360,-3.14159,3.14159,hit->GetPosition().Phi());
      Fill("MBpos_theta",180,0,3.14159,hit->GetPosition().Theta());
      Fill("MBpos_phi_theta",360,-3.14159,3.14159,hit->GetPosition().Phi(),180,0,3.14159,hit->GetPosition().Theta());
      Fill("MBposshift_phi_theta",360,-3.14159,3.14159,hitd.Phi(),180,0,3.14159,hitd.Theta());
    }//gamma events

    for(UShort_t g=0;g<mb->GetMultAB();g++){ // looping over gamma events
      if(mb->GetHitAB(g)->GetEnergy()<1 || mb->GetHitAB(g)->GetEnergy()>6000)//useful events
	continue;

      MBHitCalc* hit = mb->GetHitAB(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      Fill("egamAB",8000,0,8000,energy);
      Fill("egamABdc",8000,0,8000,energy_dc);
      Fill("egamAB_mult",40,0,40,mb->GetMult(),8000,0,8000,energy);
      Fill("egamABdc_mult",40,0,40,mb->GetMult(),8000,0,8000,energy_dc);
      Fill("egamAB_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("egamAB_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("egamAB_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("egamAB_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamABdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("egamAB_summary", 60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy);
      Fill("egamAB_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy);
      Fill("egamABdc_summary", 60,0,60,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),4000,0,4000,energy_dc);
      Fill("egamABdc_detector",20,0,20,fSett->Clu2Det(hit->GetCluster()),4000,0,4000,energy_dc);

      Fill("MBegamAB",8000,0,8000,energy);
      Fill("MBegamABdc",8000,0,8000,energy_dc);
      Fill("MBegamAB_mult",40,0,40,mb->GetMult(),8000,0,8000,energy);
      Fill("MBegamABdc_mult",40,0,40,mb->GetMult(),8000,0,8000,energy_dc);
      Fill("MBegamAB_highestemitted_detected",1000,0,4000,highestemitted,1000,0,4000,energy);     
      for(unsigned short i=0;i<emitted_energies.size();i++){
	Fill("MBegamAB_emitted_diff",5000,-1000,4000,emitted_energies.at(i)-energy);
	Fill("MBegamAB_emitted_detected",1000,0,4000,emitted_energies.at(i),1000,0,4000,energy);
      }
      Fill("MBegamAB_dettheta",4000,0,4000,energy,360,0,180,hit->GetPosition().Theta()*180/3.14159);
      Fill("MBegamABdc_dettheta",4000,0,4000,energy_dc,360,0,180,hit->GetPosition().Theta()*180/3.14159);
    }//gamma events

    Fill("mult",30,-0.5,29.5,mb->GetMult());
    Fill("mult_ab",30,-0.5,29.5,mb->GetMultAB());

    Fill("MBmult",30,-0.5,29.5,mb->GetMult());
    Fill("MBmult_ab",30,-0.5,29.5,mb->GetMultAB());
    
    //coincidences GR with GR
    double enhigh = 0;
    double enlow = 0;
    for(int i=0; i<mb->GetMult(); i++){
      for(int j=i+1; j<mb->GetMult(); j++){
	enhigh = max(mb->GetHit(i)->GetEnergy(),mb->GetHit(j)->GetEnergy());
	enlow = min(mb->GetHit(i)->GetEnergy(),mb->GetHit(j)->GetEnergy());
	if(enlow<1 || enhigh>6000)//useful events
	  continue;
	Fill("egamegam",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("MBegamMBegam",1000,0,4000,enhigh,1000,0,4000,enlow);
	enhigh = max(mb->GetHit(i)->GetDCEnergy(),mb->GetHit(j)->GetDCEnergy());
	enlow = min(mb->GetHit(i)->GetDCEnergy(),mb->GetHit(j)->GetDCEnergy());
	Fill("egamegamdc",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("MBegamMBegamdc",1000,0,4000,enhigh,1000,0,4000,enlow);
      }
    }

    for(int i=0; i<mb->GetMultAB(); i++){
      for(int j=i+1; j<mb->GetMultAB(); j++){
	enhigh = max(mb->GetHitAB(i)->GetEnergy(),mb->GetHitAB(j)->GetEnergy());
	enlow = min(mb->GetHitAB(i)->GetEnergy(),mb->GetHitAB(j)->GetEnergy());
	if(enlow<1 || enhigh>6000)//useful events
	  continue;
	Fill("egamegamAB",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("MBegamMBegamAB",1000,0,4000,enhigh,1000,0,4000,enlow);
	enhigh = max(mb->GetHitAB(i)->GetDCEnergy(),mb->GetHitAB(j)->GetDCEnergy());
	enlow = min(mb->GetHitAB(i)->GetDCEnergy(),mb->GetHitAB(j)->GetDCEnergy());
	Fill("egamegamABdc",1000,0,4000,enhigh,1000,0,4000,enlow);
	Fill("MBegamMBegamABdc",1000,0,4000,enhigh,1000,0,4000,enlow);
      }
    }

    //MB and GR coincidences
    double enmb = 0;
    double engr = 0;
    if( gr->GetHit(0) && gr->GetHit(0)->GetTS() == gs->GetTimeStamp() ){
      for(int i=0; i<mb->GetMult(); i++){
	enmb = mb->GetHit(i)->GetEnergy();
	if(enmb<1 || enmb>6000)//useful events
	  continue;
	for(int j=0; j<gr->GetMult(); j++){
	  enmb = mb->GetHit(i)->GetEnergy();
	  engr = gr->GetHit(j)->GetEnergy();
	  if(engr>enmb){
	    Fill("egamegam",1000,0,4000,engr,1000,0,4000,enmb);
	    Fill("GRegamMBegam",1000,0,4000,engr,1000,0,4000,enmb);
	  }
	  else if(engr>0){
	    Fill("egamegam",1000,0,4000,enmb,1000,0,4000,engr);
	    Fill("MBegamGRegam",1000,0,4000,enmb,1000,0,4000,engr);
	  }
	  enmb = mb->GetHit(i)->GetDCEnergy();
	  engr = gr->GetHit(j)->GetDCEnergy();
	  if(engr>enmb && enmb>0){
	    Fill("egamegamdc",1000,0,4000,engr,1000,0,4000,enmb);
	    Fill("GRegamMBegamdc",1000,0,4000,engr,1000,0,4000,enmb);
	  }
	  else if(engr>0){
	    Fill("egamegamdc",1000,0,4000,enmb,1000,0,4000,engr);
	    Fill("MBegamGRegamdc",1000,0,4000,enmb,1000,0,4000,engr);
	  }
	}
      }// multiplicity
      for(int i=0; i<mb->GetMultAB(); i++){
	enmb = mb->GetHitAB(i)->GetEnergy();
	if(enmb<1 || enmb>6000)//useful events
	  continue;
	for(int j=0; j<gr->GetMultAB(); j++){
	  enmb = mb->GetHitAB(i)->GetEnergy();
	  engr = gr->GetHitAB(j)->GetEnergy();
	  if(engr>enmb){
	    Fill("egamegamAB",1000,0,4000,engr,1000,0,4000,enmb);
	    Fill("GRegamMBegamAB",1000,0,4000,engr,1000,0,4000,enmb);
	  }
	  else if(engr>0){
	    Fill("egamegamAB",1000,0,4000,enmb,1000,0,4000,engr);
	    Fill("MBegamGRegamAB",1000,0,4000,enmb,1000,0,4000,engr);
	  }
	  enmb = mb->GetHitAB(i)->GetDCEnergy();
	  engr = gr->GetHitAB(j)->GetDCEnergy();
	  if(engr>enmb && enmb>0){
	    Fill("egamegamABdc",1000,0,4000,engr,1000,0,4000,enmb);
	    Fill("GRegamMBegamABdc",1000,0,4000,engr,1000,0,4000,enmb);
	  }
	  else if(engr>0){
	    Fill("egamegamABdc",1000,0,4000,enmb,1000,0,4000,engr);
	    Fill("MBegamGRegamABdc",1000,0,4000,enmb,1000,0,4000,engr);
	  }
	}
      }// multiplicity
    }//tracking detected
  }//miniball detected

}
