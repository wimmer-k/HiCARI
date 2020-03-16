#include "Tracking.hh"
#include "perm.h"
#include "crosssections.h"
#include <functional>

void Tracking::SortInClusters(){
  fesum.resize(fgr->GetNCluster());
  fevt = new GretinaEvent();
  for(int j=0;j<fgr->GetNCluster();j++){
    if(fset->VLevel()>2)
      cout << "---------CLUSTER"<<j<<"------------------------" << endl;
    //fgr->PrintCluster(j);
    //sort interactions
    vector<HitCalc*> hits = SortInCluster(fgr->GetClusters(j));
    //find sum of interactions
    fesum[j]=0;
    for(vector<HitCalc*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
      fesum[j] += (*hit)->GetEnergy();
    }
    //cout << " ---- sum ---- " << fesum[j] << endl;
    GretinaTrack *best = new GretinaTrack();
    best->SetEsum(fesum[j]);
    fcurrentsize = hits.size();
    if(fcurrentsize>MAXNIPOINTS){
      //cout << " just one hit " << endl;
      best->SetHits(hits);
      best->SetFOM(fcurrentsize-1);
    }
    else{
      double max =0;
      int bestperm =0;
      int badctr =0;
      int jumpctr =0;
      for(UShort_t permutation =0; permutation<fnperm[fcurrentsize];permutation++){
	double fom = FigureOfMerit(hits,permutation,fesum[j]);
	if(fcurrentsize==1){
	  if(!(fset->Probabilities()&(1<<0)))
	    fom=0;
	  max = fom;
	  break;
	}
	if(fset->VLevel()>2){
	  cout << "permutations to consider: " << fnperm[fcurrentsize] << ", current permutation " << permutation << endl;
	  cout <<"FOM " << fom << endl;
	  for(UShort_t h=0;h<fcurrentsize;h++)
	    cout << hits.at(fperm[fcurrentsize][permutation][h])->GetEnergy() << endl;
	  cout << "---------------------------------" << endl;
	}
	
	//count how many bad FOM in a row
	if(fom<fset->JumpFOM()){
	  badctr++;
	  if(fset->VLevel()>2){
	    cout << "bad permutation number" << permutation << "; " << badctr << " bad one in a row " << endl;
	  }
	}	
	else
	  badctr =0;
	
	if(fom>max){
	  max = fom;
	  bestperm = permutation;
	}
	//if there were enough bad ones in a row
	if(badctr>fset->MaxBad() && badctr<fcurrentsize){
	  int ctr = permutation%fnperm[fcurrentsize-1];
	  if(fset->VLevel()>2){
	    cout << "too many bad FOM in a row at permutation " << permutation << " out of " << fnperm[fcurrentsize]<< endl;
	    cout << ctr << " steps done with this first interaction point, which would give a total of " << fnperm[fcurrentsize-1]<< endl;
	  }
	  permutation -= ctr;
	  permutation += fnperm[fcurrentsize-1];
	  if(fset->VLevel()>2){
	    cout << "jumping to " << permutation << endl;
	  }
	  jumpctr++;
	  badctr =0;
	  //return;
	}

      }
      if(fset->VLevel()>2){
	cout << "best was "<< max << " at permutation " << bestperm << " out of  " << fnperm[fcurrentsize] <<" permutations (currentsize = " << fcurrentsize << ")" << endl; 
	for(UShort_t h=0;h<fcurrentsize;h++)
	  cout << hits.at(fperm[fcurrentsize][bestperm][h])->GetEnergy() << endl;
      }
      // bool print = false;
      // for(UShort_t h=0;h<fcurrentsize;h++){
      // 	if(hits.at(h)->GetEnergy()>954.9 && hits.at(h)->GetEnergy()<955.0)
      // 	  print = true;	
      // }
      // if(print){
      // 	cout << "best was "<< max << " at permutation " << bestperm << " out of  " << fnperm[fcurrentsize] <<" permutations (currentsize = " << fcurrentsize << ")" << endl; 
      // 	for(UShort_t h=0;h<fcurrentsize;h++)
      // 	  cout << hits.at(fperm[fcurrentsize][bestperm][h])->GetEnergy() << endl;
      // }
      vector<HitCalc*> besthits;
      besthits.resize(fcurrentsize);
      //update order
      for(int h=0;h<fcurrentsize;h++){
	if(fcurrentsize<MAXNIPOINTS+1)
	  besthits.at(h) = hits.at(fperm[fcurrentsize][bestperm][h]);
      }
      best->SetHits(besthits);
      best->SetFOM(max);
      best->SetJumps(jumpctr);
      best->SetPermutation(bestperm);
    }
    fevt->AddTrack(best);
    //fevt->PrintEvent();
    
  }//clusters
}

vector<HitCalc*> Tracking::SortInCluster(vector<HitCalc*> hits){
  sort(hits.begin(), hits.end(), HitComparer());
  // for(vector<HitCalc*>::iterator hit0=hits.begin(); hit0!=hits.end(); hit0++){
  //   cout << (*hit0)->GetEnergy() << endl;
  // }
  return hits;
}
double Tracking::FigureOfMerit(vector<HitCalc*> hits, int nperm, double esum){
  int iter = 0;
  double fom = 1;
  TVector3 in = hits.at(fperm[fcurrentsize][nperm][0])->GetPosition();
  
  if(fset->VLevel()>2)
    cout << " total energy " << esum << endl;



  //Probability of the gamma reaching the first interaction point without scattering.
  //Used if 2nd bit in probabilities is set, or if 1st bit is set and is a mult==1 event.
  if( (fset->Probabilities()&(1<<0) && fcurrentsize==1)|| fset->Probabilities()&(1<<1)){
    //single hit, calc distance in Ge
    //travel from source to hit 0;
    TVector3 p;
    p = in+fset->TargetPos();

    double alpha = p.Angle(in);
    double r = p.Mag()-fset->RadiusGE();
    r /= cos(alpha);
    double sigma = co(esum)+ph(esum)+pa(esum);
    double lambda = 72.64/(6.022e23*5.323*sigma*1e-24);
    lambda *=10; // to mm
    if(fset->Probabilities()&(1<<0))
      fom *=exp(-r/lambda);
    if(fset->VLevel()>2){
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(5) << endl;
      cout << "in.Mag " << in.Mag() << "p.Mag " << p.Mag() << " r " << r << endl;
      cout << "energy:\t" << esum  << "\ttotal:\t" << sigma << endl;
      cout << "distance travelled " << r << "\tlambda " << lambda << "\tprobability " << exp(-r/lambda) << endl;
    }
    if(fcurrentsize==1)
      return -fom;
  }

  
  for(UShort_t h=1;h<fcurrentsize;h++){
    HitCalc* hit = hits.at(fperm[fcurrentsize][nperm][h]);
    if(fset->VLevel()>2)
      cout << "esum " << esum<< " E["<<h<<"] = " << hit->GetEnergy() << " E["<<iter<<"] = " << hits.at(fperm[fcurrentsize][nperm][iter])->GetEnergy() << endl;
    TVector3 scatter = hit->GetPosition() - hits.at(fperm[fcurrentsize][nperm][iter])->GetPosition();
    double angle = scatter.Angle(in);
    if(fset->VLevel()>2)
      cout << "angle: "<< angle <<" = "<< angle *180/3.1415 << " deg "<< endl;
    double es = esum/(1+esum/511.*(1-cos(angle)));
    if(fset->VLevel()>2)
      cout << "escatter = " << es << " measured = " << esum - hits.at(fperm[fcurrentsize][nperm][iter])->GetEnergy() << endl;

    //updating etot
    esum-=hits.at(fperm[fcurrentsize][nperm][iter])->GetEnergy();
    //updating incoming vector now from point [iter] to point [h]
    in = scatter;

    //calculating probability to travel from point [iter] to point [h]
    if(fset->Probabilities()>1){
      double sigma = co(esum)+ph(esum)+pa(esum);
      double distance = in.Mag(); // in mm
      // density: 5.323 g/cm^3
      // atomic mass 72.64 g/mol
      double lambda = 72.64/(6.022e23*5.323*sigma*1e-24);
      lambda *=10; // to mm
      if(fset->VLevel()>2){
	cout << "energy for cross sections " << esum << endl;
	cout << "compton:\t" << co(esum) << endl;
	cout << "photopeak:\t" << ph(esum) << endl;
	cout << "pair crea:\t" << pa(esum) << endl;
	cout << "energy:\t" << esum  << "\ttotal:\t" << sigma << endl;
	cout << "distance travelled " << distance << "\tlambda " << lambda << "\tprobability " << exp(-distance/lambda) << endl;
      }
      if(fset->Probabilities()&(1<<2)){
	if(h<fcurrentsize-1){
	  fom*=co(esum)/sigma;
	  //cout << "prob for compton " << co(esum)/sigma << endl;
	}
	else{
	  fom*=ph(esum)/sigma;
	  //cout << "prob for photo " << ph(esum)/sigma << endl;
	}
      }
      if(fset->Probabilities()&(1<<1)){
	fom*=exp(-distance/lambda);
      }
    }
    
    //cout <<  "fom += exp(-2*(es-esum)*(es-esum)) = exp(-2*"<<(es-esum)<<"*"<<(es-esum)<<")= exp("<<-2*(es-esum)*(es-esum)<<")"<< endl;
    if(iter==0)
      fom *= exp(-2*(es-esum)*(es-esum)/1e6);
    else
      fom *= exp(-(es-esum)*(es-esum)/1e6);
    if(fset->VLevel()>2)
      cout << "FOM " << fom << endl;
    iter++;
  }
  fom = pow(fom,1./(2*iter-1));
  if(fset->VLevel()>2)
    cout << " this hits FOM " << fom << endl;
  return fom;
}
void Tracking::GeneratePermutations(){
  int current[MAXNIPOINTS+1];
  for(int j=0;j<MAXNIPOINTS+1;j++){
    int n=0;
    for(UShort_t i=0;i<j+1;i++){
      current[i] = i;
      fperm[j][n][i] =  current[i];
    }
    //cout << n << "\t"<< j << "\t";
    //printlist(fperm[j][n], j);
    
    while(get_next_perm(current, j)){
      n++;
      for(UShort_t i=0;i<j+1;i++){
	fperm[j][n][i] =  current[i];
      }
      //cout << n << "\t"<< j << "\t";
      //printlist(fperm[j][n], j);
    }
    fnperm[j]=n+1;
    cout << "generated " <<fnperm[j] << " permutations for N = " << j << endl;
  }//possible number of hits
  fcurrentsize = 0;
}
