char *cutfname =NULL;
bool filecalled = false;
void setFile(char* fname){
  cout << "setting filename to " << fname << endl;
  cutfname = fname;
  filecalled = true;
}
void makeCuts(){
  cout << "================================================================================" << endl;
  cout << "Draw and save 2D cuts in ROOT" << endl;
  cout << "     K. Wimmer (NSCL) " << endl;
  cout << "--------------------------------------------------------------------------------" << endl;  
  cout << "- set output filename by calling \"setFile(\"YOURFILE.root\")\" "<<endl;
  cout << "  default name is \"cutfile.root\" " <<endl;
  cout << "- start cut by calling \"cut(\"YOURNAME\")\" " << endl;
  cout << "================================================================================" << endl;
 
}
void cut(char* name){  
  TCutG *cut;
  cout << "You can draw a cut by clicking View->Toolbar->Graphical Cut" << endl;
  cout << "Double-click to close cut" << endl;
  cut = (TCutG*)gPad->WaitPrimitive("CUTG");
  cut->SetName(name); 
  cout << "Cut " << name << " with "<< cut->GetN() << " points "<< endl;
  vector<double> x;
  vector<double> y;
  x.resize(cut->GetN());
  y.resize(cut->GetN());
  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x[n],y[n]);
    cout << x[n] << "\t" << y[n] << endl;
  }
  TFile cutf*;
  if(filecalled==false||cutfname ==NULL)
    cutf = new TFile("cutfile.root","update");
  else
    cutf = new TFile(cutfname,"update");

  cout << "Cut " << name << " will be saved to file " << cutf->GetName() << endl;
  cutf->cd();
  cut->Write(name,TObject::kOverwrite);
  cutf->ls();
  cutf->Close();
}
