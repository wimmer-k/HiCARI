using namespace std;

void rootlogonKW() {
  gSystem->Load("libHiCARI");
  gROOT->ProcessLine(".L ~/HiCARI/scripts/FitWindow/fitwindow.C+");
  cout<<"loaded mouseclick fitting / "<<flush;
  gROOT->ProcessLine("window()");


  gROOT->ForceStyle(1);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  gStyle->SetHistLineColor(1);
  gStyle->SetOptStat(1111111);
  TColor *orange = gROOT->GetColor(5);
  orange->SetRGB(1.0, 0.612, 0.002); 
  
  TColor *green = gROOT->GetColor(3);
  green->SetRGB(0.15, 0.7, 0.15);
  
  cout<<"set style"<<endl;
  system("pwd");

}
