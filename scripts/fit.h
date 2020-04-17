Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t f2gammagaussbg(Double_t *x, Double_t *par);
Double_t f2gammabg(Double_t *x, Double_t *par);
Double_t f2gammastep(Double_t *x, Double_t *par);
Double_t f2gammabgstep(Double_t *x, Double_t *par);
Double_t f2gammagaus0(Double_t *x, Double_t *par);
Double_t f2gammagaus1(Double_t *x, Double_t *par);

Double_t flinear(Double_t *x, Double_t *par);
using namespace std;

// one peak fitting function
Double_t fgammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t fgammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  //  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];
  arg = (x[0]-mean)/(sqrt2*sigma);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];


  arg = (x[0]-mean)/(sqrt2*sigma);
  Double_t result = 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  return result;

}
// two peak fitting function
Double_t f2gammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];

  Double_t step = par[8];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t f2gammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];

  Double_t step = par[8];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammabgstep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0]; 
  
  
  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];

  Double_t step = par[8];

  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += step/pow(1+exp(sqrt2*arg),2);
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t f2gammagaus0(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  double result = 0;
  arg = (x[0]-mean1)/(sqrt2*sigma1);
  result += 1/(sqrt2pi*sigma1) * norm1 * exp(-arg*arg);

  return result;

}
Double_t f2gammagaus1(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm1  = par[2];
  Double_t mean1  = par[3];
  Double_t sigma1 = par[4];
  Double_t norm2  = par[5];
  Double_t mean2  = par[6];
  Double_t sigma2 = par[7];
  double result = 0;
  arg = (x[0]-mean2)/(sqrt2*sigma2);
  result += 1/(sqrt2pi*sigma2) * norm2 * exp(-arg*arg);

  return result;

}
Double_t flinear(Double_t *x, Double_t *par){
  return x[0]*par[0] + par[1];

}
Double_t freso(Double_t *x, Double_t *par){
  return par[0]*sqrt(1+par[1]*x[0]) + par[2]*x[0];

}
Double_t fthresh(Double_t *x, Double_t *par){
  return par[2]*(1.0 + tanh( (x[0]-par[0])/par[1] ));
}
Double_t fthreshbg(Double_t *x, Double_t *par){
  return par[2]*(1.0 + tanh( (x[0]-par[0])/par[1] )) + par[3];
}
