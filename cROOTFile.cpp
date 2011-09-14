
#include <iostream>

#include "TApplication.h"
#include "TObject.h"
#include "TROOT.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCut.h"
#include "TSpectrum.h"
#include "TStyle.h"

extern "C" { 
#include "monarch.h"
#include "cFFTransform.h"
}
#include "cROOTFile.h"

using namespace std;

cROOTFile::cROOTFile(char *filename) {
  //constructor goes here
  TString rootName = TString(filename);
  if (rootName.EndsWith(".egg")) {
    rootName.Remove(rootName.Length()-4,4);
  }
  rootName+=".root";
  fFile = new TFile(rootName, "recreate");
}

cROOTFile::~cROOTFile() {
  //destructor goes here
  if(fFile) delete fFile;
}

void cROOTFile::writeTimeTrees(struct FFTinput *input, int event) {
  if ( fFile==NULL) {
    cout << "Error! Root File does not exist;" << endl; 
  }  else if ( fFile->IsOpen() ) {
    //time series data
    TTree *tsTree = new TTree(Form("tsTree_%d",event), "time series tree");
    float t_val;
    double volt;
    double pow_t;
    tsTree->Branch("t", &t_val);
    tsTree->Branch("volt", &volt);
    tsTree->Branch("pow_t", &pow_t);
    int i;
    for(i = 0; i < input->in_size; i++) {
      
      t_val = ( (float)i * (input->step_size));
      
      volt = input->in[i];
      pow_t= volt * volt;
      tsTree->Fill();
    }
    tsTree->Write("", TObject::kOverwrite);
    tsTree->Delete(); 
    
  } else {
    cout << "Root File not opened yet!" << endl; 
  }
}
void cROOTFile::fitTimeNoise(int event = 0)
{
  
  //cout << "Fitting Time Noise " << endl;
  TTree *tsTree = (TTree*)fFile->Get(Form("tsTree_%d",event));
  tsTree->Draw("volt>>htemp","","goff");
  TH1F *h1 = (TH1F*)gDirectory->Get("htemp");
  h1->Fit("gaus", "q0","");
  h1->GetFunction("gaus")->SetLineColor(2);
  h1->SetName(Form("hNT_%d",event));
  h1->SetTitle("Noise in Time Domain");
  h1->SetXTitle("voltage (arb units)");
  h1->Write("", TObject::kOverwrite);
  delete h1; 
}

void cROOTFile::writeFreqTrees(struct FFTinput *input, int event) 
{
  if ( fFile==NULL) {
    cout << "Error! Root File does not exist;" << endl; 
  }  else if ( fFile->IsOpen() ) {
    //freq. space data
    TTree *fsTree = new TTree(Form("fsTree_%d",event), "frequency space tree");
    Double_t f_val;
    double a;
    double b;
    double aSquared;
    double bSquared;
    Double_t pow_f;
    
    fsTree->Branch("f", &f_val);
    fsTree->Branch("pow_f", &pow_f);
    
    int i;
    for(i = 0; i < input->out_size; i++) {
      
      f_val = ( (float)i / input->in_size) * (1 / (input->step_size));
      a = input->out[i][0];
      b = input->out[i][1];
      
      aSquared = a * a;
      bSquared = b * b;
      
      
      pow_f = aSquared + bSquared;
      fsTree->Fill();
    }
    fsTree->Write("", TObject::kOverwrite);
    fsTree->Delete();
  } else {
    cout << "Root File not opened yet!" << endl; 
  }
  
}
void cROOTFile::createPS(int event)
{
  //create ps histogram
  TTree *fsTree = (TTree*)fFile->Get(Form("fsTree_%d",event));
  fsTree->SetEstimate(fsTree->GetEntries());
  fsTree->Draw("pow_f:f","","goff");
  Int_t nBins = fsTree->GetSelectedRows();
  Double_t *x = fsTree->GetV2();
  Double_t *y = fsTree->GetV1();
  TH1F *hPS = new TH1F(Form("hPS_%d",event), "Power Spectrum",nBins-1, x); 
  hPS->SetContent(y);
  hPS->Write("", TObject::kOverwrite);
  hPS->Delete();
  fsTree->Delete();


}
void cROOTFile::findPeaks(int event, double* ia, int &nPeaks, const int maxP = 10 )
{
  //get ps histogram
  TH1F *hPS = (TH1F*)fFile->Get(Form("hPS_%d",event));
  //search for peaks
  TSpectrum *s = new TSpectrum(maxP, 1);
  Int_t nfound = s->Search(hPS, 1e-7, "", 5e-5);
  Float_t *xpeaks = s->GetPositionX();
  
  //set found values
  nPeaks = nfound;
  for (int i = 0; i < nfound; i++) {
    ia[i] = (double)xpeaks[i];
  }
  hPS->Write("", TObject::kOverwrite);
  hPS->Delete();

}
void cROOTFile::fitFreqNoise(int event = 0, double *ia = NULL, int nPeaks = 0)
{
  //cout << "Fitting Freq Noise " << endl;
  //cut out identified peaks
  double peakWidth = 0.02; //MHz
  TCut freqCut = Form("f>%f",peakWidth);//always cut out DC
  for (int i = 0; i< nPeaks; i++ ) {
    freqCut+= Form("(f<%f||f>%f)",ia[i]-peakWidth, ia[i]+peakWidth);
    //cout << "Cutting out " << ia[i] << endl; 

  }
  TTree *fsTree = (TTree*)fFile->Get(Form("fsTree_%d",event));
  fsTree->SetEstimate(fsTree->GetEntries());
  fsTree->Draw("pow_f>>htemp(1000)",freqCut,"goff");
  TH1F *h1 = (TH1F*)gDirectory->Get("htemp");
  
  //set up exponential fit function 
  TF1 exp("exp", "[1]/[0]*TMath::Exp(-x/[0])",h1->GetMinimum(), h1->GetMaximum());
  exp.SetParameter(0, h1->GetMean());
  exp.SetParName(0, "kT");
  exp.SetParameter(1, h1->GetEntries());
  exp.SetParName(1, "Norm");
  exp.SetLineColor(2);
  h1->Fit("exp", "q0","");
  h1->SetName(Form("hNF_%d",event));
  h1->SetTitle("Noise in Freq Domain");
  h1->SetXTitle("power (arb units)");

  //plot PS with peaks cut out 
  freqCut+= "pow_f>0.4e9&&pow_f<1.6e9";//range above noise
  fsTree->Draw("pow_f:f>>graph",freqCut,"goff");
  TGraph *graph = new TGraph(fsTree->GetSelectedRows(),fsTree->GetV2(), fsTree->GetV1());    
  graph->SetName(Form("gPS_Cut_%d", event));
  graph->SetTitle("PS after cuts");
  graph->GetHistogram()->SetXTitle("Freq (MHz)");

  h1->Write("", TObject::kOverwrite);
  graph->Write("", TObject::kOverwrite);
  delete h1; 
  fsTree->Delete(); 
}

void cROOTFile::fitFreqPeak(int event = 0, double *ia = NULL, int nPeaks = 0)
{
  //cout << "Fitting Peaks in Freq" << endl;
  double peakWidth = 0.2; //MHz
  TTree *fsTree = (TTree*)fFile->Get(Form("fsTree_%d",event));
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1);
  fsTree->GetEntry(1);
  double delf = fsTree->GetLeaf("f")->GetValue(); //MHz
  
  //build graphs around identified peaks
  for (Int_t i = 0; i< nPeaks; i++ ) {
    double f_low = ia[i]-peakWidth;
    double f0 = ia[i];
    double f_high = ia[i]+peakWidth;
    fsTree->Draw("pow_f:f",Form("(f>%f&&f<%f)",f_low, f_high), "goff");
    TGraph *graph = new TGraph(fsTree->GetSelectedRows(),fsTree->GetV2(), fsTree->GetV1());    
    graph->SetName(Form("sincFit%d_%d", event, i));
    graph->SetTitle(Form("Peak %d in Freq Domain", i));
    graph->GetHistogram()->SetXTitle("Freq (MHz)");
    
    //appropriate for rectangular windowing
    TF1 *sinc = new TF1("sinc", "[2]*TMath::Sin(2*TMath::Pi()*(x-[0])*[1]/2)^2/(2*TMath::Pi()*(x-[0])*[1]/2)^2", f_low, f_high);
    sinc->SetParameters(f0, 1 / delf, graph->GetHistogram()->GetMaximum());
    //cout << "Initializing " << f0 << " " <<  1 / delf << " ";
    //cout <<  graph->GetHistogram()->GetMaximum() << endl;
    sinc->SetParLimits(0, f_low, f_high);
    sinc->SetParLimits(1, 0, 1000 / delf);
    sinc->SetParNames("mean", "duration", "amplitude");
    sinc->SetLineColor(2);
    sinc->SetNpx(1000);
    int status = graph->Fit(sinc,"q","",f0-0.002, f0+0.002);
    if (!status == 0) {
      status = graph->Fit(sinc, "Mq","",f0-0.002, f0+0.002);
    }
    graph->Write("", TObject::kOverwrite);
    delete sinc; 
    delete graph; 
  }
  fsTree->Delete(); 
}

void cROOTFile::closeFile() {
  if ( !(fFile==NULL) ) {
    if (fFile->IsOpen() ) {
      fFile->Close();
    }
  }
}
