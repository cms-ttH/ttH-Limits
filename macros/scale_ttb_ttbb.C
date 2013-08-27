//Root includes                                   
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TAxis.h"
#include "TKey.h"
#include "TList.h"

void scale_ttb_ttbb(TString inputFileName = "input.root", TString outputFileName = "output.root") {


  TFile *outFile = TFile::Open(outputFileName,"RECREATE");

  TFile *histFile = TFile::Open(inputFileName);

  std::cout << "Copying over all histograms..." << std::endl;

  //Now copy all other histograms in the file as is (preserving names...)
  TList *keys = histFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();

    double scale = 1.0;
    if( name.Contains("ttbarPlusBBbar") )  scale = 1.66646;
    else if( name.Contains("ttbarPlusB") ) scale = 1.26073;

    TH1 *hist = 0;
    hist = (TH1 *)histFile->Get(name);

    hist->Scale(scale);

    outFile->cd();
    hist->Write(name);
  }

  outFile->Close();

  std::cout << "Done." << std::endl;

}
