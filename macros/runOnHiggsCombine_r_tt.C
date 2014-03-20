#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TList.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include <sstream>


void runOnHiggsCombine_r_tt( int maxNentries=-1, TString inputFile="higgsCombine.mH125.root", TString postfix = "", double rMin=0, double rMax=6 ){

  TFile *treefile = new TFile( inputFile );

  TTree *tree = (TTree*)treefile->Get("limit");


  Float_t r_tt, quantileExpected, deltaNLL;

  tree->SetBranchAddress("r_tt", &r_tt );
  tree->SetBranchAddress("quantileExpected", &quantileExpected );
  tree->SetBranchAddress("deltaNLL", &deltaNLL );


  TString histofilename = "HistoFiles/runOnHiggsCombine"+postfix+".root";

  TFile histofile(histofilename,"recreate");
  histofile.cd();


  std::vector<double> vec_r_tt, vec_deltaNLL;


  TH1D* h_quantileExpected = new TH1D("h_quantileExpected",";quantileExpected",10000,-1,1);
  TH1D* h_deltaNLL = new TH1D("h_deltaNLL",";deltaNLL",10000,rMin,rMax);

  int nentries = tree->GetEntries();

  std::cout << "===> Number of entries = " << nentries << std::endl;

  for (Long64_t ievt=0; ievt<tree->GetEntries();ievt++) {    //Long64_t

    if (ievt%1000 == 0){
      std::cout << "--- ... Processing event: " << ievt << std::endl;
    }

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%100000==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					    << int(double(ievt)/double(nentries)*100) << "% done" << std::endl;

    if( ievt==maxNentries ) break;

    tree->GetEntry(ievt);

    h_quantileExpected->Fill( quantileExpected );

    if( quantileExpected<0 ) continue;
    h_deltaNLL->Fill( deltaNLL );
    vec_r_tt.push_back( r_tt );
    vec_deltaNLL.push_back( deltaNLL );

  }


  int n = int(vec_r_tt.size());

  Int_t idx[n];
  double r_tt_sorted[n];
  for( int j=0; j<n; j++ ) r_tt_sorted[j] = vec_r_tt[j];

  TMath::Sort(n,r_tt_sorted,idx,false);


  double x[n], y[n];

  for( int i=0; i<n; i++ ){
    x[i] = vec_r_tt[idx[i]];
    y[i] = vec_deltaNLL[idx[i]];
  }

  TGraph* gr = new TGraph(n, x, y);
  gr->Write("gr_r_tt_deltaNLL");

  treefile->Close();

  histofile.Write();
  histofile.Close();

  std::cout << " Done! " << std::endl;

}
