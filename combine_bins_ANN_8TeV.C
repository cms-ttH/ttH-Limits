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

void combine_bins_ANN_8TeV(TString histFileName = "ttH_ljets_7TeV_2012_06_06.root",
		           TString outFileName = "ttH_ljets_scale_7TeV_2012_06_06.root") {


  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TFile *histFile = TFile::Open(histFileName);

  std::cout << "Copying over all histograms..." << std::endl;

  int numChannels = 9;
  TString channelNames[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    "ljets_j4_t3", 		//	 "j4_3t",
    "ljets_j4_t4", 		//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    "ljets_j5_t3", 		//	 "j5_3t",
    "ljets_j5_tge4", 		//	 "j5_4t",
    "ljets_jge6_t2", 		//	 "j6_2t",
    "ljets_jge6_t3", 		//	 "j6_3t",
    "ljets_jge6_tge4",    	//	 "j6_3t",
    "ge3t",	       		// dilep >= 3tags
    "e2je2t"   			// dilep 2 jets 2 tags
  };

  int nbins_new[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    19, 		//	 "j4_3t",
    8, 	 	//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    18, 		//	 "j5_3t",
    7, 		//	 "j5_4t",
    19, 		//	 "j6_2t",
    18, 		//	 "j6_3t",
    9,    	//	 "j6_3t",
    9,	       		// dilep >= 3tags
    9   			// dilep 2 jets 2 tags
  };

  int firstbin_new[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    2, 		//	 "j4_3t",
    2, 		//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    2, 		//	 "j5_3t",
    2, 		//	 "j5_4t",
    2, 		//	 "j6_2t",
    2, 		//	 "j6_3t",
    2,    	//	 "j6_4t",
    2,	       		// dilep >= 3tags
    2   			// dilep 2 jets 2 tags
  };

  int lastbin_new[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    20, 		//	 "j4_3t",
    9, 	 	//	 "j4_4t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    19, 		//	 "j5_3t",
    8, 		//	 "j5_4t",
    20, 		//	 "j6_2t",
    19, 		//	 "j6_3t",
    10,    	//	 "j6_4t",
    10,	       		// dilep >= 3tags
    10   			// dilep 2 jets 2 tags
  };

  //Now copy all other histograms in the file as is (preserving names...)
  TList *keys = histFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();
    TString name_temp = name + "_Clone";

    TString new_name = name.Copy();

    //if( !(name.Contains("ttbar_")&&name.Contains("ljets_jge6_tge4")) ) continue;

    //std::cout << " ==> name = " << name << std::endl;

    TH1 *hist_temp = 0;
    hist_temp = (TH1 *)histFile->Get(name)->Clone(name_temp);

    double integral_old = hist_temp->Integral();

    TH1F *hist = 0;

    bool doRebin = true;

    if( doRebin ){
      bool someone_caught_me = false;
      for (int iChan = 0; iChan < numChannels; ++iChan) {
	if( !name.Contains(channelNames[iChan]) ) continue;
      
	someone_caught_me = true;

	int firstbin = firstbin_new[iChan];
	int lastbin  = lastbin_new[iChan];
	int nbins = nbins_new[iChan];

	double xmin = hist_temp->GetBinLowEdge(firstbin);
	double xmax = hist_temp->GetBinLowEdge(lastbin) + hist_temp->GetBinWidth(lastbin);

	int nbins_old = hist_temp->GetNbinsX();

	//std::cout << " nbins_old = " << nbins_old << " nbins = " << nbins << " xmin = " << xmin << " xmax = " << xmax << std::endl; 
	hist = new TH1F("hist","", nbins, xmin, xmax);
	//std::cout << " first bin = " << lastbin << std::endl;
	if( true ){
	  double sum = 0;
	  double err2 = 0;
	  for( int iBin=1; iBin<=firstbin; iBin++ ){
	    sum  += hist_temp->GetBinContent(iBin);
	    err2 += hist_temp->GetBinError(iBin) * hist_temp->GetBinError(iBin);
	    //std::cout << "\t iBin = " << iBin << " sum = " << sum << " err2 = " << err2 << std::endl;
	  }
	  hist->SetBinContent(1,sum);
	  hist->SetBinError(1,sqrt(err2));
	}

	//std::cout << " last bin = " << lastbin << std::endl;
	if( true ){
	  double sum = 0;
	  double err2 = 0;
	  for( int iBin=lastbin; iBin<=nbins_old; iBin++ ){
	    sum  += hist_temp->GetBinContent(iBin);
	    err2 += hist_temp->GetBinError(iBin) * hist_temp->GetBinError(iBin);
	    //std::cout << "\t iBin = " << iBin << " sum = " << sum << " err2 = " << err2 << std::endl;
	  }
	  hist->SetBinContent(nbins,sum);
	  hist->SetBinError(nbins,sqrt(err2));
	}

	for( int iBin=firstbin+1; iBin<lastbin; iBin++ ){
	  hist->SetBinContent(iBin-firstbin+1,hist_temp->GetBinContent(iBin));
	  hist->SetBinError(iBin-firstbin+1,hist_temp->GetBinError(iBin));
	}
      }
    }

    if( !doRebin || !someone_caught_me ) hist = (TH1F*)hist_temp->Clone(new_name);

    //if( !someone_caught_me ) continue;

    double integral_new = hist->Integral();

    outFile->cd();
    hist->Write(new_name);

    delete hist;
  }

  outFile->Close();

  std::cout << "Done." << std::endl;

}
