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

void separate_btag_rate_shape(TString histFileName = "ttH_ljets_7TeV_2012_06_06.root",
			      TString outFileName =  "ttH_ljets_scale_7TeV_2012_06_06.root") {


  TFile *outFile = TFile::Open(outFileName,"RECREATE");

  TFile *histFile = TFile::Open(histFileName);


  int numChannels = 15;
  TString channelNames[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    "ljets_j4_t3", 		//	 "j4_3t",
    "ljets_j4_t4", 		//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    "ljets_j5_t3", 		//	 "j5_3t",
    "ljets_j5_tge4", 		//	 "j5_3t",
    "ljets_jge6_t2", 		//	 "j6_2t",
    "ljets_jge6_t3", 		//	 "j6_3t",
    "ljets_jge6_tge4",    	//	 "j6_3t",
    "ge3t",	       		// dilep >= 3tags
    "e2je2t",   		// dilep 2 jets 2 tags
    "e3je2t",   		// dilep 3 jets 2 tags
    "ge4je2t",   		// dilep >=4 jets 2 tags
    "SS_e3je1t",   		// dilep SS 3 jets 1 tags
    "SS_ge4je1t",   		// dilep SS >=4 jets 1 tags
    "SS_e3jge2t",   		// dilep SS 3 jets >=2 tags
    "SS_ge4jge2t",   		// dilep SS >=4 jets >=2 tags
  };

  int numSamples = 16;
  TString sampleNames[] = {"ttbar","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson",
			   "ttH110","ttH115","ttH120","ttH125","ttH130","ttH135","ttH140"};
  
  // std::cout << "channel  ";
  // for (int iSample = 0; iSample < numSamples; ++iSample) {
  //   std::cout << " & " << sampleNames[iSample];
  // }
  // std::cout << " \\\\ " << std::endl;

  //std::cout << " channel & signal & background & S/B & S/sqrt(B) \\\\ " << std::endl;

  TH1::SetDefaultSumw2();

  for (int iChan = 0; iChan < numChannels; ++iChan) {

    std::cout << channelNames[iChan] << " ";

    for (int iSample = 0; iSample < numSamples; ++iSample) {

      TString histName = sampleNames[iSample];
      histName += "_CFMlpANN_";
      histName += channelNames[iChan];

      std::cout << histName << " " << std::endl;

      TString histName_eff_bUp = histName;
      histName_eff_bUp += "_CMS_eff_bUp";

      TString histName_eff_bDown = histName;
      histName_eff_bDown += "_CMS_eff_bDown";

      TString histName_fake_bUp = histName;
      histName_fake_bUp += "_CMS_fake_bUp";

      TString histName_fake_bDown = histName;
      histName_fake_bDown += "_CMS_fake_bDown";

      ///// new names
      TString new_histName_eff_bUp = histName;
      new_histName_eff_bUp += "_CMS_eff_bShapeUp";
      TString new_histName_eff_bDown = histName;
      new_histName_eff_bDown += "_CMS_eff_bShapeDown";
      TString new_histName_fake_bUp = histName;
      new_histName_fake_bUp += "_CMS_fake_bShapeUp";
      TString new_histName_fake_bDown = histName;
      new_histName_fake_bDown += "_CMS_fake_bShapeDown";

      TString new_histName_eff_bRateUp = histName;
      new_histName_eff_bRateUp += "_CMS_eff_bRateUp";
      TString new_histName_eff_bRateDown = histName;
      new_histName_eff_bRateDown += "_CMS_eff_bRateDown";
      TString new_histName_fake_bRateUp = histName;
      new_histName_fake_bRateUp += "_CMS_fake_bRateUp";
      TString new_histName_fake_bRateDown = histName;
      new_histName_fake_bRateDown += "_CMS_fake_bRateDown";


      TH1D* hist = (TH1D*)histFile->Get(histName);

      TH1D* hist_eff_bRateUp = (TH1D*)hist->Clone(new_histName_eff_bRateUp);
      TH1D* hist_eff_bRateDown = (TH1D*)hist->Clone(new_histName_eff_bRateDown);
      TH1D* hist_fake_bRateUp = (TH1D*)hist->Clone(new_histName_fake_bRateUp);
      TH1D* hist_fake_bRateDown = (TH1D*)hist->Clone(new_histName_fake_bRateDown);

      TH1D* hist_eff_bUp_temp = (TH1D*)histFile->Get(histName_eff_bUp);
      TH1D* hist_eff_bDown_temp = (TH1D*)histFile->Get(histName_eff_bDown);
      TH1D* hist_fake_bUp_temp = (TH1D*)histFile->Get(histName_fake_bUp);
      TH1D* hist_fake_bDown_temp = (TH1D*)histFile->Get(histName_fake_bDown);

      TH1D* hist_eff_bUp = (TH1D*)hist_eff_bUp_temp->Clone(new_histName_eff_bUp);
      TH1D* hist_eff_bDown = (TH1D*)hist_eff_bDown_temp->Clone(new_histName_eff_bDown);
      TH1D* hist_fake_bUp = (TH1D*)hist_fake_bUp_temp->Clone(new_histName_fake_bUp);
      TH1D* hist_fake_bDown = (TH1D*)hist_fake_bDown_temp->Clone(new_histName_fake_bDown);

      double sum = hist->Integral();

      double sum_eff_bUp = hist_eff_bUp->Integral();
      double sum_eff_bDown = hist_eff_bDown->Integral();
      double sum_fake_bUp = hist_fake_bUp->Integral();
      double sum_fake_bDown = hist_fake_bDown->Integral();

      double scale_eff_bRateUp = ( sum>0. ) ? sum_eff_bUp/sum : 1.;
      double scale_eff_bRateDown = ( sum>0. ) ? sum_eff_bDown/sum : 1.;
      double scale_fake_bRateUp = ( sum>0. ) ? sum_fake_bUp/sum : 1.;
      double scale_fake_bRateDown = ( sum>0. ) ? sum_fake_bDown/sum : 1.;

      double scale_eff_bUp = ( sum_eff_bUp>0. ) ? sum/sum_eff_bUp : 1.;
      double scale_eff_bDown = ( sum_eff_bDown>0. ) ? sum/sum_eff_bDown : 1.;
      double scale_fake_bUp = ( sum_fake_bUp>0. ) ? sum/sum_fake_bUp : 1.;
      double scale_fake_bDown = ( sum_fake_bDown>0. ) ? sum/sum_fake_bDown : 1.;

      hist_eff_bUp->Scale(scale_eff_bUp);
      hist_eff_bDown->Scale(scale_eff_bDown);
      hist_fake_bUp->Scale(scale_fake_bUp);
      hist_fake_bDown->Scale(scale_fake_bDown);

      hist_eff_bRateUp->Scale(scale_eff_bRateUp);
      hist_eff_bRateDown->Scale(scale_eff_bRateDown);
      hist_fake_bRateUp->Scale(scale_fake_bRateUp);
      hist_fake_bRateDown->Scale(scale_fake_bRateDown);


      outFile->cd();

      hist_eff_bUp->Write(new_histName_eff_bUp);

      hist_eff_bUp->Write(new_histName_eff_bUp);
      hist_eff_bDown->Write(new_histName_eff_bDown);
      hist_fake_bUp->Write(new_histName_fake_bUp);
      hist_fake_bDown->Write(new_histName_fake_bDown);

      hist_eff_bRateUp->Write(new_histName_eff_bRateUp);
      hist_eff_bRateDown->Write(new_histName_eff_bRateDown);
      hist_fake_bRateUp->Write(new_histName_fake_bRateUp);
      hist_fake_bRateDown->Write(new_histName_fake_bRateDown);
    }
  }



  //Now copy all other histograms in the file as is (preserving names...)
  TList *keys = histFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();

    TH1 *hist2 = 0;
    hist2 = (TH1 *)histFile->Get(name);

    outFile->cd();
    hist2->Write(name);
  }


  outFile->Close();

  std::cout << "Done." << std::endl;

}
