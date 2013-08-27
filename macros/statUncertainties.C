#include "TString.h"

void statUncertainties(TString inFileName, TString outFileName, bool is8TeV=true ) {

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

  //Open the root file
  TFile *inFile = TFile::Open(inFileName);
  TFile *outFile = TFile::Open(outFileName,"RECREATE");

  int numSamples = 16;
  TString sampleNames[] = {"ttbar","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson",
			   "ttH110","ttH115","ttH120","ttH125","ttH130","ttH135","ttH140"};


  for (int iChan = 0; iChan < numChannels; ++iChan) {

    for (int iSample = 0; iSample < numSamples; ++iSample) {

      TString histName = sampleNames[iSample];
      histName += "_CFMlpANN_";
      histName += channelNames[iChan];

      std::cout << histName << " " << std::endl;

      TH1D* hist = (TH1D*)inFile->Get(histName);

      int numBins = hist->GetNbinsX();

      TString era = ( is8TeV ) ? "8TeV" : "7TeV";

      for( int bin=0; bin<numBins; bin++ ){
	TString bin_name;
	bin_name.Form("%d",bin);

	TString histName_Up = histName;
	histName_Up += "_" + sampleNames[iSample] + "_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name + "Up";

	TString histName_Down = histName;
	histName_Down += "_" + sampleNames[iSample] + "_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name + "Down";

	if( sampleNames[iSample].Contains("ttH") ){
	  histName_Up = histName;
	  histName_Up += "_ttH_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name + "Up";

	  histName_Down = histName;
	  histName_Down += "_ttH_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name + "Down";
	}

	TH1D* hist_Up   = (TH1D*)hist->Clone(histName_Up);
	TH1D* hist_Down = (TH1D*)hist->Clone(histName_Down);

	hist_Up->SetBinContent(bin+1, hist->GetBinContent(bin+1) + hist->GetBinError(bin+1) );
	hist_Down->SetBinContent(bin+1, hist->GetBinContent(bin+1) - hist->GetBinError(bin+1) );

	outFile->cd();

	hist_Up->Write(histName_Up);
	hist_Down->Write(histName_Down);
      }
    }
  }


  //Now copy all other histograms in the file as is (preserving names...)
  TList *keys = inFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();

    TH1 *hist2 = 0;
    hist2 = (TH1 *)inFile->Get(name);

    outFile->cd();
    hist2->Write(name);
  }


  outFile->Close();

  std::cout << "Done." << std::endl;

}
