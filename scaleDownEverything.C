#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TKey.h"
#include "TString.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>

using namespace std;



void scaleDown (TH1* inputHisto,  TString cat, TString sample,  bool doScale) {
  // if (cat != "MVA_ljets_j5_t3" || sample != "ttH125" ||  !doScale)
  //     return;

  cout << "Inside Scale down...." << endl;

  if (!inputHisto) {

    cout << "You gave me a null pointer for histo " << cat << " " << sample << endl;

  }
  
  if (!doScale)
    return;

  cout << "--- Scaling: i = " << inputHisto->Integral();
               

  inputHisto->Scale(5.1/19.5);
  
  //inputHisto->Rebin(inputHisto->GetNbinsX());
  
  cout << " f = " << inputHisto->Integral()
       << endl;
  
  //updateFile->cd();
  inputHisto->Write(inputHisto->GetName(), TFile::kOverwrite);
  //cout << "Wrote out a new histo " << endl;


}



class AccTableRow {

public:   
  string label;
  
  double signalYield;
  double backgroundYield;

  double ratio;

  int bins;

  void calcRatio () {

    if (backgroundYield == 0) {
      cout << "Oops, background yield is 0... set ratio to zero" << endl;
      ratio = 0.0;
      return;
    }
    ratio = signalYield/backgroundYield;
    
  };

  AccTableRow (): label(""), signalYield(0.0), backgroundYield(0.0), ratio(0.0), bins(0)  {};

};




void scaleDownSignal (TString inputName, bool doScale = false) {

  
  //TFile * myFile = new TFile("base_OS_lessSignal.root", "UPDATE");
  //TFile * myFile = new TFile("oneBin_OS.root", "UPDATE");

  TFile * myFile = new TFile (inputName, "UPDATE");

  vector<TString> samples;
  samples.push_back("ttH125");
  samples.push_back("ttbar");
  samples.push_back("ttbarPlusB");
  samples.push_back("ttbarPlusBBbar");
  samples.push_back("ttbarPlusCCbar");
  samples.push_back("singlet");
  samples.push_back("wjets");
  samples.push_back("zjets");
  samples.push_back("ttbarW");
  samples.push_back("ttbarZ");
  samples.push_back("diboson");
  
  
  vector<TString> channelNames;
  //channelNames.push_back ("MVA_e3je2t");
  //channelNames.push_back ("MVA_ge4je2t");

  channelNames.push_back("MVA_TTL_1b_1nb");
  channelNames.push_back("MVA_TTL_1b_2nb");
  channelNames.push_back("MVA_TTL_1b_3+nb");
  channelNames.push_back("MVA_TTL_2b_0nb");
  channelNames.push_back("MVA_TTL_2b_1nb");
  channelNames.push_back("MVA_TTL_2b_2+nb");
  channelNames.push_back ("MVA_ge3t");
  channelNames.push_back ("MVA_e3je2t");
  channelNames.push_back ("MVA_ge4je2t");
  channelNames.push_back ("MVA_ljets_j4_t3");
  channelNames.push_back ("MVA_ljets_j4_t4");
  channelNames.push_back ("MVA_ljets_j5_t3");
  channelNames.push_back ("MVA_ljets_j5_tge4");
  channelNames.push_back ("MVA_ljets_jge6_t2");
  channelNames.push_back ("MVA_ljets_jge6_t3");
  channelNames.push_back ("MVA_ljets_jge6_tge4");
  

  vector<TString> suffixes;
  suffixes.push_back("CMS_ttH_CSVHFStats1");
  suffixes.push_back("CMS_ttH_CSVHFStats2");

  suffixes.push_back("CMS_ttH_CSVLFStats1");
  suffixes.push_back("CMS_ttH_CSVLFStats2");

  suffixes.push_back("CMS_ttH_CSVCErr1");
  suffixes.push_back("CMS_ttH_CSVCErr2");
  

  suffixes.push_back("CMS_ttH_CSVHF");
  suffixes.push_back("CMS_ttH_CSVLF");
  
  suffixes.push_back("CMS_scale_j");
  
  suffixes.push_back("CMS_ttH_topPtcorr");

  suffixes.push_back("Q2scale_ttH_ttbar0p");
  suffixes.push_back("Q2scale_ttH_ttbar1p");
  suffixes.push_back("Q2scale_ttH_ttbar2p");
  suffixes.push_back("Q2scale_ttH_ttbar_b");
  suffixes.push_back("Q2scale_ttH_ttbar_bb");
  suffixes.push_back("Q2scale_ttH_ttbar_cc");

  suffixes.push_back("TES");
  suffixes.push_back("jetTauFake");
  suffixes.push_back("eTauFake");
  suffixes.push_back("tauIdEff");
  
  
  
  //suffixes.push_back("
  
  //suffixes.push_back("PU");

  vector<TString> upDown;
  upDown.push_back("Up");
  upDown.push_back("Down");

  vector<AccTableRow> allRows;
  

    
  for (unsigned iName = 0; iName < channelNames.size(); iName++){
    AccTableRow aRow;
    aRow.label = channelNames[iName];
    
    for (unsigned iSam = 0; iSam < samples.size(); iSam++){

      if (channelNames[iName].Contains("TTL")
          && (samples[iSam] == "ttbarPlusBBbar"
              || samples[iSam] == "ttbarPlusCCbar"
              || samples[iSam] == "ttbarPlusB")) {
        cout << "Skipping " << channelNames[iName] <<  "_" << samples[iSam] << endl;
        continue;
      }
      
      TString baseName = samples[iSam]+ "_" + channelNames[iName];
      cout << "Histo base name is " << baseName << endl;
      
      TH1 * baseHist = (TH1*) myFile->Get(baseName);
      
      if (baseHist==0) {
        cout << "Could not find " << baseName << endl;
        return;
      }

      scaleDown (baseHist, channelNames[iName], samples[iSam], doScale);
      
      aRow.bins = baseHist->GetNbinsX();
      
      if (samples[iSam] != "ttH125") {
        //cout << "--- Histo is background" << endl;
        aRow.backgroundYield +=  baseHist->Integral();
        
      } else {
        //cout << "--- Histo is signal " << endl;
        aRow.signalYield += baseHist->Integral();
      }

      
          
      
      for (unsigned iSuf = 0; iSuf < suffixes.size(); iSuf++){        
        for (unsigned iOther = 0; iOther < upDown.size(); iOther++){
          
          if ( (suffixes[iSuf] == "Q2scale_ttH_ttbar0p")) {
            if ( samples[iSam] != "ttbar" ) {            
              continue;
            }
            if (channelNames[iName] != "MVA_ljets_j4_t3"
                && channelNames[iName] != "MVA_ljets_j4_t4"
                && channelNames[iName] != "MVA_TTL_1b_1nb"
                && channelNames[iName] != "MVA_TTL_2b_0nb" ){              
              continue;
            }
              
          }

          if ( suffixes[iSuf] == "Q2scale_ttH_ttbar1p" ) {

            if ( samples[iSam] != "ttbar" ) {            
              continue;
            }
            if (channelNames[iName] != "MVA_ljets_j5_t3"
                && channelNames[iName] != "MVA_ljets_j5_t4"
                && channelNames[iName] != "MVA_TTL_1b_2nb"
                && channelNames[iName] != "MVA_TTL_2b_1nb"
                && channelNames[iName] != "MVA_ge3t"
                && channelNames[iName] != "MVA_e3je2t"
                ){              
              continue;
            }
            
          }

          if ( suffixes[iSuf] == "Q2scale_ttH_ttbar2p" ) {

            if ( samples[iSam] != "ttbar" ) {            
              continue;
            }
            if (channelNames[iName] != "MVA_ljets_jge6_t2"
                && channelNames[iName] != "MVA_ljets_jge6_t3"
                && channelNames[iName] != "MVA_ljets_jge6_tge4"
                && channelNames[iName] != "MVA_TTL_1b_3+nb"
                && channelNames[iName] != "MVA_TTL_2b_2+nb"
                && channelNames[iName] != "MVA_ge4je2t"){              
              continue;
            }
            
          }

          if ( (suffixes[iSuf] == "CMS_ttH_CSVCErr1" || suffixes[iSuf] == "CMS_ttH_CSVCErr2")
               && samples[iSam] != "ttbarPlusCCbar") {
            continue;
          }
          

          if ( suffixes[iSuf] == "Q2scale_ttH_ttbar_bb"
               && samples[iSam] != "ttbarPlusBBbar" )
            continue;

          if ( suffixes[iSuf] == "Q2scale_ttH_ttbar_b"
               && samples[iSam] != "ttbarPlusB" )
            continue;

          if ( suffixes[iSuf] == "Q2scale_ttH_ttbar_cc"
               && samples[iSam] != "ttbarPlusCCbar" )
            continue;

          if ( suffixes[iSuf] == "TES"
               && !channelNames[iName].Contains("TTL") )
            continue;

          if ( suffixes[iSuf] == "jetTauFake"
               && !channelNames[iName].Contains("TTL"))
            continue;

          if ( suffixes[iSuf] == "tauIdEff"
               && !channelNames[iName].Contains("TTL"))
            continue;

          if ( suffixes[iSuf] == "eTauFake"
               && !channelNames[iName].Contains("TTL"))
            continue;

          TString thisName = samples[iSam]+ "_"
            + channelNames[iName] + "_" + suffixes[iSuf] + upDown[iOther];

          cout << "--- This is histo " << thisName
               << endl;
    
          TH1* pHisto = (TH1*) myFile->Get(thisName);


          scaleDown(pHisto, channelNames[iName], samples[iSam], doScale);
          
          
          delete pHisto;
      
        }
      }
    }
    aRow.calcRatio();
    allRows.push_back(aRow);
  }
  cout << "Done with loop over names" << endl;
  cout << "==================== Table ===================" << endl;

  for (unsigned iRow = 0; iRow < allRows.size(); iRow++){
    cout << allRows[iRow].label
         << " signal = " << allRows[iRow].signalYield
         << " bkg = " << allRows[iRow].backgroundYield
         << " bins = " << allRows[iRow].bins
         << endl
         << "---ratio = " << allRows[iRow].ratio
         << endl;
  }
  
  myFile->Close();
}
