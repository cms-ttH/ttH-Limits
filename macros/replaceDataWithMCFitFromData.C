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
#include "TString.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooRealVar.h"
#include "RooRandom.h"

/*

Starting with the combined root file
LJ_OSDIL_TAU_MVA.root 

Use as

cp LJ_OSDIL_TAU_MVA.root LJ_OSDIL_TAU_MVA_replaceDataWithMCFitFromData_postFitS_ttH.root

root -b -q head.C replaceDataWithMCFitFromData.C+'("LJ_OSDIL_TAU_MVA_replaceDataWithMCFitFromData_postFitS_ttH.root","postFitS","mlfit_LJ_OSDIL_TAU.root","wsTest_LJ_OSDIL_TAU.root",1.)'

Your new file to use with datacard maker and limits and fits is
LJ_OSDIL_TAU_MVA_replaceDataWithMCFitFromData_postFitS_ttH.root


*/

class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

void replaceDataWithMCFitFromData( TString dataFileName = "", TString fitTypeName = "postFitS", TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root", double scaleTTH=1.0 ) {

  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open( wsFileName );
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open( fitFileName );

  TFile * myFile = new TFile(dataFileName, "UPDATE");

  TString fitLabel = fitTypeName;
  RooFitResult *fitRes;  
  if( fitTypeName.EqualTo("postFitB") )      fitRes = (RooFitResult*)fitFile->Get("fit_b");
  else if( fitTypeName.EqualTo("postFitS") ) fitRes = (RooFitResult*)fitFile->Get("fit_s");
  else                                       fitRes = (RooFitResult*)fitFile->Get("nuisances_prefit_res");


  //Try making "snapshots" with preFit and postFit nuisance values
  if( fitTypeName.EqualTo("postFitB") )      w->saveSnapshot("postfitB",RooArgSet(fitRes->floatParsFinal()),true);
  else if( fitTypeName.EqualTo("postFitS") ) w->saveSnapshot("postfitS",RooArgSet(fitRes->floatParsFinal()),true);
  else                                       w->saveSnapshot("prefit",RooArgSet(fitRes->floatParsFinal()),true);

  if( fitTypeName.EqualTo("postFitB") )      w->loadSnapshot("postfitB");
  else if( fitTypeName.EqualTo("postFitS") ) w->loadSnapshot("postfitS");
  else                                       w->loadSnapshot("prefit");

  myFile->cd();
  
  TString discName = "MVA";
  const int numChannels = 16;
  TString channelNames[] = {

    "ljets_j4_t3",
    "ljets_j4_t4",
	
    "ljets_j5_t3",
    "ljets_j5_tge4",

    "ljets_jge6_t2",
    "ljets_jge6_t3",
    "ljets_jge6_tge4",

    "ge3t",		
    "e3je2t",
    "ge4je2t",

    "TTL_1b_1nb",
    "TTL_1b_2nb",
    "TTL_1b_3+nb",
    "TTL_2b_0nb",
    "TTL_2b_1nb",
    "TTL_2b_2+nb"
  };

  int numSamples = 10;//11;
  // TString sampleNames[] = {"ttH125","ttbar","ttbarPlusB","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson"};
  TString sampleNames[] = {"ttbar","ttbarPlusB","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson"};

      
  RooAbsPdf* nuisancePdf =  w->pdf("nuisancePdf");

  //  first, get the data_obs histogram then delete it
  //  then create a new data_obs histogram and replace it with   
  for ( int iChan = 0; iChan < numChannels; iChan++ ){

    TString catName = channelNames[iChan];
    TString dataName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan] + TString(";*");
    TString dataRealName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan];
    TH1 * checkData = (TH1*) myFile->Get(dataRealName);
    TH1* dataHist = (TH1D*)checkData->Clone(dataRealName+"_temp");
    double theNorm = checkData->Integral();
    cout << "The data  " << channelNames[iChan]  << "  norm is  " << theNorm  << endl;
    myFile->Delete(dataName);
    //myFile->Write();
    TH1 * sumHisto;
    double channelSumTotal = 0;

    for (int iSample = 0; iSample < numSamples; iSample++) {
      TString histoName = sampleNames[iSample] + "_" + discName + "_" +  channelNames[iChan];

      TString procName = sampleNames[iSample];
      //Construct the name of the normalization variable.
      TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;

      RooAbsReal *fitN = w->function(normName);

      TH1 *bkgTempHist;
      if( fitN ){
	RooAbsPdf* bkgPdf = w->pdf("shapeBkg_"+catName+"_"+procName+"_morph");
	RooProdPdf prodTemp("prodTemp_"+catName+"_"+procName+"_"+fitLabel,"prodTemp_"+catName+"_"+procName+"_"+fitLabel,*nuisancePdf,*bkgPdf);
	bkgTempHist = prodTemp.createHistogram("CMS_th1x");
	bkgTempHist->Scale(fitN->getVal());
      }
      else {
	bkgTempHist = (TH1*)dataHist->Clone("pdf_bin"+catName+"_"+procName+"_"+fitLabel+"CloneTemp");
	bkgTempHist->Reset();
      }

      TH1* tempHisto = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+procName+"_"+fitLabel+"Clone");
      tempHisto->Reset();
      for (int iBin = 0; iBin < tempHisto->GetNbinsX(); iBin++) {
	tempHisto->SetBinContent(iBin+1,bkgTempHist->GetBinContent(iBin+1));
      }

      if ( iSample == 0 ){
        TString cloneName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan];
        sumHisto = (TH1*) tempHisto->Clone(cloneName);
        channelSumTotal += tempHisto->Integral();
        cout << histoName << "   norm =  " << tempHisto->Integral() << ", now sum is " << channelSumTotal << endl;
      } else {
        sumHisto->Add(tempHisto);
        channelSumTotal += tempHisto->Integral();
        cout << histoName << "   norm =  " << tempHisto->Integral() << ", now sum is " << channelSumTotal << endl;
      }

      delete bkgTempHist;
    }

    TString ttHhistoName = "ttH125_" + discName + "_" +  channelNames[iChan];
    TH1 * ttH_input = (TH1*) myFile->Get(ttHhistoName);
    double SM_ttH_input = ttH_input->Integral();

    TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";
    RooAbsReal *fitTTH = w->function(normName_ttH);

    TH1* sigTempHist = w->pdf("shapeSig_"+catName+"_ttH_morph")->createHistogram("CMS_th1x");
    // sigTempHist->Scale(fitTTH->getVal());
    sigTempHist->Scale( SM_ttH_input/sigTempHist->Integral() );

    TH1* h_ttH = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_ttH_"+fitLabel+"Clone");
    h_ttH->Reset();
    for (int iBin = 0; iBin < h_ttH->GetNbinsX(); iBin++) {
      h_ttH->SetBinContent(iBin+1,sigTempHist->GetBinContent(iBin+1));
    }

    h_ttH->Scale( scaleTTH );

    sumHisto->Add(h_ttH);
    channelSumTotal += h_ttH->Integral();
    cout << ttHhistoName << "   norm =  " << h_ttH->Integral() << ", now sum is " << channelSumTotal << endl;


    //Now, fix the errors so they look like data
    for (int iBin = 1; iBin <= sumHisto->GetNbinsX(); ++iBin) {
      double n = sumHisto->GetBinContent(iBin);
      sumHisto->SetBinError(iBin,sqrt(n));
    }

    cout << "Channel " << channelNames[iChan] << " norm =  " <<  sumHisto->Integral() << endl;
    sumHisto->SetDirectory(myFile);
    sumHisto->Write();

    
  }

  std::cout << "Done!" << std::endl;

}
