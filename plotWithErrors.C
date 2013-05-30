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


class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

void plotWithErrors( TString dataFileName = "", TString catName = "", bool is8TeV=true, int nToys=500, bool debug=false ) {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  //Now we need a list of categories and a list of processes
  TList processes;
  processes.Add(new LabelInfo("ttbar","$t\\bar{t}+$LF"));
  processes.Add(new LabelInfo("ttbarPlusBBbar","$t\\bar{t}+b\\bar{b}$"));
  processes.Add(new LabelInfo("ttbarPlusCCbar","$t\\bar{t}+c\\bar{c}$"));
  processes.Add(new LabelInfo("ttbarW","$t\\bar{t}W$"));
  processes.Add(new LabelInfo("ttbarZ","$t\\bar{t}Z$"));
  processes.Add(new LabelInfo("singlet","Single $t$"));
  processes.Add(new LabelInfo("wjets","$W+$jets"));
  processes.Add(new LabelInfo("zjets","$Z+$jets"));
  processes.Add(new LabelInfo("diboson","Diboson"));

  int NumBkgs = 9;

  int bin_ttjets_other = 0;
  int bin_ttjets_bbbar = 1;
  int bin_ttjets_ccbar = 2;
  int bin_ttW     = 3;
  int bin_ttZ     = 4;
  int bin_singlet = 5;
  int bin_wjets   = 6;
  int bin_zjets   = 7;
  int bin_diboson = 8;

  int bin_ewk = 9;
  int bin_ttV = 10;

  int bin_ttH = 11;


  Color_t color[22];
  color[0] = kBlack;

  color[bin_diboson] = kCyan;
  color[bin_zjets]   = kGreen+2;
  color[bin_wjets]   = kAzure+1;
  color[bin_singlet] = kOrange+7;
  color[bin_ttW]     = kBlue-9;
  color[bin_ttZ]     = color[bin_ttW];

  color[bin_ttjets_bbbar]  = kMagenta+2;
  color[bin_ttjets_ccbar]  = kOrange-7;
  color[bin_ttjets_other]  = kRed;

  color[bin_ewk] = kAzure+1;
  color[bin_ttV] = kBlue-9;

  color[bin_ttH]     = kMagenta;


  //These are the list of possible categories
  TList allCategories;
  allCategories.Add(new LabelInfo("ljets_jge6_t2","Lepton + #geq6 jets + 2 tags"));
  allCategories.Add(new LabelInfo("ljets_j4_t3","Lepton + 4 jets + 3 tags"));
  allCategories.Add(new LabelInfo("ljets_j5_t3","Lepton + 5 jets + 3 tags"));
  allCategories.Add(new LabelInfo("ljets_jge6_t3","Lepton + #geq6 jets + 3 tags"));
  allCategories.Add(new LabelInfo("ljets_j4_t4","Lepton + 4 jets + 4 tags"));
  allCategories.Add(new LabelInfo("ljets_j5_tge4","Lepton + 5 jets + #geq4 tags"));
  allCategories.Add(new LabelInfo("ljets_jge6_tge4","Lepton + #geq6 jets + #geq4 tags"));
  allCategories.Add(new LabelInfo("e2je2t","Dilepton + 2 jets + 2 tags"));
  allCategories.Add(new LabelInfo("ge3t","Dilepton + #geq3 tags"));


  //These are the list of possible categories
  TList fitResults;
  fitResults.Add(new LabelInfo("nuisances_prefit_res","preFit"));
  fitResults.Add(new LabelInfo("fit_b","postFitB"));
  fitResults.Add(new LabelInfo("fit_s","postFitS"));


  TString era = ( is8TeV ) ? "8TeV" : "7TeV";

  std::string cmsinfo = "CMS Preliminary,  #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}";
  if( is8TeV ) cmsinfo = "CMS Preliminary,  #sqrt{s} = 8 TeV, L = 5.1 fb^{-1}";
  TLatex CMSInfoLatex(0.47, 0.91, cmsinfo.c_str());
  CMSInfoLatex.SetNDC();
  CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.035);

  //Suppress the printout from making plots
  //gErrorIgnoreLevel = 2000;

  //RooFit::RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling) ;
 
  //gSystem->Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so");

  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open("wsTest.root");
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open("results_mlfit.root");
  TFile *dataFile = TFile::Open(dataFileName);


  double ratioMax = 2.3;
  double ratioMin = 0.0;

  TIter nextFit(&fitResults);
  LabelInfo *fitRes = 0;

  while ((fitRes = (LabelInfo *)nextFit())) {

    TString fitName = fitRes->name;
    TString fitLabel = fitRes->label;

    std::cout << "  ==> Category = " << catName << std::endl;
    std::cout << "\t doing " << fitLabel << std::endl;

    RooFitResult *fitFR = (RooFitResult*)fitFile->Get( fitName );

    std::vector<TString> nuisance_names;
    nuisance_names.clear();

    RooArgSet myArgs = fitFR->floatParsFinal();

    TIterator *nextArg(myArgs.createIterator());
    for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
      //RooRealVar *rrv = a;      
      if( debug ) std::cout << " name = " << rrv->GetName() << ":\t old value = " << rrv->getVal() << ",\t error = " << rrv->getError() << std::endl;
      TString nuis_name = rrv->GetName();
      nuisance_names.push_back( nuis_name );
    }


    std::cout << " Number of nuisances = " << int(nuisance_names.size()) << std::endl;

    int NumPars = int(nuisance_names.size());

    std::string fitinfo = "";
    if( fitLabel.EqualTo("preFit") )        fitinfo = "Pre-Fit";
    else if( fitLabel.EqualTo("postFitB") ) fitinfo = "Post-Fit (B)";
    else if( fitLabel.EqualTo("postFitS") ) fitinfo = "Post-Fit (S+B)";


    TH1D* h_nuis[NumPars];
    for( int iPar=0; iPar<NumPars; iPar++ ){
      TString nuis_name = nuisance_names[iPar];
      h_nuis[iPar] = new TH1D( "h_nuis_"+nuis_name+"_"+fitLabel,"", 50, -4, 4 );
      h_nuis[iPar]->GetXaxis()->SetTitle(nuis_name);
    }


    TLatex FitinfoLatex(0.45, 0.95, fitinfo.c_str());
    FitinfoLatex.SetNDC();
    FitinfoLatex.SetTextFont(42);
    FitinfoLatex.SetTextSize(0.035);

    TIter nextCat(&allCategories);
    LabelInfo *category = 0;
    TString catLabel = "";
    while ((category = (LabelInfo *)nextCat())) {
      TString temp_catName = category->name;
      if( catName.EqualTo(temp_catName) ){
	catLabel = category->label;
	break;
      }
    }


    TString selectioninfo = catLabel;//"Lepton + #geq6 jets + 2 tags";
    TLatex SELECTIONInfoLatex(0.14, 0.91, selectioninfo);
    SELECTIONInfoLatex.SetNDC();
    SELECTIONInfoLatex.SetTextFont(42);
    SELECTIONInfoLatex.SetTextSize(0.035);


    TH1 *dataHist = (TH1 *)dataFile->Get("data_obs_CFMlpANN_"+catName);
    dataHist->SetLineColor(1);
    dataHist->SetLineWidth(2);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(0.75);

    int numBins = dataHist->GetNbinsX();  

    for( int iBin=0; iBin<numBins; iBin++ ){
      if( dataHist->GetBinContent(iBin+1)<0.9 ) dataHist->SetBinContent(iBin+1,0.);
    }

    //Now get the fit central value
    w->saveSnapshot(fitLabel,RooArgSet(fitFR->floatParsFinal()),true);
    w->loadSnapshot(fitLabel);


    //Calculate the total in the fit...
    RooArgSet fitArgs;
   
    TH1* h_bkg[NumBkgs];

    RooAbsPdf* nuisancePdf =  w->pdf("nuisancePdf");

    TIter nextProcess(&processes);
    LabelInfo *process = 0;
    int nProc = 0;
    std::cout << "\t Filling individual sample components" << std::endl;
    while ((process = (LabelInfo *)nextProcess())) {

      TString procName = process->name;
      //Construct the name of the normalization variable.
      TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;

      RooAbsReal *fitN = w->function(normName);

      if (fitN) {
	fitArgs.add(*fitN);
      }

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

      h_bkg[nProc] = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+procName+"_"+fitLabel+"Clone");
      h_bkg[nProc]->Reset();
      for (int iBin = 1; iBin <= numBins; ++iBin) {
	h_bkg[nProc]->SetBinContent(iBin,bkgTempHist->GetBinContent(iBin));
      }

      delete bkgTempHist;
      nProc++;
    }

    TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";
    RooAbsReal *fitTTH = w->function(normName_ttH);


    TH1* sigTempHist = w->pdf("shapeSig_"+catName+"_ttH_morph")->createHistogram("CMS_th1x");
    sigTempHist->Scale(fitTTH->getVal());

    TH1* h_ttH = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_ttH_"+fitLabel+"Clone");
    h_ttH->Reset();
    for (int iBin = 1; iBin <= numBins; ++iBin) {
      h_ttH->SetBinContent(iBin,sigTempHist->GetBinContent(iBin));
    }


    delete sigTempHist;

    h_ttH->SetLineColor(color[bin_ttH]);
    h_ttH->SetLineWidth(2);


    TH1D* h_ewk_bkg = (TH1D*)h_bkg[bin_wjets]->Clone("h_ewk_"+catName+"_"+fitLabel);
    h_ewk_bkg->Add(h_bkg[bin_zjets]);
    h_ewk_bkg->Add(h_bkg[bin_diboson]);

    TH1D* h_ttbarV_bkg = (TH1D*)h_bkg[bin_ttW]->Clone("h_ttbarV_"+catName+"_"+fitLabel);
    h_ttbarV_bkg->Add(h_bkg[bin_ttZ]);

    h_ewk_bkg->SetFillColor(color[bin_ewk]);
    h_ttbarV_bkg->SetFillColor(color[bin_ttV]);
    h_bkg[bin_singlet]->SetFillColor(color[bin_singlet]);
    h_bkg[bin_ttjets_other]->SetFillColor(color[bin_ttjets_other]);
    h_bkg[bin_ttjets_ccbar]->SetFillColor(color[bin_ttjets_ccbar]);
    h_bkg[bin_ttjets_bbbar]->SetFillColor(color[bin_ttjets_bbbar]);

    h_ewk_bkg->SetLineColor(color[bin_ewk]);
    h_ttbarV_bkg->SetLineColor(color[bin_ttV]);
    h_bkg[bin_singlet]->SetLineColor(color[bin_singlet]);
    h_bkg[bin_ttjets_other]->SetLineColor(color[bin_ttjets_other]);
    h_bkg[bin_ttjets_ccbar]->SetLineColor(color[bin_ttjets_ccbar]);
    h_bkg[bin_ttjets_bbbar]->SetLineColor(color[bin_ttjets_bbbar]);


    THStack *hs = new THStack("hs"+catName+"_"+fitLabel,"");
    hs->Add(h_ewk_bkg);
    hs->Add(h_ttbarV_bkg);
    hs->Add(h_bkg[bin_singlet]);
    hs->Add(h_bkg[bin_ttjets_bbbar]);
    hs->Add(h_bkg[bin_ttjets_ccbar]);
    hs->Add(h_bkg[bin_ttjets_other]);


    RooAddition fitTotal("fitTotal","fitTotal",fitArgs);

    //Make a normalized plot...
    TH1 *fitTempHist = w->pdf("pdf_bin"+catName+"_bonly")->createHistogram("CMS_th1x");
    fitTempHist->Scale(fitTotal.getVal());

    //Roofit plots histogram vs "bin number;" turn this into a histogram over ANN output
    TH1 *fitHist = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+fitLabel+"Clone");
    fitHist->Reset();
    for (int iBin = 1; iBin <= numBins; ++iBin) {
      fitHist->SetBinContent(iBin,fitTempHist->GetBinContent(iBin));
    }
  

    fitHist->SetLineColor(kGreen+1);
    fitHist->SetLineWidth(2);
    fitHist->SetTitle(";;Events");
  
    //Just to be safe, ditch the temp histogram
    delete fitTempHist;
  
    //Now we need to fluctuate the nuisance parameters and get the fluctuated histograms
//     TH2F *fluctHist = new TH2F("fluctHist","",
// 			       numBins,dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax(),
// 			       1000, fitHist->GetMinimum()/2, fitHist->GetMaximum()*2);
    // TH2F *fluctHist = new TH2F("fluctHist","",
    // 			       numBins,dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax(),
    // 			       100000, 0, fitHist->GetMaximum()*2);

    TH1D* fluctHist[numBins];
    //for( int iBin=0; iBin<numBins; iBin++ ) fluctHist[iBin] = new TH1D(Form("fluctHist_%d",iBin),"",300,0,5*fitHist->GetBinContent(iBin+1));
    for( int iBin=0; iBin<numBins; iBin++ ) fluctHist[iBin] = new TH1D(Form("fluctHist_%d_",iBin)+catName+"_"+fitLabel,"",10000,0,4*fitHist->GetBinContent(iBin+1));

    std::cout << "\t Starting toy creation to find uncertainty band" << std::endl;

    int nSample = nToys;
    for (int i = 0; i < nSample; ++i) {
      RooArgSet myArgs_temp = fitFR->randomizePars();

      int iPar=0;
      TIterator *nextArg_temp(myArgs_temp.createIterator());
      for (TObject *a = nextArg_temp->Next(); a != 0; a = nextArg_temp->Next()) {
	RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
	TString nuis_name = rrv->GetName();
	if( !nuis_name.EqualTo(nuisance_names[iPar]) ) assert(0);
	h_nuis[iPar]->Fill(rrv->getVal());
	iPar++;
      }

      w->saveSnapshot("fluctuation"+fitLabel,RooArgSet(myArgs_temp),true);
      w->loadSnapshot("fluctuation"+fitLabel);
    
      TH1 *temp = w->pdf("pdf_bin"+catName)->createHistogram("CMS_th1x");
      temp->Scale(fitTotal.getVal());
    
      // for (int iBin=0; iBin<numBins; ++iBin) {
      // 	double x = dataHist->GetBinCenter(iBin);
      // 	double y = temp->GetBinContent(iBin);
      
      // 	fluctHist->Fill(x,y);
      // }

      for (int iBin=0; iBin<numBins; ++iBin) fluctHist[iBin]->Fill(std::min(temp->GetBinContent(iBin+1),fluctHist[iBin]->GetXaxis()->GetXmax()-0.001));
    
      delete temp;
    }
  
    //Now, extract the upper and lower error bands.  Define the error
    //band as the shortest interval containing the central value as well
    //as 68% of the samples.
    TH1 *errHist = (TH1 *)fitHist->Clone("errHist_"+catName+"_"+fitLabel);
    errHist->Reset();
    errHist->SetFillColor(1);
    errHist->SetFillStyle(3654);
    errHist->SetLineColor(0);
    errHist->SetMarkerColor(0);
    errHist->SetMarkerStyle(1);
    errHist->SetMarkerSize(0.);

    std::cout << "\t Finished toy creation.  Now finding uncertainty band" << std::endl;

    for (int xBin=0; xBin<numBins; ++xBin) {

      int yBinMin = fluctHist[xBin]->FindBin(fitHist->GetBinContent(xBin+1));

      int yBinMax = yBinMin; //To start...

      double sum = fluctHist[xBin]->GetBinContent(yBinMin);

      while (sum < 0.68269*nSample && yBinMin > 0 && yBinMax <= fluctHist[xBin]->GetNbinsX()) {

	double contMin = fluctHist[xBin]->GetBinContent(yBinMin-1);
	double contMax = fluctHist[xBin]->GetBinContent(yBinMax+1);

	if (contMin > contMax) {
	  sum += contMin;        
	  --yBinMin;
	} else if (contMax > contMin) {
	  sum += contMax;
	  ++yBinMax;
	} else {
	  sum += contMax;
	  sum += contMin;
	  --yBinMin;
	  ++yBinMax;
	}

      }

      if (yBinMin == 0) std::cerr << "WARNING: error extends to underflow!" << std::endl;
      if (yBinMax == fluctHist[xBin]->GetNbinsX()+1) std::cerr << "WARNING: error extends to overflow!" << std::endl;

      double yMin = fluctHist[xBin]->GetBinCenter(yBinMin);
      double yMax = fluctHist[xBin]->GetBinCenter(yBinMax);
      double yAve = (yMin+yMax)/2;

      errHist->SetBinContent(xBin+1,yAve);
      errHist->SetBinError(xBin+1,yMax-yAve);
    }



    std::cout << "\t Finished uncertainty band creation.  Making plots" << std::endl;

    //Now calculate the ratio
    TH1 *ratio = (TH1 *)dataHist->Clone("ratio"+catName+"_"+fitLabel);
    ratio->SetMarkerStyle(1);

    TH1 *ratioErr = (TH1 *)ratio->Clone("ratioErr"+catName+"_"+fitLabel);
    ratioErr->Reset();
    ratioErr->SetFillColor(kGreen);
    ratioErr->SetLineColor(0);
    ratioErr->SetLineWidth(0);
    ratioErr->SetMarkerStyle(0);
    ratioErr->SetMarkerColor(0);
    ratioErr->SetMarkerSize(0.);

    for (int iBin = 1; iBin <= numBins; ++iBin) {

      double dataCentral = dataHist->GetBinContent(iBin);
      double mcCentral = fitHist->GetBinContent(iBin);
      double mcUp = errHist->GetBinContent(iBin)+errHist->GetBinError(iBin);
      double mcDown = errHist->GetBinContent(iBin)-errHist->GetBinError(iBin);

      double ratioUp = mcUp/mcCentral;
      double ratioDown = mcDown/mcCentral;
      double ratioAve = (ratioUp+ratioDown)/2;

      ratioErr->SetBinContent(iBin,ratioAve);
      ratioErr->SetBinError(iBin, ratioUp-ratioAve);

      double myratio = dataCentral/mcCentral;
      double myratio_err = sqrt( dataCentral ) / mcCentral;

      ratio->SetBinContent(iBin,myratio);
      ratio->SetBinError(iBin,myratio_err);

      if( (myratio>ratioMax) && ((myratio - myratio_err)<ratioMax) ){
	double minner = myratio - myratio_err;
	ratio->SetBinContent(iBin,ratioMax-0.0001);
	ratio->SetBinError(iBin,ratioMax-0.0001-minner);
      }
    }


    if( debug ) {
      for( int xBin = 1; xBin <= numBins; xBin++ ) {
	std::cout << " xBin = " << xBin << ",\t dataHist = " << dataHist->GetBinContent(xBin) << ",\t fitHist = " << fitHist->GetBinContent(xBin) << ",\t errHist = " << errHist->GetBinContent(xBin) << " +/- " << errHist->GetBinError(xBin) << ",\t ratio = " << ratio->GetBinContent(xBin) << ",\t ratioErr = " << ratioErr->GetBinContent(xBin) << " +/- " << ratioErr->GetBinError(xBin) << std::endl;
      }
    }

    //Hack to get it plotted with ratio plot
    TCanvas* myC = new TCanvas("myC", "myC", 600,700);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1.e-5;
    myC->Divide(1,2,small,small);
    const float padding=1e-5; const float ydivide=0.3;
    myC->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
    myC->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    myC->GetPad(1)->SetLeftMargin(.11);
    myC->GetPad(2)->SetLeftMargin(.11);
    myC->GetPad(1)->SetRightMargin(.05);
    myC->GetPad(2)->SetRightMargin(.05);
    myC->GetPad(1)->SetBottomMargin(.3);
    myC->GetPad(2)->SetBottomMargin(.3);
    myC->GetPad(1)->Modified();
    myC->GetPad(2)->Modified();
    myC->cd(1);
    gPad->SetBottomMargin(small);
    gPad->Modified();

    ratioErr->SetMinimum(ratioMin);
    ratioErr->SetMaximum(ratioMax);
    ratioErr->GetYaxis()->SetNdivisions(50000+404);
    ratioErr->GetYaxis()->SetLabelSize(0.1); //make y label bigger
    ratioErr->GetXaxis()->SetLabelSize(0.1); //make y label bigger
    ratioErr->GetXaxis()->SetTitleOffset(1.1);
    ratioErr->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle()); //make y label bigger
    ratioErr->GetXaxis()->SetLabelSize(0.12);
    ratioErr->GetXaxis()->SetLabelOffset(0.04);
    ratioErr->GetXaxis()->SetTitleSize(0.12);
    ratioErr->GetYaxis()->SetTitle("Data/MC");
    ratioErr->GetYaxis()->SetTitleSize(0.1);
    ratioErr->GetYaxis()->SetTitleOffset(.45);
    myC->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();


    ratioErr->SetTitle(";ANN output;Data/MC");


    dataHist->SetTitle(";;Events");
    dataHist->GetYaxis()->SetTitleSize(0.05);
    dataHist->GetYaxis()->SetTitleOffset(1.);

    double xmin = dataHist->GetBinLowEdge(1);
    double xmax = dataHist->GetBinLowEdge(numBins) + dataHist->GetBinWidth(numBins);


    h_ttH->Scale( fitHist->Integral() / h_ttH->Integral() );


    myC->cd(1);
    //myC->GetPad(1)->SetLogy(1);

    int maxBin_data = dataHist->GetMaximumBin();
    int maxBin_bkg  = errHist->GetMaximumBin();

    double max_data = dataHist->GetBinContent(maxBin_data) + dataHist->GetBinError(maxBin_data);
    double max_bkg  = errHist->GetBinContent(maxBin_bkg) + errHist->GetBinError(maxBin_bkg);
    double max_sig  = h_ttH->GetMaximum();

    double histMax = TMath::Max( max_data, max_bkg );
    histMax = TMath::Max( histMax, max_sig );

    dataHist->SetMaximum(1.1*histMax);
    fitHist->SetMaximum(1.1*histMax);
    errHist->SetMaximum(1.1*histMax);


    dataHist->Draw("e1");
    hs->Draw("histsame");
    errHist->Draw("e2same");
    //fitHist->Draw("HISTSAME");
    h_ttH->Draw("histsame");
    dataHist->Draw("e1SAME");


    CMSInfoLatex.Draw();
    SELECTIONInfoLatex.Draw();
    FitinfoLatex.Draw();


    myC->cd(2);
    ratioErr->Draw("e2");
    ratio->Draw("e1SAME");

    TLine* myLine;
    myLine = new TLine(xmin, 1.000, xmax, 1.000);
    myLine->Draw();

    myC->SaveAs("Images/dataToMC_"+era+"_plotWithErrors_"+fitLabel+"_"+catName+".png");
    myC->SaveAs("Images/dataToMC_"+era+"_plotWithErrors_"+fitLabel+"_"+catName+".pdf");


    std::cout << "\t Finished making (non-debug) plots" << std::endl;


    if( debug ){
      std::cout << "\t\t Beginning to debug " << std::endl;


      TCanvas *c2 = new TCanvas("c2","",900,800);


      TH2D* h_corr = new TH2D("h_corr","", NumPars, 0, NumPars, NumPars, 0, NumPars );
      for( int iPar=0; iPar<NumPars; iPar++ ){
	h_corr->GetXaxis()->SetBinLabel(iPar+1,nuisance_names[iPar]);
	h_corr->GetYaxis()->SetBinLabel(iPar+1,nuisance_names[iPar]);
	for( int jPar=0; jPar<NumPars; jPar++ ){
	  double correlation = fitFR->correlation( nuisance_names[iPar], nuisance_names[jPar] );
	  h_corr->SetBinContent(iPar+1,jPar+1,correlation);
	}
      }
      h_corr->Draw("colz");
      c2->Print("Images/correlations_nuisances_2D_"+fitLabel+"_"+catName+".png");


      TCanvas *c1 = new TCanvas("c1");
      for( int iPar=0; iPar<NumPars; iPar++ ){
	double mean = h_nuis[iPar]->GetMean();
	double rms  = h_nuis[iPar]->GetRMS();

	double meanErr = h_nuis[iPar]->GetMeanError();
	double rmsErr  = h_nuis[iPar]->GetRMSError();

	TString s_mean = Form("Mean = %.2f +/- %.2f, RMS = %.2f +/- %.2f",mean,meanErr,rms,rmsErr);
	TLatex MeanInfoLatex(0.15, 0.91, s_mean );
	MeanInfoLatex.SetNDC();
	MeanInfoLatex.SetTextFont(42);
	MeanInfoLatex.SetTextSize(0.035);

	h_nuis[iPar]->Draw();
	MeanInfoLatex.Draw();
	FitinfoLatex.Draw();
	c1->Print("Images/nuis_"+nuisance_names[iPar]+"_"+fitLabel+"_"+catName+".png");
      }


      TLine* l_nom = new TLine( 0,0,1,1 );
      TLine* l_ave = new TLine( 0,0,1,1 );
      TLine* l_ave_p1 = new TLine( 0,0,1,1 );
      TLine* l_ave_m1 = new TLine( 0,0,1,1 );

      l_nom->SetLineColor(kBlack);
      l_ave->SetLineColor(kGreen);
      l_ave_p1->SetLineColor(kRed);
      l_ave_m1->SetLineColor(kBlue);

      l_nom->SetLineWidth(2);
      l_ave->SetLineWidth(2);
      l_ave_p1->SetLineWidth(2);
      l_ave_m1->SetLineWidth(2);

      for( int iBin=0; iBin<numBins; iBin++ ){

	TString s_bin = Form("%d",iBin+1);
	TString proj_name = "projection_bin"+s_bin+"_"+fitLabel+"_"+catName;

	TH1D* h_proj = (TH1D*)fluctHist[iBin]->Clone(proj_name);

	int nProjBins = h_proj->GetNbinsX();
	std::cout << "  h_proj nbins = " << nProjBins << std::endl;

	double ave = errHist->GetBinContent(iBin+1);
	double aveErr = errHist->GetBinError(iBin+1);
	double nom = fitHist->GetBinContent(iBin+1);

	double xminProj = std::max( h_proj->GetXaxis()->GetXmin(), ave-3*aveErr );
	double xmaxProj = std::min( h_proj->GetXaxis()->GetXmax(), ave+3*aveErr );

	double ymaxProj = 1.05 * h_proj->GetMaximum();

	h_proj->GetXaxis()->SetRangeUser(xminProj,xmaxProj);
	h_proj->GetYaxis()->SetRangeUser(0,ymaxProj);

	h_proj->SetLineColor(kBlack);
	h_proj->Draw("hist");

	l_nom->DrawLine(nom,0,nom,ymaxProj);
	l_ave->DrawLine(ave,0,ave,ymaxProj);
	l_ave_p1->DrawLine(ave+aveErr,0,ave+aveErr,ymaxProj);
	l_ave_m1->DrawLine(ave-aveErr,0,ave-aveErr,ymaxProj);

	printf("  bin = %d,\t nom = %.2f,\t ave = %.2f,\t ave - err = %.2f,\t ave+err = %.2f \n",iBin+1,nom,ave,ave-aveErr,ave+aveErr );

	c1->Print("Images/projection_"+fitLabel+"_"+catName+"_bin"+s_bin+".png");
      }

      delete c1;
      delete c2;
      delete h_corr;

      std::cout << "\t\t Finished debugging " << std::endl;

    }

    //delete [] fluctHist;
    delete myLine;
    delete myC;

  }

  std::cout << "Completely done with plotWithErrrors.C" << std::endl;

}
