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

/*

Script to produce plots pre- and post-fit with uncertainty bands

to use, you need a datacard (e.g. datacard.dat), and a root file (e.g. file.root)

Do this once

combine -M MaxLikelihoodFit -m 125 --rMin -10 --rMax 10 --minos all datacard.dat

text2workspace.py -m 125 -D data_obs datacard.dat -b -o wsTest.root



Do this each time you want plots

root -b -q head.C plotWithErrors.C+'("file.root")'

*/

class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

void plotWithErrors( TString dataFileName = "", int nToys=500, bool blind=true, bool useLegend=true, bool debug=false, TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root" ) {


  TString imageDir = "Images";
  struct stat st;
  if( stat(imageDir.Data(),&st) != 0 )  mkdir(imageDir.Data(),0777);

  TString debugDir = imageDir + "/debug";
  if( debug && stat(debugDir.Data(),&st) != 0 )  mkdir(debugDir.Data(),0777);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  //Now we need a list of categories and a list of processes
  TList processes;
  processes.Add(new LabelInfo("ttbar","$t\\bar{t}+$LF"));
  processes.Add(new LabelInfo("ttbarPlusB","$t\\bar{t}+b$"));
  processes.Add(new LabelInfo("ttbarPlusBBbar","$t\\bar{t}+b\\bar{b}$"));
  processes.Add(new LabelInfo("ttbarPlusCCbar","$t\\bar{t}+c\\bar{c}$"));
  processes.Add(new LabelInfo("ttbarW","$t\\bar{t}W$"));
  processes.Add(new LabelInfo("ttbarZ","$t\\bar{t}Z$"));
  processes.Add(new LabelInfo("singlet","Single $t$"));
  processes.Add(new LabelInfo("wjets","$W+$jets"));
  processes.Add(new LabelInfo("zjets","$Z+$jets"));
  processes.Add(new LabelInfo("diboson","Diboson"));

  int NumBkgs = 10;

  int bin_ttjets_other = 0;
  int bin_ttjets_b     = 1;
  int bin_ttjets_bbbar = 2;
  int bin_ttjets_ccbar = 3;
  int bin_ttW     = 4;
  int bin_ttZ     = 5;
  int bin_singlet = 6;
  int bin_wjets   = 7;
  int bin_zjets   = 8;
  int bin_diboson = 9;

  int bin_ewk = 10;
  int bin_ttV = 11;

  int bin_ttH = 12;

  int bin_data = 13;

  Color_t color[22];
  color[0] = kBlack;
  //old
  color[bin_diboson] = kCyan;
  color[bin_zjets]   = kGreen+2;
  color[bin_wjets]   = kAzure+1;
  color[bin_ttW]     = kBlue-9;
  color[bin_ttZ]     = color[bin_ttW];

  //new
  color[bin_ttH]     = kBlue;

  color[bin_ttjets_bbbar]  = kRed+3;
  color[bin_ttjets_b]      = kRed-2;
  color[bin_ttjets_ccbar]  = kRed-7;
  color[bin_ttjets_other]  = kRed+1;

  color[bin_singlet] = kMagenta;
  color[bin_ttV] = kBlue-10;

  color[bin_ewk] = kAzure+2;



  std::vector<TString> label(22);
  label[bin_data] = "Data";

  label[bin_singlet] = "single t";
  label[bin_ttjets_bbbar]  = "t#bar{t} + b#bar{b}";
  label[bin_ttjets_b]  = "t#bar{t} + b";
  label[bin_ttjets_ccbar]  = "t#bar{t} + c#bar{c}";
  label[bin_ttjets_other]  = "t#bar{t} + lf";

  label[bin_ttH] = "t#bar{t}H(125)";

  label[bin_ewk] = "EWK";
  label[bin_ttV] = "t#bar{t} + W,Z";

  label[bin_zjets] = "Z + jets";
  label[bin_wjets] = "W + jets";
  label[bin_diboson] = "WW, WZ, ZZ";



  //These are the list of possible categories
  TList allCategories;
  allCategories.Add(new LabelInfo("ljets_jge6_t2","Lepton + #geq6 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("ljets_j4_t3","Lepton + 4 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("ljets_j5_t3","Lepton + 5 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("ljets_jge6_t3","Lepton + #geq6 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("ljets_j4_t4","Lepton + 4 jets + 4 b-tags"));
  allCategories.Add(new LabelInfo("ljets_j5_tge4","Lepton + 5 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("ljets_jge6_tge4","Lepton + #geq6 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("e3je2t","Dilepton + 3 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("ge4je2t","Dilepton + #geq4 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("ge3t","Dilepton + #geq3 b-tags"));

  // Get real tau category titles
  allCategories.Add(new LabelInfo("TTL_1b_1nb","Lep + #tau_{h}#tau_{h} + 4 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_1b_2nb","Lep + #tau_{h}#tau_{h} + 5 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_1b_3+nb","Lep + #tau_{h}#tau_{h} + #geq6 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_2b_0nb","Lep + #tau_{h}#tau_{h} + 4 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("TTL_2b_1nb","Lep + #tau_{h}#tau_{h} + 5 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("TTL_2b_2+nb","Lep + #tau_{h}#tau_{h} + #geq6 jets + 2 b-tags"));

  // Add SS because why not
  allCategories.Add(new LabelInfo("SS_ge4je1t","SS_ge4je1t"));
  allCategories.Add(new LabelInfo("SS_e3jge2t","SS_e3jge2t"));
  allCategories.Add(new LabelInfo("SS_ge4jge2t","SS_ge4jge2t"));


  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open( wsFileName );
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open( fitFileName );
  TFile *dataFile = TFile::Open(dataFileName);


  TList categories;
  TIter nextCat(&allCategories);
  LabelInfo *cat = 0;
  while (( cat = (LabelInfo *)nextCat())) {
    TString name = cat->name;
    if( w->pdf("pdf_bin"+name) ) categories.Add(cat);
  }

  int numCats = categories.GetEntries();//.GetSize();


  //These are the list of possible categories
  TList fitResults;
  fitResults.Add(new LabelInfo("nuisances_prefit_res","preFit"));
  fitResults.Add(new LabelInfo("fit_b","postFitB"));
  //// uncomment out below to see S+B fit
  //fitResults.Add(new LabelInfo("fit_s","postFitS"));


  TString era = "8TeV";

  std::string cmsinfo = "CMS Preliminary,  #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}";
  TLatex CMSInfoLatex(0.47, 0.91, cmsinfo.c_str());
  CMSInfoLatex.SetNDC();
  CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.035);

  //Suppress the printout from making plots
  //gErrorIgnoreLevel = 2000;

  //RooFit::RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling) ;
 
  //gSystem->Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so");



  double ratioMax = 2.3;
  double ratioMin = 0.0;

  TIter nextFit(&fitResults);
  LabelInfo *fitRes = 0;

  while ((fitRes = (LabelInfo *)nextFit())) {

    TString fitName = fitRes->name;
    TString fitLabel = fitRes->label;

    std::cout << " ===> Fit type = " << fitLabel << std::endl;

    RooFitResult *fitFR = (RooFitResult*)fitFile->Get( fitName );

    std::vector<TString> nuisance_names;
    nuisance_names.clear();

    std::vector<TString> nuisance_names_noMCstat;
    nuisance_names_noMCstat.clear();
    RooArgSet myArgs = fitFR->floatParsFinal();

    TIterator *nextArg(myArgs.createIterator());
    for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
      //RooRealVar *rrv = a;      
      if( debug ) std::cout << " name = " << rrv->GetName() << ":\t old value = " << rrv->getVal() << ",\t error = " << rrv->getError() << std::endl;
      TString nuis_name = rrv->GetName();
      nuisance_names.push_back( nuis_name );
      if( !nuis_name.Contains("ANNbin") ) nuisance_names_noMCstat.push_back( nuis_name );
    }


    std::cout << " Number of nuisances = " << int(nuisance_names.size()) << std::endl;
    std::cout << " Number of nuisances = " << int(nuisance_names_noMCstat.size()) << " (not including MCstats)" << std::endl;

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




    //Now get the fit central value
    w->saveSnapshot(fitLabel,RooArgSet(fitFR->floatParsFinal()),true);
    w->loadSnapshot(fitLabel);


    int numBins_most = 20;
    TH1D* fluctHist[numCats][numBins_most];
    RooAddition* fitTotal_cat[numCats];
    int numBins[numCats];

    TIter nextCat1(&categories);
    LabelInfo *category1 = 0;
    int iCat1=0;
    while ((category1 = (LabelInfo *)nextCat1())) {
      TString catLabel = category1->label;
      TString catName  = category1->name;

      //Roofit plots histogram vs "bin number;" turn this into a histogram over ANN output
      TH1 *dataHist = (TH1 *)dataFile->Get("data_obs_MVA_"+catName);
      TH1 *fitHist = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+fitLabel+"Clone_temp");
      numBins[iCat1] = fitHist->GetNbinsX();

      //Calculate the total in the fit...
      RooArgSet fitArgs;
      TIter nextProcess(&processes);
      LabelInfo *process = 0;
      //std::cout << "\t Filling individual sample components" << std::endl;
      while ((process = (LabelInfo *)nextProcess())) {
	TString procName = process->name;
	//Construct the name of the normalization variable.
	TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;
	RooAbsReal *fitN = w->function(normName);
	if (fitN) {
	  fitArgs.add(*fitN);
	}
      }
      fitTotal_cat[iCat1] = new RooAddition(Form("fitTotal_cat%d",iCat1),"fitTotal",fitArgs);

      //Now we need to fluctuate the nuisance parameters and get the fluctuated histograms
      //for( int iBin=0; iBin<numBins; iBin++ ) fluctHist[iBin] = new TH1D(Form("fluctHist_%d",iBin),"",300,0,5*fitHist->GetBinContent(iBin+1));
      for( int iBin=0; iBin<numBins[iCat1]; iBin++ ) fluctHist[iCat1][iBin] = new TH1D(Form("fluctHist_%d_%d_",iCat1,iBin)+catName+"_"+fitLabel,"",10000,0,4*fitHist->GetBinContent(iBin+1));
      iCat1++;
    }


    std::cout << " \t Begin making toys (nToys = " << nToys << ")" << std::endl;
    for (int iToy = 0; iToy < nToys; ++iToy) {
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

      TIter nextCat2(&categories);
      LabelInfo *category2 = 0;
      int iCat2=0;
      while ((category2 = (LabelInfo *)nextCat2())) {
	TString catLabel = category2->label;
	TString catName  = category2->name;

	TH1 *temp = w->pdf("pdf_bin"+catName)->createHistogram("CMS_th1x");
	temp->Scale(fitTotal_cat[iCat2]->getVal());

	for (int iBin=0; iBin<numBins[iCat2]; iBin++) fluctHist[iCat2][iBin]->Fill(std::min(temp->GetBinContent(iBin+1),fluctHist[iCat2][iBin]->GetXaxis()->GetXmax()-0.001));

	iCat2++;
	delete temp;
      }
    }
    std::cout << "\t Finished toy creation" << std::endl;


    if( debug ){
      std::cout << "\t\t Beginning to debug " << std::endl;

      TCanvas *c1 = new TCanvas("c1","",900,800);

      int NumPars_noMCstat = int(nuisance_names_noMCstat.size());
      TH2D* h_corr = new TH2D("h_corr","", NumPars_noMCstat, 0, NumPars_noMCstat, NumPars_noMCstat, 0, NumPars_noMCstat );
      for( int iPar=0; iPar<NumPars_noMCstat; iPar++ ){
	if( nuisance_names_noMCstat[iPar].Contains("ANNbin") ) continue;

	TString nuisanceName = nuisance_names_noMCstat[iPar];
	//// // If you want to change the names
	//nuisanceName = nuisanceName.ReplaceAll("CMS_ttH_","");
	//nuisanceName = nuisanceName.ReplaceAll("ttH_","");

	h_corr->GetXaxis()->SetBinLabel(iPar+1,nuisanceName);
	h_corr->GetYaxis()->SetBinLabel(iPar+1,nuisanceName);
	for( int jPar=0; jPar<NumPars_noMCstat; jPar++ ){
	  if( nuisance_names_noMCstat[jPar].Contains("ANNbin") ) continue;
	  double correlation = fitFR->correlation( nuisance_names_noMCstat[iPar], nuisance_names_noMCstat[jPar] );
	  h_corr->SetBinContent(iPar+1,jPar+1,correlation);
	}
      }
      h_corr->SetMaximum(1.);
      h_corr->SetMinimum(-1.);
      h_corr->LabelsOption("v","X");
      h_corr->Draw("colz");
      c1->GetPad(0)->SetLeftMargin(0.23);
      c1->GetPad(0)->SetBottomMargin(0.27);
      c1->GetPad(0)->SetTopMargin(0.05);

      c1->Print(debugDir+"/correlations_nuisances_2D_"+fitLabel+".png");


      TCanvas *c2 = new TCanvas("c2");
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
	c2->Print(debugDir+"/nuis_"+nuisance_names[iPar]+"_"+fitLabel+".png");
      }


      delete h_corr;
      delete c1;
      delete c2;
    }


    std::cout << " \t Begin making plots for each category " << std::endl;
    w->loadSnapshot(fitLabel);

    TIter nextCat3(&categories);
    LabelInfo *category3 = 0;
    int iCat3=0;
    while ((category3 = (LabelInfo *)nextCat3())) {
      TString catLabel = category3->label;
      TString catName  = category3->name;

      bool isTauCat = ( catName.Contains("TTL_") );
      std::cout << "\n\t\t category = " << catLabel << std::endl;
      
      TString selectioninfo = catLabel;//"Lepton + #geq6 jets + 2 tags";
      TLatex SELECTIONInfoLatex(0.12, 0.91, selectioninfo);
      SELECTIONInfoLatex.SetNDC();
      SELECTIONInfoLatex.SetTextFont(42);
      SELECTIONInfoLatex.SetTextSize(0.035);


      TH1 *dataHist = (TH1 *)dataFile->Get("data_obs_MVA_"+catName);
      dataHist->SetLineColor(1);
      dataHist->SetLineWidth(2);
      dataHist->SetMarkerStyle(20);
      dataHist->SetMarkerSize(0.75);
    



      TH1* h_bkg[NumBkgs];

      RooAbsPdf* nuisancePdf =  w->pdf("nuisancePdf");

      TIter nextProcess(&processes);
      LabelInfo *process = 0;
      int nProc = 0;
      //std::cout << "\t Filling individual sample components" << std::endl;
      while ((process = (LabelInfo *)nextProcess())) {

	TString procName = process->name;
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

	h_bkg[nProc] = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+procName+"_"+fitLabel+"Clone");
	h_bkg[nProc]->Reset();
	for (int iBin = 0; iBin < numBins[iCat3]; iBin++) {
	  h_bkg[nProc]->SetBinContent(iBin+1,bkgTempHist->GetBinContent(iBin+1));
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
      for (int iBin = 0; iBin < numBins[iCat3]; iBin++) {
	h_ttH->SetBinContent(iBin+1,sigTempHist->GetBinContent(iBin+1));
      }


      delete sigTempHist;

      h_ttH->SetLineColor(color[bin_ttH]);
      h_ttH->SetLineWidth(4);


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
      h_bkg[bin_ttjets_b]->SetFillColor(color[bin_ttjets_b]);

      h_bkg[bin_zjets]->SetFillColor(color[bin_zjets]);
      h_bkg[bin_wjets]->SetFillColor(color[bin_wjets]);
      h_bkg[bin_diboson]->SetFillColor(color[bin_diboson]);


      h_ewk_bkg->SetLineColor(color[bin_ewk]);
      h_ttbarV_bkg->SetLineColor(color[bin_ttV]);
      h_bkg[bin_singlet]->SetLineColor(color[bin_singlet]);
      h_bkg[bin_ttjets_other]->SetLineColor(color[bin_ttjets_other]);
      h_bkg[bin_ttjets_ccbar]->SetLineColor(color[bin_ttjets_ccbar]);
      h_bkg[bin_ttjets_bbbar]->SetLineColor(color[bin_ttjets_bbbar]);
      h_bkg[bin_ttjets_b]->SetLineColor(color[bin_ttjets_b]);

      h_bkg[bin_zjets]->SetLineColor(color[bin_zjets]);
      h_bkg[bin_wjets]->SetLineColor(color[bin_wjets]);
      h_bkg[bin_diboson]->SetLineColor(color[bin_diboson]);


      THStack *hs = new THStack("hs"+catName+"_"+fitLabel,"");

      if( isTauCat ){
	hs->Add(h_bkg[bin_diboson]);
	hs->Add(h_bkg[bin_ttjets_other]);
	hs->Add(h_bkg[bin_wjets]);
	hs->Add(h_bkg[bin_zjets]);
	hs->Add(h_bkg[bin_singlet]);
	hs->Add(h_ttbarV_bkg);
      }
      else {
	hs->Add(h_ewk_bkg);
	hs->Add(h_ttbarV_bkg);
	hs->Add(h_bkg[bin_singlet]);
	hs->Add(h_bkg[bin_ttjets_bbbar]);
	hs->Add(h_bkg[bin_ttjets_b]);
	hs->Add(h_bkg[bin_ttjets_ccbar]);
	hs->Add(h_bkg[bin_ttjets_other]);
      }

      //Make a normalized plot...
      TH1 *fitTempHist = w->pdf("pdf_bin"+catName+"_bonly")->createHistogram("CMS_th1x");
      fitTempHist->Scale(fitTotal_cat[iCat3]->getVal());

      //Roofit plots histogram vs "bin number;" turn this into a histogram over ANN output
      TH1 *fitHist = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+fitLabel+"Clone");
      fitHist->Reset();
      for (int iBin = 0; iBin < numBins[iCat3]; iBin++) {
	fitHist->SetBinContent(iBin+1,fitTempHist->GetBinContent(iBin+1));
      }
  

      fitHist->SetLineColor(kGreen+1);
      fitHist->SetLineWidth(2);
      fitHist->SetTitle(";;Events");
  
      //Just to be safe, ditch the temp histogram
      delete fitTempHist;
    

      //Now, extract the upper and lower error bands.  Define the error
      //band as the shortest interval containing the central value as well
      //as 68% of the samples.
      TH1 *errHist_1sig = (TH1 *)fitHist->Clone("errHist_1sig_"+catName+"_"+fitLabel);
      errHist_1sig->Reset();
      errHist_1sig->SetFillColor(1);
      errHist_1sig->SetFillStyle(3654);
      errHist_1sig->SetLineColor(1);
      errHist_1sig->SetMarkerColor(0);
      errHist_1sig->SetMarkerStyle(1);
      errHist_1sig->SetMarkerSize(0.);


      for (int xBin=0; xBin<numBins[iCat3]; ++xBin) {

	int yBinMin = fluctHist[iCat3][xBin]->FindBin(fitHist->GetBinContent(xBin+1));
	int yBinMax = yBinMin; //To start...

	double sum = fluctHist[iCat3][xBin]->GetBinContent(yBinMin);

	while (sum < 0.68269*nToys && yBinMin > 0 && yBinMax <= fluctHist[iCat3][xBin]->GetNbinsX()) {

	  double contMin = fluctHist[iCat3][xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist[iCat3][xBin]->GetBinContent(yBinMax+1);

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

	if (yBinMin == 0) std::cerr << "WARNING (1sig): error extends to underflow for bin " << xBin+1 << "!" << std::endl;
	if (yBinMax == fluctHist[iCat3][xBin]->GetNbinsX()+1) std::cerr << "WARNING (1sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist[iCat3][xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist[iCat3][xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_1sig->SetBinContent(xBin+1,yAve);
	errHist_1sig->SetBinError(xBin+1,yMax-yAve);
      }


      TH1 *errHist_2sig = (TH1 *)fitHist->Clone("errHist_2sig_"+catName+"_"+fitLabel);
      errHist_2sig->Reset();
      errHist_2sig->SetFillColor(1);
      errHist_2sig->SetFillStyle(3654);
      errHist_2sig->SetLineColor(0);
      errHist_2sig->SetMarkerColor(0);
      errHist_2sig->SetMarkerStyle(1);
      errHist_2sig->SetMarkerSize(0.);

      //std::cout << "\t Now finding 2 sigma uncertainty band" << std::endl;

      for (int xBin=0; xBin<numBins[iCat3]; ++xBin) {

	int yBinMin = fluctHist[iCat3][xBin]->FindBin(fitHist->GetBinContent(xBin+1));

	int yBinMax = yBinMin; //To start...

	double sum = fluctHist[iCat3][xBin]->GetBinContent(yBinMin);

	while (sum < 0.954499736*nToys && yBinMin > 0 && yBinMax <= fluctHist[iCat3][xBin]->GetNbinsX()) {

	  double contMin = fluctHist[iCat3][xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist[iCat3][xBin]->GetBinContent(yBinMax+1);

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

	if (yBinMin == 0) std::cerr << "WARNING (2sig): error extends to underflow for bin " << xBin+1 << "!" << std::endl;
	if (yBinMax == fluctHist[iCat3][xBin]->GetNbinsX()+1) std::cerr << "WARNING (2sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist[iCat3][xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist[iCat3][xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_2sig->SetBinContent(xBin+1,yAve);
	errHist_2sig->SetBinError(xBin+1,yMax-yAve);
      }



      //std::cout << "\t Finished uncertainty band creation.  Making plots" << std::endl;
      if( blind ){
	double threshold = 0.03;
	bool aboveThreshold = false;
	for (int iBin=0; iBin<numBins[iCat3]; iBin++) {
	  double bkgInBin = fitHist->GetBinContent(iBin+1);
	  double sigInBin = h_ttH->GetBinContent(iBin+1);
	  double SoverB = ( bkgInBin>0. ) ? sigInBin/bkgInBin : 0.;
	  aboveThreshold = ( aboveThreshold || (SoverB>threshold) );
	  if( aboveThreshold ) dataHist->SetBinContent(iBin+1,0.);
	}
      }

      //Now calculate the ratio
      TH1 *ratio = (TH1 *)dataHist->Clone("ratio"+catName+"_"+fitLabel);
      ratio->SetMarkerStyle(1);

      TH1 *ratioErr_1sig = (TH1 *)ratio->Clone("ratioErr_1sig"+catName+"_"+fitLabel);
      ratioErr_1sig->Reset();
      ratioErr_1sig->SetFillColor(kGreen);
      ratioErr_1sig->SetLineColor(0);
      ratioErr_1sig->SetLineWidth(0);
      ratioErr_1sig->SetMarkerStyle(0);
      ratioErr_1sig->SetMarkerColor(0);
      ratioErr_1sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numBins[iCat3]; iBin++) {

	double dataCentral = dataHist->GetBinContent(iBin);
	double mcCentral = fitHist->GetBinContent(iBin);
	double mcUp = errHist_1sig->GetBinContent(iBin)+errHist_1sig->GetBinError(iBin);
	double mcDown = errHist_1sig->GetBinContent(iBin)-errHist_1sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_1sig->SetBinContent(iBin,ratioAve);
	ratioErr_1sig->SetBinError(iBin, ratioUp-ratioAve);

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


      TH1 *ratioErr_2sig = (TH1 *)ratio->Clone("ratioErr_2sig"+catName+"_"+fitLabel);
      ratioErr_2sig->Reset();
      ratioErr_2sig->SetFillColor(kYellow);
      ratioErr_2sig->SetLineColor(0);
      ratioErr_2sig->SetLineWidth(0);
      ratioErr_2sig->SetMarkerStyle(0);
      ratioErr_2sig->SetMarkerColor(0);
      ratioErr_2sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numBins[iCat3]; iBin++) {

	//double dataCentral = dataHist->GetBinContent(iBin);
	double mcCentral = fitHist->GetBinContent(iBin);
	double mcUp = errHist_2sig->GetBinContent(iBin)+errHist_2sig->GetBinError(iBin);
	double mcDown = errHist_2sig->GetBinContent(iBin)-errHist_2sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_2sig->SetBinContent(iBin,ratioAve);
	ratioErr_2sig->SetBinError(iBin, ratioUp-ratioAve);
      }



      if( debug ) {
	for( int xBin = 1; xBin <= numBins[iCat3]; xBin++ ) {
	  std::cout << " xBin = " << xBin << ",\t dataHist = " << dataHist->GetBinContent(xBin) << ",\t fitHist = " << fitHist->GetBinContent(xBin) << ",\t errHist_1sig = " << errHist_1sig->GetBinContent(xBin) << " +/- " << errHist_1sig->GetBinError(xBin) << ",\t ratio = " << ratio->GetBinContent(xBin) << ",\t ratioErr_1sig = " << ratioErr_1sig->GetBinContent(xBin) << " +/- " << ratioErr_1sig->GetBinError(xBin) << std::endl;
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

      ratioErr_1sig->SetMinimum(ratioMin);
      ratioErr_1sig->SetMaximum(ratioMax);
      ratioErr_1sig->GetYaxis()->SetNdivisions(50000+404);
      ratioErr_1sig->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratioErr_1sig->GetXaxis()->SetLabelSize(0.1); //make y label bigger
      ratioErr_1sig->GetXaxis()->SetTitleOffset(1.1);
      ratioErr_1sig->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle()); //make y label bigger
      ratioErr_1sig->GetXaxis()->SetLabelSize(0.12);
      ratioErr_1sig->GetXaxis()->SetLabelOffset(0.04);
      ratioErr_1sig->GetXaxis()->SetTitleSize(0.12);
      ratioErr_1sig->GetYaxis()->SetTitle("Data/MC");
      ratioErr_1sig->GetYaxis()->SetTitleSize(0.1);
      ratioErr_1sig->GetYaxis()->SetTitleOffset(.45);
      myC->cd(2);
      gPad->SetTopMargin(small);
      gPad->SetTickx();
      gPad->Modified();

      
      TLegend *legend = new TLegend(0.14,0.75,0.94,0.89);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);

      legend->SetNColumns(3);

      double scale_ttH = ( fitLabel.EqualTo("postFitS") && !blind ) ? 1. : 30.; // fitHist->Integral() / h_ttH->Integral()


      if( isTauCat ){
	scale_ttH = ( fitLabel.EqualTo("postFitS") && !blind ) ? 1. : 100.;

	legend->AddEntry(h_bkg[bin_ttjets_other],Form("%s (%.1f)",label[bin_ttjets_other].Data(),h_bkg[bin_ttjets_other]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_diboson],Form("%s (%.1f)",label[bin_diboson].Data(),h_bkg[bin_diboson]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_wjets],Form("%s (%.1f)",label[bin_wjets].Data(),h_bkg[bin_wjets]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_zjets],Form("%s (%.1f)",label[bin_zjets].Data(),h_bkg[bin_zjets]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_singlet],Form("%s (%.1f)",label[bin_singlet].Data(),h_bkg[bin_singlet]->Integral()),"f");
	legend->AddEntry(h_ttbarV_bkg,Form("%s (%.1f)",label[bin_ttV].Data(),h_ttbarV_bkg->Integral()),"f");
      }
      else{
	legend->AddEntry(h_bkg[bin_ttjets_other],Form("%s (%.1f)",label[bin_ttjets_other].Data(),h_bkg[bin_ttjets_other]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_ttjets_ccbar],Form("%s (%.1f)",label[bin_ttjets_ccbar].Data(),h_bkg[bin_ttjets_ccbar]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_ttjets_b],Form("%s (%.1f)",label[bin_ttjets_b].Data(),h_bkg[bin_ttjets_b]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_ttjets_bbbar],Form("%s (%.1f)",label[bin_ttjets_bbbar].Data(),h_bkg[bin_ttjets_bbbar]->Integral()),"f");
	legend->AddEntry(h_bkg[bin_singlet],Form("%s (%.1f)",label[bin_singlet].Data(),h_bkg[bin_singlet]->Integral()),"f");

	legend->AddEntry(h_ttbarV_bkg,Form("%s (%.1f)",label[bin_ttV].Data(),h_ttbarV_bkg->Integral()),"f");
	legend->AddEntry(h_ewk_bkg,Form("%s (%.1f)",label[bin_ewk].Data(),h_ewk_bkg->Integral()),"f");
      }

      legend->AddEntry(dataHist,Form("%s (%.1f)",label[bin_data].Data(),dataHist->Integral()),"lpe");

      legend->AddEntry(errHist_1sig,Form("Sum MC (%.1f)",fitHist->Integral()),"f");

      if( fitLabel.EqualTo("postFitS") && !blind ) legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f)",h_ttH->Integral()),"l");
      else                                         legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f#times%.0f)",h_ttH->Integral(),scale_ttH),"l");




      //TLegend *legend_ratio = new TLegend(0.14,0.15,0.89,0.19);
      TLegend *legend_ratio = new TLegend(0.14,0.89,0.34,0.99);

      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetShadowColor(kWhite);
      legend_ratio->SetTextFont(42);
      legend_ratio->SetTextSize(0.055);

      if( !fitLabel.EqualTo("preFit") ){
	legend_ratio->SetNColumns(2);

	legend_ratio->AddEntry(ratioErr_1sig,"Fit #pm1#sigma","f");
	legend_ratio->AddEntry(ratioErr_2sig,"Fit #pm2#sigma","f");
      }
      else {
	legend_ratio->AddEntry(ratioErr_1sig,"Stat+Syst Uncertainty","f");
      }

      ratioErr_1sig->SetTitle(";MVA output;Data/MC");


      dataHist->SetTitle(";;Events");
      dataHist->GetYaxis()->SetTitleSize(0.05);
      dataHist->GetYaxis()->SetTitleOffset(1.);

      double xmin = dataHist->GetBinLowEdge(1);
      double xmax = dataHist->GetBinLowEdge(numBins[iCat3]) + dataHist->GetBinWidth(numBins[iCat3]);


      if( fitLabel.EqualTo("postFitS") && !blind ) hs->Add(h_ttH);
      else                                         h_ttH->Scale( scale_ttH );


      myC->cd(1);
      //myC->GetPad(1)->SetLogy(1);

      int maxBin_data = dataHist->GetMaximumBin();
      int maxBin_bkg  = errHist_1sig->GetMaximumBin();

      double max_data = dataHist->GetBinContent(maxBin_data) + dataHist->GetBinError(maxBin_data);
      double max_bkg  = errHist_1sig->GetBinContent(maxBin_bkg) + errHist_1sig->GetBinError(maxBin_bkg);
      double max_sig  = h_ttH->GetMaximum();

      double histMax = TMath::Max( max_data, max_bkg );
      histMax = TMath::Max( histMax, max_sig );

      if( useLegend ) histMax = 1.4 * histMax;
      else            histMax = 1.1 * histMax;

      dataHist->SetMaximum(histMax);
      fitHist->SetMaximum(histMax);
      errHist_1sig->SetMaximum(histMax);

      dataHist->SetMinimum(1.5E-1);
      dataHist->Draw("e1");
      hs->Draw("histsame");
      errHist_1sig->Draw("e2same");
      //fitHist->Draw("HISTSAME");
      if( !fitLabel.EqualTo("postFitS") || blind ) h_ttH->Draw("histsame");
      dataHist->Draw("e1SAME");

      if( useLegend ) legend->Draw();

      CMSInfoLatex.Draw();
      SELECTIONInfoLatex.Draw();
      FitinfoLatex.Draw();


      myC->cd(2);
      ratioErr_1sig->Draw("e2");
      if( !fitLabel.EqualTo("preFit") ){
	ratioErr_2sig->Draw("e2same");
	ratioErr_1sig->Draw("e2same");
      }
      ratio->Draw("e1SAME");

      myC->GetPad(1)->RedrawAxis();
      myC->GetPad(2)->RedrawAxis();


      TLine* myLine;
      myLine = new TLine(xmin, 1.000, xmax, 1.000);
      myLine->Draw();
      legend_ratio->Draw();

      myC->SaveAs(imageDir+"/dataToMC_"+era+"_plotWithErrors_"+fitLabel+"_"+catName+".png");
      myC->SaveAs(imageDir+"/dataToMC_"+era+"_plotWithErrors_"+fitLabel+"_"+catName+".pdf");

      iCat3++;
      delete myC;
      delete legend;
      delete legend_ratio;
    } // end loop over categories
  }// end loop over fitRes


  /*

  // debug info that is not currently implemented
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

	TH1D* h_proj = (TH1D*)fluctHist[iCat][iBin]->Clone(proj_name);

	int nProjBins = h_proj->GetNbinsX();
	std::cout << "  h_proj nbins = " << nProjBins << std::endl;

	double ave = errHist_1sig->GetBinContent(iBin+1);
	double aveErr = errHist_1sig->GetBinError(iBin+1);
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
  */
  std::cout << "Done!" << std::endl;

}
