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
#include "TVectorD.h"

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

/*

Script to produce yield tables per category

to use, you need a datacard (e.g. datacard.dat), and a root file (e.g. file.root)

combine -M MaxLikelihoodFit -m 125 --rMin -10 --rMax 10 --minos all datacard.dat -t -1

text2workspace.py -m 125 -D data_obs datacard.dat -b -o wsTest.root

root -b -q head.C printNorms.C'("file.root")'

*/


////
Double_t getPropagatedErrorMinos(const RooAbsReal &var, const RooFitResult &fitRes);
////

void printNorms(TString dataFileName = "", TString dataFileName7TeV = "", int fitType=0, TString prefix_ttH = "ttH125", TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root") {

  gStyle->SetOptStat(0);

  //Suppress the printout from making plots
  gErrorIgnoreLevel = 2000;
 
  //gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so");
  
  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open( wsFileName );
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open( fitFileName );
  RooFitResult *preFitFR = (RooFitResult*)fitFile->Get("nuisances_prefit_res");  
  if( fitType==1 )      preFitFR = (RooFitResult*)fitFile->Get("fit_b");
  else if( fitType==2 ) preFitFR = (RooFitResult*)fitFile->Get("fit_s");

  //Try making "snapshots" with preFit and postFit nuisance values
  if( fitType==1 )      w->saveSnapshot("postfitB",RooArgSet(preFitFR->floatParsFinal()),true);
  else if( fitType==2 ) w->saveSnapshot("postfitS",RooArgSet(preFitFR->floatParsFinal()),true);
  else                  w->saveSnapshot("prefit",RooArgSet(preFitFR->floatParsFinal()),true);

  //And if we want it, the data file
  TFile *dataFile = 0;
  if (dataFileName != "") {
    dataFile = TFile::Open(dataFileName);
  }
  TFile *dataFile7TeV = 0;
  if (dataFileName7TeV != "") {
    dataFile7TeV = TFile::Open(dataFileName7TeV);
  }


  //Now we need a list of categories and a list of processes
  TList processes;
  processes.Add(new LabelInfo("ttbar","$\\ttbar+$lf"));
  processes.Add(new LabelInfo("ttbarPlusB","$\\ttbar+$b"));
  processes.Add(new LabelInfo("ttbarPlusBBbar","$\\ttbar+\\bbbar$"));
  processes.Add(new LabelInfo("ttbarPlusCCbar","$\\ttbar+\\ccbar$"));
  processes.Add(new LabelInfo("ttbarW","$\\ttbar$W"));
  processes.Add(new LabelInfo("ttbarZ","$\\ttbar$Z"));
  processes.Add(new LabelInfo("singlet","Single t"));
  processes.Add(new LabelInfo("wjets","W+jets"));
  processes.Add(new LabelInfo("zjets","Z+jets"));
  processes.Add(new LabelInfo("diboson","Diboson"));

  int nRows = processes.GetEntries() + 1;
  if (dataFile) ++nRows;

  //These are the list of possible categories
  TList allCategories;
  // 7 TeV LJ + OSDIL
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_jge6_t2","Lepton + #geq6 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_j4_t3","Lepton + 4 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_j5_t3","Lepton + 5 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_jge6_t3","Lepton + #geq6 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_j4_t4","Lepton + 4 jets + 4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_j5_tge4","Lepton + 5 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ljets_jge6_tge4","Lepton + #geq6 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_e2je2t","Dilepton + 2 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("hbb_7TeV_ge3t","Dilepton + #geq3 b-tags"));

  // 8 TeV LJ + OSDIL
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_jge6_t2","Lepton + #geq6 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_j4_t3","Lepton + 4 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_j5_t3","Lepton + 5 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_jge6_t3","Lepton + #geq6 jets + 3 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_j4_t4","Lepton + 4 jets + 4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_j5_tge4","Lepton + 5 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ljets_jge6_tge4","Lepton + #geq6 jets + #geq4 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_e3je2t","Dilepton + 3 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ge4je2t","Dilepton + #geq4 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("hbb_8TeV_ge3t","Dilepton + #geq3 b-tags"));

  // 8 TeV TAU
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_1b_1nb","Lep + #tau_{h}#tau_{h} + 2 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_1b_2nb","Lep + #tau_{h}#tau_{h} + 3 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_1b_3+nb","Lep + #tau_{h}#tau_{h} + #geq4 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_2b_0nb","Lep + #tau_{h}#tau_{h} + 2 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_2b_1nb","Lep + #tau_{h}#tau_{h} + 3 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("htt_8TeV_TTL_2b_2+nb","Lep + #tau_{h}#tau_{h} + #geq4 jets + 2 b-tags"));

  //Now look in the postfit workspace to see which ones were actually used
  TList categories;

  TIter nextCat(&allCategories);
  LabelInfo *c = 0;

  while (( c = (LabelInfo *)nextCat())) categories.Add(c);

  const int NumSamples = 10+3;
  const int NumCats = 16 + 7 + 2;

  int NumCatsLJ = 7; // for 7 and 8 TeV
  int NumCatsDIL7TeV = 2;
  int NumCatsDIL8TeV = 3;

  double table_value[NumSamples][NumCats];
  double table_value_syst_err[NumSamples][NumCats];
  double table_value_stat_err[NumSamples][NumCats];
  

  std::vector<TString> table_labels;

  //Now loop over all categories, and within the category, loop over all processes:  

  TIter nextCategory(&categories);
  LabelInfo *category = 0;
  int iCat = 0;
  while ((category = (LabelInfo *)nextCategory())) {

    TString catName = category->name;
    TString catLabel = category->label;

    TString dataCatName = catName;
    if( dataCatName.Contains("ch1_") ) dataCatName.ReplaceAll("ch1_","");
    if( dataCatName.Contains("hbb_7TeV_") ) dataCatName.ReplaceAll("hbb_7TeV_","");
    if( dataCatName.Contains("hbb_8TeV_") ) dataCatName.ReplaceAll("hbb_8TeV_","");
    if( dataCatName.Contains("htt_8TeV_") ) dataCatName.ReplaceAll("htt_8TeV_","");

    TIter nextProcess(&processes);
    LabelInfo *process = 0;

    RooArgSet fitArgs;

    //bool firstProc = true;

    double ttH_val = -1;
    double ttH_err = -1;

    double total_sumw2 = 0.;

    TString ttHHistName = ( catName.BeginsWith("hbb_7TeV_") ) ? prefix_ttH + "_CFMlpANN_" + dataCatName : prefix_ttH + "_MVA_" + dataCatName;
    TH1 *ttHHist = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get(ttHHistName) : (TH1 *)dataFile->Get(ttHHistName);
    ttH_val = ttHHist->Integral();

    int nbins_ttH = ttHHist->GetNbinsX();
    double sumw2_ttH = 0.;
    for( int b=0; b<nbins_ttH; b++ ) sumw2_ttH += ttHHist->GetBinError(b+1)*ttHHist->GetBinError(b+1);

    TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";


    if( catName.BeginsWith("hbb_7TeV_") ){
      
      RooAbsReal *fitTTH = w->function(normName_ttH);
      if( fitTTH ){
	//w->loadSnapshot("prefit");        
	if( fitType==1 )      w->loadSnapshot("postfitB");
	else if( fitType==2 ) w->loadSnapshot("postfitS");
	else                  w->loadSnapshot("prefit");

	ttH_val = fitTTH->getVal();
	//ttH_err = fitTTH->getPropagatedError(*preFitFR);
	ttH_err = getPropagatedErrorMinos(*fitTTH,*preFitFR);
      }
    }
    else{

      RooArgSet fitSigArgs;

      // hbb
      TString normName_ttH_bb = "n_exp_final_bin"+catName+"_proc_ttH_hbb";
      RooAbsReal *fitTTH_bb = w->function(normName_ttH_bb);
      if( fitTTH_bb ) fitSigArgs.add( *fitTTH_bb );

      // hcc
      TString normName_ttH_cc = "n_exp_final_bin"+catName+"_proc_ttH_hcc";
      RooAbsReal *fitTTH_cc = w->function(normName_ttH_cc);
      if( fitTTH_cc ) fitSigArgs.add( *fitTTH_cc );

      // hgg
      TString normName_ttH_gg = "n_exp_final_bin"+catName+"_proc_ttH_hgg";
      RooAbsReal *fitTTH_gg = w->function(normName_ttH_gg);
      if( fitTTH_gg ) fitSigArgs.add( *fitTTH_gg );

      // hgluglu
      TString normName_ttH_gluglu = "n_exp_final_bin"+catName+"_proc_ttH_hgluglu";
      RooAbsReal *fitTTH_gluglu = w->function(normName_ttH_gluglu);
      if( fitTTH_gluglu ) fitSigArgs.add( *fitTTH_gluglu );

      // htt
      TString normName_ttH_tt = "n_exp_final_bin"+catName+"_proc_ttH_htt";
      RooAbsReal *fitTTH_tt = w->function(normName_ttH_tt);
      if( fitTTH_tt ) fitSigArgs.add( *fitTTH_tt );

      // hww
      TString normName_ttH_ww = "n_exp_final_bin"+catName+"_proc_ttH_hww";
      RooAbsReal *fitTTH_ww = w->function(normName_ttH_ww);
      if( fitTTH_ww ) fitSigArgs.add( *fitTTH_ww );

      // hzg
      TString normName_ttH_zg = "n_exp_final_bin"+catName+"_proc_ttH_hzg";
      RooAbsReal *fitTTH_zg = w->function(normName_ttH_zg);
      if( fitTTH_zg ) fitSigArgs.add( *fitTTH_zg );

      // hzz
      TString normName_ttH_zz = "n_exp_final_bin"+catName+"_proc_ttH_hzz";
      RooAbsReal *fitTTH_zz = w->function(normName_ttH_zz);
      if( fitTTH_zz ) fitSigArgs.add( *fitTTH_zz );




      //Sum up the category, pre and post fit.
      RooAddition fitSigTotal("fitSigTotal","fitSigTotal",fitSigArgs);

      //w->loadSnapshot("prefit");
      if( fitType==1 )      w->loadSnapshot("postfitB");
      else if( fitType==2 ) w->loadSnapshot("postfitS");
      else                  w->loadSnapshot("prefit");
      ttH_val = fitSigTotal.getVal();
      //     double preFitTotalErr = fitTotal.getPropagatedError(*preFitFR);
      ttH_err = getPropagatedErrorMinos(fitSigTotal,*preFitFR);
    }


    int iSample = 0;

    if( iCat==0 ) table_labels.push_back("$\\ttbar\\PH(125.6)$");
    table_value[iSample][iCat] = ttH_val;
    table_value_syst_err[iSample][iCat] = ttH_err;
    table_value_stat_err[iSample][iCat] = sqrt(sumw2_ttH);
    iSample++;

    while ((process = (LabelInfo *)nextProcess())) {

      TString procName = process->name;
      TString procLabel = process->label;

      //Construct the name of the normalization variable.
      TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;
      TString tempHistName = ( catName.BeginsWith("hbb_7TeV_") ) ? procName + "_CFMlpANN_"+dataCatName : procName + "_MVA_"+dataCatName;
      TH1 *procHist = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get(tempHistName) : (TH1 *)dataFile->Get(tempHistName);
    
      double sumw2 = 0.;

      if (procHist) {
        int nbins = procHist->GetNbinsX();
        for( int b=0; b<nbins; b++ ) sumw2 += procHist->GetBinError(b+1)*procHist->GetBinError(b+1);
      }

      //Extract normalization before fit
      double preFitNorm = -9e20, preFitErr = -9e10;

      RooAbsReal *fitN = w->function(normName);

      if (fitN) {
        fitArgs.add(*fitN);

        //Extract normalization and error on nomralization (before fit)
        //w->loadSnapshot("prefit");
        if( fitType==1 )      w->loadSnapshot("postfitB");
        else if( fitType==2 ) w->loadSnapshot("postfitS");
        else                  w->loadSnapshot("prefit");
        preFitNorm = fitN->getVal();
//         preFitErr = fitN->getPropagatedError(*preFitFR);
        preFitErr = getPropagatedErrorMinos(*fitN,*preFitFR);
      }

      total_sumw2 += sumw2;


      if (fitN) {

        if( iCat==0 ) table_labels.push_back(procLabel);
        table_value[iSample][iCat] = preFitNorm;
        table_value_syst_err[iSample][iCat] = preFitErr;
        table_value_stat_err[iSample][iCat] = sqrt(sumw2);
        iSample++;

      } else {

        if( iCat==0 ) table_labels.push_back(procLabel);
        table_value[iSample][iCat] = 0;
        table_value_syst_err[iSample][iCat] = 0;
        table_value_stat_err[iSample][iCat] = 0;
        iSample++;
      }

    }

    //Sum up the category, pre and post fit.
    RooAddition fitTotal("fitTotal","fitTotal",fitArgs);

    //w->loadSnapshot("prefit");
    if( fitType==1 )      w->loadSnapshot("postfitB");
    else if( fitType==2 ) w->loadSnapshot("postfitS");
    else                  w->loadSnapshot("prefit");
    double preFitTotalN = fitTotal.getVal();
//     double preFitTotalErr = fitTotal.getPropagatedError(*preFitFR);
    double preFitTotalErr = getPropagatedErrorMinos(fitTotal,*preFitFR);


    if( iCat==0 ) table_labels.push_back("Total bkg");
    table_value[iSample][iCat] = preFitTotalN;
    table_value_syst_err[iSample][iCat] = preFitTotalErr;
    table_value_stat_err[iSample][iCat] = sqrt(total_sumw2);
    iSample++;


    if (dataFile) {
      //Get the data normalization
      TString dataHistName = ( catName.BeginsWith("hbb_7TeV_") ) ? "data_obs_CFMlpANN_" + dataCatName : "data_obs_MVA_" + dataCatName;
      TH1 *dataHist = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get(dataHistName) : (TH1 *)dataFile->Get(dataHistName);
      double nData = dataHist->Integral();

      if( iCat==0 ) table_labels.push_back("Data");
      table_value[iSample][iCat] = nData;
      table_value_syst_err[iSample][iCat] = 0;
      table_value_stat_err[iSample][iCat] = 0;
      
    }

    ++iCat;

    
  } // while loop over categories



  std::cout << " ********************* " << std::endl;
  std::cout << "    LJ yields 7 TeV    " << std::endl;
  std::cout << " ********************* " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|c|c|c|c|c|} \\hline" << std::endl;
  std::cout << "& $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets \\\\" << std::endl;
  std::cout << "& 2 b-tags & 3 b-tags & 3 b-tags & 3 b-tags & 4 b-tags & $\\geq$4 b-tags & $\\geq$4 b-tags \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$\\ttbar$Z") || table_labels[iSample].EqualTo("Z+jets") ) continue;

    //These rows will never be filled for 7 TeV
    if( table_labels[iSample].EqualTo("$\\ttbar+$b") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$\\ttbar$W") )  std::cout << "$\\ttbar$+W/Z";
    else if( table_labels[iSample].EqualTo("W+jets") ) std::cout << "W/Z+jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=0; iCat2<NumCatsLJ; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$\\ttbar$W") || table_labels[iSample].EqualTo("W+jets") ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;

        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);
      }
      double err = sqrt( syst_err2 + stat_err2);

      if( err>100. ){
        val = 10 * double( int( val/10.+ 0.8 ) ); 
        err = 10 * double( int( err/10.+ 0.8 ) ); 
      }

      if( iSample<NumSamples-1 ){
        if( err>10. ) std::cout << Form(" & %.0f $\\pm$ %.0f", val, err);
        else          std::cout << Form(" & %.1f $\\pm$ %.1f", val, err);
      }
      else                       std::cout << Form(" & %.0f", val);
    }
    std::cout << " \\\\" << std::endl;
    if( table_labels[iSample].Contains("\\ttbar\\PH") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;




  std::cout << " ********************* " << std::endl;
  std::cout << "    DIL yields 7 TeV   " << std::endl;
  std::cout << " ********************* " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
  std::cout << "& 3 jets + 2 b-tags & $\\geq$4 jets + 2 b-tags & $\\geq$3 b-tags  \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$\\ttbar$Z") || table_labels[iSample].EqualTo("Z+jets") ) continue;

    //These rows will never be filled for 7 TeV
    if( table_labels[iSample].EqualTo("$\\ttbar+$b") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$\\ttbar$W") )  std::cout << "$\\ttbar$+W/Z";
    else if( table_labels[iSample].EqualTo("W+jets") ) std::cout << "W/Z+jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=NumCatsLJ; iCat2<NumCatsLJ+NumCatsDIL7TeV; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$\\ttbar$W") || table_labels[iSample].EqualTo("W+jets") ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;

        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);
      }
      double err = sqrt( syst_err2 + stat_err2);

      if( err>100. ){
        val = 10 * double( int( val/10.+ 0.8 ) ); 
        err = 10 * double( int( err/10.+ 0.8 ) ); 
      }

      if( iSample<NumSamples-1 ){
        if( err>10. ) std::cout << Form(" & %.0f $\\pm$ %.0f", val, err);
        else          std::cout << Form(" & %.1f $\\pm$ %.1f", val, err);
      }
      else                       std::cout << Form(" & %.0f", val);
    }
    std::cout << " \\\\" << std::endl;
    if( table_labels[iSample].Contains("\\ttbar\\PH") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;




  std::cout << " ********************* " << std::endl;
  std::cout << "    LJ yields 8 TeV    " << std::endl;
  std::cout << " ********************* " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|c|c|c|c|c|} \\hline" << std::endl;
  std::cout << "& $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets \\\\" << std::endl;
  std::cout << "& 2 b-tags & 3 b-tags & 3 b-tags & 3 b-tags & 4 b-tags & $\\geq$4 b-tags & $\\geq$4 b-tags \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$\\ttbar$Z") || table_labels[iSample].EqualTo("Z+jets") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$\\ttbar$W") )  std::cout << "$\\ttbar$+W/Z";
    else if( table_labels[iSample].EqualTo("W+jets") ) std::cout << "W/Z+jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=NumCatsLJ+NumCatsDIL7TeV; iCat2<NumCatsLJ+NumCatsDIL7TeV+NumCatsLJ; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$\\ttbar$W") || table_labels[iSample].EqualTo("W+jets") ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;

        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);
      }
      double err = sqrt( syst_err2 + stat_err2);

      if( err>100. ){
        val = 10 * double( int( val/10.+ 0.8 ) ); 
        err = 10 * double( int( err/10.+ 0.8 ) ); 
      }

      if( iSample<NumSamples-1 ){
        if( err>10. ) std::cout << Form(" & %.0f $\\pm$ %.0f", val, err);
        else          std::cout << Form(" & %.1f $\\pm$ %.1f", val, err);
      }
      else                       std::cout << Form(" & %.0f", val);
    }
    std::cout << " \\\\" << std::endl;
    if( table_labels[iSample].Contains("\\ttbar\\PH") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;




  std::cout << " ********************* " << std::endl;
  std::cout << "    DIL yields 8 TeV   " << std::endl;
  std::cout << " ********************* " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
  std::cout << "& 3 jets + 2 b-tags & $\\geq$4 jets + 2 b-tags & $\\geq$3 b-tags  \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$\\ttbar$Z") || table_labels[iSample].EqualTo("Z+jets") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$\\ttbar$W") )  std::cout << "$\\ttbar$+W/Z";
    else if( table_labels[iSample].EqualTo("W+jets") ) std::cout << "W/Z+jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=NumCatsLJ+NumCatsDIL7TeV+NumCatsLJ; iCat2<NumCatsLJ+NumCatsDIL7TeV+NumCatsLJ+NumCatsDIL8TeV; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$\\ttbar$W") || table_labels[iSample].EqualTo("W+jets") ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;

        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);
      }
      double err = sqrt( syst_err2 + stat_err2);

      if( err>100. ){
        val = 10 * double( int( val/10.+ 0.8 ) ); 
        err = 10 * double( int( err/10.+ 0.8 ) ); 
      }

      if( iSample<NumSamples-1 ){
        if( err>10. ) std::cout << Form(" & %.0f $\\pm$ %.0f", val, err);
        else          std::cout << Form(" & %.1f $\\pm$ %.1f", val, err);
      }
      else                       std::cout << Form(" & %.0f", val);
    }
    std::cout << " \\\\" << std::endl;
    if( table_labels[iSample].Contains("\\ttbar\\PH") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;



  std::cout << " ********************* " << std::endl;
  std::cout << "    Tau yields 8 TeV   " << std::endl;
  std::cout << " ********************* " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|c|c|c|c|} \\hline" << std::endl;
  std::cout << "& 2 jets & 3 jets & $\\geq$4 jets & 2 jets & 3 jets & $\\geq$4 jets \\\\" << std::endl;
  std::cout << "& 1 b-tag & 1 b-tag & 1 b-tag & 2 b-tags & 2 b-tags & 2 b-tags \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$\\ttbar$Z") || table_labels[iSample].EqualTo("Z+jets") ) continue;

    //These rows will never be filled for tau
    if( table_labels[iSample].EqualTo("$\\ttbar+\\bbbar$") || 
        table_labels[iSample].EqualTo("$\\ttbar+$b") || 
        table_labels[iSample].EqualTo("$\\ttbar+\\ccbar$")) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$\\ttbar$W") )  std::cout << "$\\ttbar$+W/Z";
    else if( table_labels[iSample].EqualTo("W+jets") ) std::cout << "W/Z+jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=NumCatsLJ+NumCatsDIL7TeV+NumCatsLJ+NumCatsDIL8TeV; iCat2<NumCats; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$\\ttbar$W") || table_labels[iSample].EqualTo("W+jets") ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;

        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);
      }
      double err = sqrt( syst_err2 + stat_err2);

      if( err>100. ){
        val = 10 * double( int( val/10.+ 0.8 ) ); 
        err = 10 * double( int( err/10.+ 0.8 ) ); 
      }

      if( iSample<NumSamples-1 ){
        if( err>10. ) std::cout << Form(" & %.0f $\\pm$ %.0f", val, err);
        else          std::cout << Form(" & %.1f $\\pm$ %.1f", val, err);
      }
      else                       std::cout << Form(" & %.0f", val);
    }
    std::cout << " \\\\" << std::endl;
    if( table_labels[iSample].Contains("\\ttbar\\PH") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;



}


Double_t getPropagatedErrorMinos(const RooAbsReal &var, const RooFitResult &fr) {

  // This code is stolen from RooAbsReal::getPropagatedError(), but
  // modified to use the average of the asymmetric errors (e.g. those
  // from MINOS) instead of the symmetric error (which ends up coming
  // from HESSE).
  //
  // Calculate error on [var] by propagated errors on
  // parameters with correlations as given by fit result The linearly
  // propagated error is calculated as follows:
  //
  // error(x) = F_a(x) * Corr(a,a') F_a'(x)
  //
  // where F_a(x) = [ f(x,a+da) - f(x,a-da) ] / 2, with f(x) as [var]
  // function and 'da' taken from the fit result (symmetrized MINOS
  // errors if available).  Corr(a,a') = the correlation matrix from the
  // fit result
  //


  // Clone var for internal use
  RooAbsReal* cloneFunc = (RooAbsReal*) var.cloneTree() ;
  RooArgSet* errorParams = cloneFunc->getObservables(fr.floatParsFinal()) ;
  RooArgSet* nset = cloneFunc->getParameters(*errorParams) ;
    
  // Make list of parameter instances of cloneFunc in order of error matrix
  RooArgList paramList ;
  const RooArgList& fpf = fr.floatParsFinal() ;
  vector<int> fpf_idx ;
  for (Int_t i=0 ; i<fpf.getSize() ; i++) {
    RooAbsArg* par = errorParams->find(fpf[i].GetName()) ;
    if (par) {
      paramList.add(*par) ;
      fpf_idx.push_back(i) ;
    }
  }

  vector<Double_t> plusVar, minusVar ;    
  
  // Create vector of plus,minus variations for each parameter  
  TMatrixDSym V(paramList.getSize()==fr.floatParsFinal().getSize()?
		fr.covarianceMatrix():
		fr.reducedCovarianceMatrix(paramList)) ;
  
  for (Int_t ivar=0 ; ivar<paramList.getSize() ; ivar++) {
    
    RooRealVar& rrv = (RooRealVar&)fpf[fpf_idx[ivar]] ;
    
    Double_t cenVal = rrv.getVal() ;

    const double thres = 0.00001;
    Double_t errValUp = fabs(rrv.getAsymErrorHi());
    Double_t errValDown = fabs(rrv.getAsymErrorLo());
    errValUp = ( errValUp<thres*rrv.getError() ) ? rrv.getError() : errValUp;
    errValDown = ( errValDown<thres*rrv.getError() ) ? rrv.getError() : errValDown;
    
    // Make Plus variation
    ((RooRealVar*)paramList.at(ivar))->setVal(cenVal+errValUp) ;
    plusVar.push_back(cloneFunc->getVal(nset)) ;
    
    // Make Minus variation
    ((RooRealVar*)paramList.at(ivar))->setVal(cenVal-errValDown) ;
    minusVar.push_back(cloneFunc->getVal(nset)) ;
    
    ((RooRealVar*)paramList.at(ivar))->setVal(cenVal) ;
  }
  
  TMatrixDSym C(paramList.getSize()) ;      
  for (int i=0 ; i<paramList.getSize() ; i++) {
    for (int j=i ; j<paramList.getSize() ; j++) {
      C(i,j) = V(i,j)/sqrt(V(i,i)*V(j,j)) ;
      C(j,i) = C(i,j) ;
    }
  }
  
  // Make vector of variations
  TVectorD F(plusVar.size()) ;
  for (unsigned int k=0 ; k<plusVar.size() ; k++) {
    F[k] = (plusVar[k]-minusVar[k])/2 ;
  }

  // Calculate error in linear approximation from variations and correlation coefficient
  Double_t sum = F*(C*F) ;

  delete cloneFunc ;
  delete errorParams ;
  delete nset ;

  return sqrt(sum) ;
}



