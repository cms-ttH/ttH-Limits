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
#include "TMatrix.h"
#include "TVector.h"
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
#include "TMath.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooRealVar.h"
#include "RooRandom.h"

/*

Script to produce plots pre- and post-fit with uncertainty bands

to use, you need a datacard (e.g. datacard.dat), and a root file (e.g. file.root)

Do this once

combine -M MaxLikelihoodFit -m 125 --rMin -10 --rMax 10 --minos all datacard.dat

text2workspace.py -m 125 -D data_obs datacard.dat -b -o wsTest.root



Do this each time you want plots

root -b -q head.C generateToys.C+'("file.root")'

*/

class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

//*****************************************************************************

TMatrix getUpperTriangularMatrix( RooFitResult *fitRes );

//*****************************************************************************


void generateToys( TString dataFileName = "", TString fitTypeName = "preFit", int nToys=500, int nJobs=1, int jobN=1, bool debug=false, TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root" ) {

  RooRandom::randomGenerator()->SetSeed(jobN);

  TString str_jobN = Form("%d",jobN);
  TString str_nJobs = Form("%d",nJobs);

  TString histoDir = "HistoFiles";
  struct stat st;
  if( stat(histoDir.Data(),&st) != 0 )  mkdir(histoDir.Data(),0777);

  TString histofilename = histoDir + "/generateToys_" + fitTypeName + "_histo_" + str_jobN + "_of_" + str_nJobs + ".root";
  TFile histofile(histofilename,"recreate");
  histofile.cd();


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
  allCategories.Add(new LabelInfo("TTL_1b_1nb","Lep + #tau_{h}#tau_{h} + 2 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_1b_2nb","Lep + #tau_{h}#tau_{h} + 3 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_1b_3+nb","Lep + #tau_{h}#tau_{h} + #geq4 jets + 1 b-tag"));
  allCategories.Add(new LabelInfo("TTL_2b_0nb","Lep + #tau_{h}#tau_{h} + 2 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("TTL_2b_1nb","Lep + #tau_{h}#tau_{h} + 3 jets + 2 b-tags"));
  allCategories.Add(new LabelInfo("TTL_2b_2+nb","Lep + #tau_{h}#tau_{h} + #geq4 jets + 2 b-tags"));

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
  if( fitTypeName.EqualTo("preFit") )        fitResults.Add(new LabelInfo("nuisances_prefit_res","preFit"));
  else if( fitTypeName.EqualTo("postFitB") ) fitResults.Add(new LabelInfo("fit_b","postFitB"));
  else if( fitTypeName.EqualTo("postFitS") ) fitResults.Add(new LabelInfo("fit_s","postFitS"));
  else {
    assert(0);
  }


  //Suppress the printout from making plots
  //gErrorIgnoreLevel = 2000;

  //RooFit::RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling) ;
 
  //gSystem->Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so");




  TIter nextFit(&fitResults);
  LabelInfo *fitRes = 0;

  while ((fitRes = (LabelInfo *)nextFit())) {

    TString fitName = fitRes->name;
    TString fitLabel = fitRes->label;

    std::cout << " ===> Fit type = " << fitLabel << std::endl;

    RooFitResult *fitFR = (RooFitResult*)fitFile->Get( fitName );

    std::vector<TString> nuisance_names;
    nuisance_names.clear();

    std::vector<double> nuisance_bestfit;
    nuisance_bestfit.clear();

    std::vector<TString> nuisance_names_noMCstat;
    nuisance_names_noMCstat.clear();
    RooArgSet myArgs = fitFR->floatParsFinal();

    TIterator *nextArg(myArgs.createIterator());
    for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
      TString nuis_name = rrv->GetName();
      nuisance_names.push_back( nuis_name );
      nuisance_bestfit.push_back( rrv->getVal() );
      if( !nuis_name.Contains("ANNbin") ) nuisance_names_noMCstat.push_back( nuis_name );
      if( debug ){
	printf(" name = %s: val = %.4f,\t getError = %.3f,\t getAsymErrorLo = %.3f,\t getAsymErrorHi = %.3f \n",
 	   nuis_name.Data(), rrv->getVal(), rrv->getError(), rrv->getAsymErrorLo(), rrv->getAsymErrorHi() );
      }
    }


    std::cout << " Number of nuisances = " << int(nuisance_names.size()) << std::endl;
    std::cout << " Number of nuisances = " << int(nuisance_names_noMCstat.size()) << " (not including MCstats)" << std::endl;

    int NumPars = int(nuisance_names.size());


    TH1D** h_nuis = new TH1D*[NumPars];
    for( int iPar=0; iPar<NumPars; iPar++ ){
      TString nuis_name = nuisance_names[iPar];
      h_nuis[iPar] = new TH1D( "h_nuis_"+nuis_name+"_"+fitLabel,"", 50, -4, 4 );
      h_nuis[iPar]->GetXaxis()->SetTitle(nuis_name);
    }





    //Now get the fit central value
    w->saveSnapshot(fitLabel,RooArgSet(fitFR->floatParsFinal()),kTRUE);
    w->loadSnapshot(fitLabel);


    int numBins_most = 20;
    TH1D* fluctHist[numCats][numBins_most];
    RooAddition** fitTotal_cat = new RooAddition*[numCats];
    int numBins[numCats];

    double useSoverB[numCats][numBins_most];

    double minLogSoverB = -3.5, maxLogSoverB = -1.0;
    int numBinsSoverB = 20;
    TH1D* fluctHist_SoverB[numBinsSoverB];

    TH1D* h_data_SoverB = new TH1D(Form("h_data_SoverB_%s",fitLabel.Data()),";log_{10}(S/B)",numBinsSoverB, minLogSoverB, maxLogSoverB);
    TH1D* h_bkg_SoverB = new TH1D(Form("h_bkg_SoverB_%s",fitLabel.Data()),";log_{10}(S/B)",numBinsSoverB, minLogSoverB, maxLogSoverB);
    TH1D* h_sig_SoverB = new TH1D(Form("h_sig_SoverB_%s",fitLabel.Data()),";log_{10}(S/B)",numBinsSoverB, minLogSoverB, maxLogSoverB);


    TH1D* fluctHist_category_yield[numCats];
    TH1D* h_category_yield_data = new TH1D(Form("h_category_yield_data_%s",fitLabel.Data()),";category",numCats,0,numCats);
    TH1D* h_category_yield_bkg  = new TH1D(Form("h_category_yield_bkg_%s",fitLabel.Data()),";category",numCats,0,numCats);
    TH1D* h_category_yield_sig  = new TH1D(Form("h_category_yield_sig_%s",fitLabel.Data()),";category",numCats,0,numCats);



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


      //Make a normalized plot...
      TH1 *fitTempHist = w->pdf("pdf_bin"+catName+"_bonly")->createHistogram("CMS_th1x");
      fitTempHist->Scale(fitTotal_cat[iCat1]->getVal());
      for (int iBin = 0; iBin < numBins[iCat1]; iBin++) {
	if( fitHist->GetBinContent(iBin+1)<1 ) fitHist->SetBinContent(iBin+1,fitTempHist->GetBinContent(iBin+1));
      }


      TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";
      RooAbsReal *fitTTH = w->function(normName_ttH);


      TH1* sigTempHist = w->pdf("shapeSig_"+catName+"_ttH_morph")->createHistogram("CMS_th1x");
      sigTempHist->Scale(fitTTH->getVal());

      for (int iBin = 0; iBin < numBins[iCat1]; iBin++) {
	double signal = sigTempHist->GetBinContent(iBin+1);
	double background = fitTempHist->GetBinContent(iBin+1);

	double logSoverB = ( background>0. ) ? TMath::Log10(signal/background) : -99.;
	useSoverB[iCat1][iBin] = logSoverB;
	double fillLogSoverB = TMath::Max( logSoverB, minLogSoverB );
	fillLogSoverB = TMath::Min( fillLogSoverB, maxLogSoverB );
	h_data_SoverB->Fill(fillLogSoverB,dataHist->GetBinContent(iBin+1));
	h_bkg_SoverB->Fill(fillLogSoverB,fitTempHist->GetBinContent(iBin+1));
	h_sig_SoverB->Fill(fillLogSoverB,sigTempHist->GetBinContent(iBin+1));

	if( signal/background>0.1 ) printf("%s\t%s\t%d\t S/B = %.2f\t S = %.1f\t B = %.1f\n",fitLabel.Data(),catName.Data(),iBin,signal/background,signal,background);
      }


      double integral_data = dataHist->Integral();
      double integral_bkg  = fitTempHist->Integral();
      double integral_sig  = sigTempHist->Integral();

      h_category_yield_data->Fill(iCat1,integral_data);
      h_category_yield_bkg->Fill(iCat1,integral_bkg);
      h_category_yield_sig->Fill(iCat1,integral_sig);

      

      //Now we need to fluctuate the nuisance parameters and get the fluctuated histograms
      for( int iBin=0; iBin<numBins[iCat1]; iBin++ ) fluctHist[iCat1][iBin] = new TH1D(Form("fluctHist_%s_%s_%d",fitLabel.Data(),catName.Data(),iBin),"",10000,0,5*fitHist->GetBinContent(iBin+1));

      iCat1++;

      delete fitTempHist;
      delete sigTempHist;
    }

    for( int iBin=0; iBin<numBinsSoverB; iBin++ ) fluctHist_SoverB[iBin] = new TH1D(Form("fluctHist_SoverB_%s_%d",fitLabel.Data(),iBin),"",100000,0,4*h_bkg_SoverB->GetBinContent(iBin+1));

    for( int iBin=0; iBin<numCats; iBin++ ) fluctHist_category_yield[iBin] = new TH1D(Form("fluctHist_category_yield_%s_%d",fitLabel.Data(),iBin),"",100000,0,4*h_category_yield_bkg->GetBinContent(iBin+1));


    TMatrix myL = getUpperTriangularMatrix( fitFR );
    TMatrix* _Lt = new TMatrix(TMatrix::kTransposed,myL);

    std::cout << " \t Begin making toys (nToys = " << nToys << ")" << std::endl;
    for (int iToy = 0; iToy < nToys; ++iToy) {
      //RooArgSet myArgs_temp = fitFR->randomizePars();
      RooArgSet myArgs_temp = fitFR->floatParsFinal();

      // create a vector of unit Gaussian variables
      TVector g(NumPars);
      for( int k=0; k<NumPars; k++ ) g(k)= RooRandom::gaussian();
      // multiply this vector by Lt to introduce the appropriate correlations
      g *= (*_Lt);

      if( iToy%10==0 ) std::cout << "\t\t iToy " << iToy << std::endl;
      int iPar=0;
      TIterator *nextArg_temp(myArgs_temp.createIterator());
      for (TObject *a = nextArg_temp->Next(); a != 0; a = nextArg_temp->Next()) {
	RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
	TString nuis_name = rrv->GetName();
	double oldValue = nuisance_bestfit[iPar];
	// add the mean value offsets and store the results
	rrv->setVal(oldValue + g(iPar));
	///// For debugging
	//double newValue = rrv->getVal();
	//printf(" name = %s: old val = %.4f,\t new val = %.4f,\t getError = %.3f,\t getAsymErrorLo = %.3f,\t getAsymErrorHi = %.3f \n", nuis_name.Data(), oldValue, newValue, rrv->getError(), rrv->getAsymErrorLo(), rrv->getAsymErrorHi() );
	if( !nuis_name.EqualTo(nuisance_names[iPar]) ) assert(0);
	h_nuis[iPar]->Fill(rrv->getVal());
	iPar++;
      }

      w->saveSnapshot("fluctuation"+fitLabel,RooArgSet(myArgs_temp),kTRUE);
      w->loadSnapshot("fluctuation"+fitLabel);

      TH1D* temp_SoverB = (TH1D*)h_data_SoverB->Clone(Form("temp_SoverB_%s_%d",fitLabel.Data(),iToy));
      temp_SoverB->Reset();

      TIter nextCat2(&categories);
      LabelInfo *category2 = 0;
      int iCat2=0;
      while ((category2 = (LabelInfo *)nextCat2())) {
	TString catLabel = category2->label;
	TString catName  = category2->name;

	TH1 *temp = w->pdf("pdf_bin"+catName)->createHistogram("CMS_th1x");
	temp->Scale(fitTotal_cat[iCat2]->getVal());

	fluctHist_category_yield[iCat2]->Fill(std::min(temp->Integral(),fluctHist_category_yield[iCat2]->GetXaxis()->GetXmax()-0.001));

	for (int iBin=0; iBin<numBins[iCat2]; iBin++) fluctHist[iCat2][iBin]->Fill(std::min(temp->GetBinContent(iBin+1),fluctHist[iCat2][iBin]->GetXaxis()->GetXmax()-0.001));

	for (int iBin=0; iBin<numBins[iCat2]; iBin++){
	  double logSoverB = useSoverB[iCat2][iBin];
	  double fillLogSoverB = TMath::Max( logSoverB, minLogSoverB );
	  fillLogSoverB = TMath::Min( fillLogSoverB, maxLogSoverB );
	  temp_SoverB->Fill(fillLogSoverB,temp->GetBinContent(iBin+1));
	  //if( temp_SoverB->FindBin(fillLogSoverB)==numBinsSoverB ) printf("%d\t%s\t%s\t%d\t log(S/B) = %.2f,\t content = %.1f,\t total = %.1f \n",iToy,fitLabel.Data(),catName.Data(),iBin,logSoverB,temp->GetBinContent(iBin+1),temp_SoverB->GetBinContent(temp_SoverB->FindBin(fillLogSoverB)));
	}

	iCat2++;
	delete temp;
      }

      for (int iBin=0; iBin<numBinsSoverB; iBin++) fluctHist_SoverB[iBin]->Fill(std::min(temp_SoverB->GetBinContent(iBin+1),fluctHist_SoverB[iBin]->GetXaxis()->GetXmax()-0.001));

    }
    std::cout << "\t Finished toy creation" << std::endl;


    histofile.cd();
    TIter nextCat5(&categories);
    LabelInfo *category5 = 0;
    int iCat5=0;
    while ((category5 = (LabelInfo *)nextCat5())) {
      TString catName  = category5->name;
      for( int iBin=0; iBin<numBins[iCat5]; iBin++ ) fluctHist[iCat5][iBin]->Write(Form("fluctHist_%s_%s_%d",fitLabel.Data(),catName.Data(),iBin));
      iCat5++;
    }
    for( int iPar=0; iPar<NumPars; iPar++ ) h_nuis[iPar]->Write("h_nuis_"+nuisance_names[iPar]+"_"+fitLabel);

    for( int iBin=0; iBin<numBinsSoverB; iBin++ ) fluctHist_SoverB[iBin]->Write(Form("fluctHist_SoverB_%s_%d",fitLabel.Data(),iBin));

    if( jobN!=1 ){
      h_data_SoverB->Reset();
      h_bkg_SoverB->Reset();
      h_sig_SoverB->Reset();

      h_category_yield_data->Reset();
      h_category_yield_bkg->Reset();
      h_category_yield_sig->Reset();
    }

    h_data_SoverB->Write(Form("data_SoverB_%s",fitLabel.Data()));
    h_bkg_SoverB->Write(Form("bkg_SoverB_%s",fitLabel.Data()));
    h_sig_SoverB->Write(Form("sig_SoverB_%s",fitLabel.Data()));




    for( int iBin=0; iBin<numCats; iBin++ ) fluctHist_category_yield[iBin]->Write(Form("fluctHist_category_yield_%s_%d",fitLabel.Data(),iBin));
    h_category_yield_data->Write(Form("data_category_yield_%s",fitLabel.Data()));
    h_category_yield_bkg->Write(Form("bkg_category_yield_%s",fitLabel.Data()));
    h_category_yield_sig->Write(Form("sig_category_yield_%s",fitLabel.Data()));



    TH1D* h_numToys = new TH1D("h_numToys","",3,0,3);
    if( fitTypeName.EqualTo("preFit") )        h_numToys->SetBinContent(1,nToys);
    else if( fitTypeName.EqualTo("postFitB") ) h_numToys->SetBinContent(2,nToys);
    else if( fitTypeName.EqualTo("postFitS") ) h_numToys->SetBinContent(3,nToys);
    h_numToys->Write("h_numToys");


    // delete quantities defined for each fit result
    for( int iPar=0; iPar<NumPars; iPar++ ) delete h_nuis[iPar];
    delete h_nuis;

    for( int iCat=0; iCat<numCats; iCat++ ) delete fitTotal_cat[iCat];
    delete fitTotal_cat;

  }// end loop over fitRes


  
  histofile.Write();
  histofile.Close();

  std::cout << "Done!" << std::endl;

}


TMatrix getUpperTriangularMatrix( RooFitResult *fitRes ){

  RooArgSet myArgs = fitRes->floatParsFinal();

  std::vector<TString> nuisance_names;
  nuisance_names.clear();

  std::vector<double> nuisance_errors;
  nuisance_errors.clear();

  TIterator *nextArg(myArgs.createIterator());
  for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      

    double thres = 0.00001;

    double errHi = fabs( rrv->getAsymErrorHi() );
    double errLo = fabs( rrv->getAsymErrorLo() );
    errHi = ( errHi<thres*rrv->getError() ) ? rrv->getError() : errHi;
    errLo = ( errLo<thres*rrv->getError() ) ? rrv->getError() : errLo;

    double err = 0.5 * ( errHi + errLo );

    TString nuis_name = rrv->GetName();
    nuisance_names.push_back( nuis_name );
    nuisance_errors.push_back( err );
//     printf(" name = %s: val = %.4f,\t err = %.4f,\t getError = %.3f,\t getAsymErrorLo = %.3f,\t getAsymErrorHi = %.3f \n",
// 	   nuis_name.Data(), rrv->getVal(), err, rrv->getError(), rrv->getAsymErrorLo(), rrv->getAsymErrorHi() );
  }

  int NumPars = int( nuisance_names.size() );

  TMatrix Cov(NumPars,NumPars);
  for( int iPar=0; iPar<NumPars; iPar++ ){
    for( int jPar=0; jPar<NumPars; jPar++ ){

      double correlation = fitRes->correlation( nuisance_names[iPar], nuisance_names[jPar] );
      double covariance  = nuisance_errors[iPar] * nuisance_errors[jPar] * correlation;

      Cov(iPar,jPar) = covariance;
      Cov(jPar,iPar) = covariance;
    }
  }


  TMatrix L(NumPars,NumPars);
  for( int iPar= 0; iPar<NumPars; iPar++ ){
    // calculate the diagonal term first
    L(iPar,iPar)= Cov(iPar,iPar);
    for( int kPar=0; kPar<iPar; kPar++) {
      Double_t tmp= L(kPar,iPar);
      L(iPar,iPar) -= tmp*tmp;
    }
    L(iPar,iPar)= sqrt(L(iPar,iPar));
    // then the off-diagonal terms
    for( int jPar=iPar+1; jPar<NumPars; jPar++ ){
      L(iPar,jPar)= Cov(iPar,jPar);
      for(Int_t kPar=0; kPar<iPar; kPar++) {
	L(iPar,jPar)-= L(kPar,iPar)*L(kPar,jPar);
      }
      L(iPar,jPar)/= L(iPar,iPar);
    }
  }

  return L;
}

