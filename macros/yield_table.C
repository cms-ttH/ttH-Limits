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


void yield_table(TString fileName = "result_combined.root") {

  TFile *file_input = TFile::Open(fileName);

  int numChannels = 9;
  TString channelNames[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    "ljets_jge6_t2", 		//	 "j6_2t",
    "ljets_j4_t3", 		//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    "ljets_j5_t3", 		//	 "j5_3t",
    "ljets_jge6_t3", 		//	 "j6_3t",
    "ljets_j4_t4", 		//	 "j4_3t",
    "ljets_j5_tge4", 		//	 "j5_3t",
    "ljets_jge6_tge4",    	//	 "j6_3t",
    "e2je2t", 			// dilep 2 jets 2 tags
    "ge3t"	       		// dilep >= 3tags
  };

  int numSamples = 9;
  TString sampleNames[] = {"ttbar","ttbarPlusBBbar","ttbarPlusCCbar","ttbarW","ttbarZ","singlet","wjets","zjets","diboson","ttH125"};
  TString sampleNames_latex[] = {"$t\\bar{t}+$lf","$t\\bar{t}+b\\bar{b}$","$t\\bar{t}+c\\bar{c}$","$t\\bar{t}W$","$t\\bar{t}Z$","Single $t$","$W+$jets","$Z+$jets","Diboson"};

  int numSyst = 11;
  TString systNames[] = {"CMS_fake_bRate","CMS_eff_bRate","CMS_ttH_PUcorr","Q2scale_ttH_ttbar0p","Q2scale_ttH_ttbar_bb","Q2scale_ttH_ttbar_cc","CMS_scale_j","_eff_bShape","_fake_bShape","Q2scale_ttH_ttbar1p","Q2scale_ttH_ttbar2p"};



  int NumRateSys = 12;//12;


  double lumi_unc = 1.044;
  double eff_lep_unc = 1.04;
  //      { "ttbar","ttbarPlusBBbar","ttbarPlusCCbar", "ttbarW", "ttbarZ", "singlet", "wjets", "zjets", "diboson", "ttH120"};
  double syst_vect[12][10] = {
    { lumi_unc,  lumi_unc, lumi_unc, lumi_unc, lumi_unc, lumi_unc, lumi_unc, lumi_unc, lumi_unc, lumi_unc },// "lumi"
    { 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.125 },                                                          // "QCDscale_ttH"
    { 1.12, 1.12, 1.12, 1.15, 1.15, 1.02, 0., 0., 0., 0. },                                                 // "QCDscale_ttbar"
    { 1.09, 1.09, 1.09, 0., 1.09, 0., 0., 0., 0., 1.08 },                                                   // "pdf_gg"
    { 0., 0., 0., 1.07, 0., 0., 1.048,	1.042, 0.,  0. },                                                   // "pdf_qqbar"
    { 0., 0., 0., 0., 0., 1.046, 0., 0., 0., 0. },                                                          // "pdf_qg"
    { 0., 0., 0., 0., 0., 0., 1.013, 1.012, 0., 0. },                                                       // "QCDscale_V
    { 0., 0., 0., 0., 0., 0., 0., 0., 1.035, 0. },                                                          // "QCDscale_VV"
    { 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01 },                                         // "CMS_ttH_pu"
    { 1.015, 1.015, 1.015, 1.015, 1.015, 1.015, 1.015, 1.015, 1.015, 1.015 },                               // "CMS_res_j"
    { eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc, eff_lep_unc }, // "CMS_ttH_eff_lep"
    { 0., 1.5, 0., 0., 0., 0., 0., 0., 0., 0. }                                                             // "CMS_ttH_QCDscale_ttbb"
  };




  const int NumSamples = 9+3;
  const int NumCats = 9;

  double table_value[NumSamples][NumCats];
  double table_value_syst_err[NumSamples][NumCats];
  double table_value_stat_err[NumSamples][NumCats];
  std::vector<TString> table_labels;



  for( int iChan = 0; iChan < numChannels; ++iChan ){


    TString ttH_histName = "ttH125";
    ttH_histName += "_CFMlpANN_";
    ttH_histName += channelNames[iChan];

    TH1D* hist_ttH = (TH1D*)file_input->Get(ttH_histName);
    double yield_ttH = hist_ttH->Integral();

    int nbins_ttH = hist_ttH->GetNbinsX();
    double sumw2_ttH = 0.;
    for( int b=0; b<nbins_ttH; b++ ) sumw2_ttH += hist_ttH->GetBinError(b+1)*hist_ttH->GetBinError(b+1);

    int iSamp = 0;

    if( iChan==0 ) table_labels.push_back("$t\\bar{t}H(125)$");
    table_value[iSamp][iChan] = yield_ttH;

    double ttH_syst_up_2 = 0;
    double ttH_syst_down_2 = 0;
    for (int iSyst = 0; iSyst < numSyst; ++iSyst) {
      TString sysName = ( systNames[iSyst].Contains("bShape") ) ? channelNames[iChan] + systNames[iSyst] : systNames[iSyst];

      if( sysName.Contains("Q2scale_ttH_ttbar") ) continue;
      
      TString ttH_histName_up =   ttH_histName + "_" + sysName + "Up";
      TString ttH_histName_down = ttH_histName + "_" + sysName + "Down";

      TH1D* hist_up   = (TH1D*)file_input->Get(ttH_histName_up);
      TH1D* hist_down = (TH1D*)file_input->Get(ttH_histName_down);

      double syst_up   = hist_up->Integral() - yield_ttH;
      double syst_down = hist_down->Integral() - yield_ttH;

      double syst_up_2 = syst_up*syst_up;
      double syst_down_2 = syst_down*syst_down;

      ttH_syst_up_2 += syst_up_2;
      ttH_syst_down_2 += syst_down_2;
    }

    for (int iSyst = 0; iSyst < NumRateSys; ++iSyst) {
      double mysyst = ( syst_vect[iSyst][9]>1. ) ? syst_vect[iSyst][9] - 1. : 0.;
      mysyst *= yield_ttH;
      ttH_syst_up_2 += mysyst*mysyst;
      ttH_syst_down_2 += mysyst*mysyst;
    }
    

    table_value_syst_err[iSamp][iChan] = 0.5* ( sqrt(ttH_syst_up_2) + sqrt(ttH_syst_down_2) );//sqrt( ttH_syst_2 );
    table_value_stat_err[iSamp][iChan] = sqrt(sumw2_ttH);
    iSamp++;


    double sum_bkg = 0.;
    double sumw2_bkg = 0.;

    double tot_shape_syst_up[numSyst+1];
    double tot_shape_syst_down[numSyst+1];
    double tot_rate_syst_up[NumRateSys];
    double tot_rate_syst_down[NumRateSys];
    for( int iSyst = 0; iSyst < numSyst+1; ++iSyst ){ 
      tot_shape_syst_up[iSyst] = 0.;
      tot_shape_syst_down[iSyst] = 0.;
    }
    for( int iSyst = 0; iSyst < NumRateSys; ++iSyst ){ 
      tot_rate_syst_up[iSyst] = 0.;
      tot_rate_syst_down[iSyst] = 0.;
    }

    for (int iSample = 0; iSample < numSamples; ++iSample) {

	TString histName = sampleNames[iSample];
	histName += "_CFMlpANN_";
	histName += channelNames[iChan];

	TH1D* hist = (TH1D*)file_input->Get(histName);
	double sumw2 = 0.;
	for( int b=0; b<nbins_ttH; b++ ) sumw2 += hist->GetBinError(b+1)*hist->GetBinError(b+1);
	sumw2_bkg += sumw2;
    
	double sum_hist = hist->Integral();
	sum_bkg += sum_hist;


	double bkg_syst_up_2 = 0;
	double bkg_syst_down_2 = 0;
	for (int iSyst = 0; iSyst < numSyst; ++iSyst) {
	  TString sysName = ( systNames[iSyst].Contains("bShape") ) ? channelNames[iChan] + systNames[iSyst] : systNames[iSyst];

	  if( sysName.Contains("Q2scale_ttH_ttbar0p") ){
	    if( !channelNames[iChan].Contains("e2je2t")
		&& !channelNames[iChan].Contains("ljets_j4_t2")
		&& !channelNames[iChan].Contains("ljets_j4_t3")
		&& !channelNames[iChan].Contains("ljets_j4_t4") ) continue;
	  }
	  if( sysName.Contains("Q2scale_ttH_ttbar1p") ){
	    if( !channelNames[iChan].Contains("ge3t")
		&& !channelNames[iChan].Contains("ljets_j5_t2")
		&& !channelNames[iChan].Contains("ljets_j5_t3")
		&& !channelNames[iChan].Contains("ljets_j5_tge4") ) continue;
	  }
	  if( sysName.Contains("Q2scale_ttH_ttbar2p") ){
	    if( !channelNames[iChan].Contains("ljets_jge6_t2")
		&& !channelNames[iChan].Contains("ljets_jge6_t3")
		&& !channelNames[iChan].Contains("ljets_jge6_tge4") ) continue;
	  }

	  if( sysName.Contains("Q2scale_ttH_ttbar_bb") && !sampleNames[iSample].Contains("ttbarPlusBBbar") ) continue;
	  if( sysName.Contains("Q2scale_ttH_ttbar_cc") && !sampleNames[iSample].Contains("ttbarPlusCCbar") ) continue;

	  if( (sysName.Contains("Q2scale_ttH_ttbar0p") ||
	       sysName.Contains("Q2scale_ttH_ttbar1p") ||
	       sysName.Contains("Q2scale_ttH_ttbar2p")) &&
	      iSample!=0 ) continue;

      
	  TString histName_up =   histName + "_" + sysName + "Up";
	  TString histName_down = histName + "_" + sysName + "Down";

	  TH1D* hist_up   = (TH1D*)file_input->Get(histName_up);
	  TH1D* hist_down = (TH1D*)file_input->Get(histName_down);

	  double syst_up   = hist_up->Integral() - sum_hist;
	  double syst_down = hist_down->Integral() - sum_hist;

	  double syst_up_2 = syst_up*syst_up;
	  double syst_down_2 = syst_down*syst_down;

	  bkg_syst_up_2 += syst_up_2;
	  bkg_syst_down_2 += syst_down_2;

	  tot_shape_syst_up[iSyst] += syst_up;
	  tot_shape_syst_down[iSyst] += syst_down;
	}

	//if (systNames[n] == "Q2scale_ttH_V") {
	TString syst = "-";
	if( sampleNames[iSample] == "wjets" || sampleNames[iSample] == "zjets") {
	  int iJet = 0;
	  if (channelNames[iChan].Contains("e2je2t")) {
	    iJet = 2;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if (channelNames[iChan].Contains("e3je2t")) {
	    iJet = 3;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if (channelNames[iChan].Contains("ge4je2t")) {
	    iJet = 4;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if (channelNames[iChan].Contains("ge3t")) {
	    iJet = 3;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if (channelNames[iChan].Contains("j4_t")) {
	    iJet = 4;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if  (channelNames[iChan].Contains("j5_t")) {
	    iJet = 5;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else if  (channelNames[iChan].Contains("jge6_t")) {
	    iJet = 6;
	    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
	  } else {
	    syst = "-";
	  }
	}

	if( !syst.EqualTo("-") ){
	  double mysyst_val = (syst.Atof() - 1.)*sum_hist;
	  tot_shape_syst_up[numSyst] += mysyst_val;
	  tot_shape_syst_down[numSyst] += mysyst_val;
	  bkg_syst_up_2 += mysyst_val*mysyst_val;
	  bkg_syst_down_2 += mysyst_val*mysyst_val;
	}


	for (int iSyst = 0; iSyst < NumRateSys; ++iSyst) {
	  double mysyst = ( syst_vect[iSyst][iSample]>1. ) ? syst_vect[iSyst][iSample] - 1. : 0.;
	  mysyst *= sum_hist;
	  bkg_syst_up_2 += mysyst*mysyst;
	  bkg_syst_down_2 += mysyst*mysyst;
	  tot_rate_syst_up[iSyst] += mysyst;
	  tot_rate_syst_down[iSyst] += mysyst;
	}


	if( iChan==0 ) table_labels.push_back(sampleNames_latex[iSample]);
 	table_value[iSamp][iChan] = sum_hist;
 	table_value_syst_err[iSamp][iChan] = 0.5* ( sqrt(bkg_syst_up_2) + sqrt(bkg_syst_down_2) );
	table_value_stat_err[iSamp][iChan] = sqrt(sumw2);
	iSamp++;
    }


    double total_up_2 = 0;
    double total_down_2 = 0;

    for (int iSyst = 0; iSyst < numSyst+1; ++iSyst){
      total_up_2 += tot_shape_syst_up[iSyst] * tot_shape_syst_up[iSyst];
      total_down_2 += tot_shape_syst_down[iSyst] * tot_shape_syst_down[iSyst];
    }

    for (int iSyst = 0; iSyst < NumRateSys; ++iSyst){
      total_up_2 += tot_rate_syst_up[iSyst] * tot_rate_syst_up[iSyst];
      total_down_2 += tot_rate_syst_down[iSyst] * tot_rate_syst_down[iSyst];
    }


    if( iChan==0 ) table_labels.push_back("Total bkg");
    table_value[iSamp][iChan] = sum_bkg;
    table_value_syst_err[iSamp][iChan] = 0.5* ( sqrt(total_up_2) + sqrt(total_down_2) );
    table_value_stat_err[iSamp][iChan] = sqrt(sumw2_bkg);
    iSamp++;




    TString data_histName = "data_obs";
    data_histName += "_CFMlpANN_";
    data_histName += channelNames[iChan];

    TH1D* hist_data = (TH1D*)file_input->Get(data_histName);

    if( iChan==0 ) table_labels.push_back("Data");
    table_value[iSamp][iChan] = hist_data->Integral();
    table_value_syst_err[iSamp][iChan] = 0;
    table_value_stat_err[iSamp][iChan] = 0;

  }





  std::cout << " *************** " << std::endl;
  std::cout << "    LJ yields    " << std::endl;
  std::cout << " *************** " << std::endl;



  std::cout << "    \\begin{tabular}{|l|c|c|c|c|c|c|c|} \\hline" << std::endl;
  std::cout << "& $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets \\\\" << std::endl;
  std::cout << "& 2 tags & 3 tags & 3 tags & 3 tags & 4 tags & $\\geq$4 tags & $\\geq$4 tags \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$t\\bar{t}Z$") || table_labels[iSample].EqualTo("$Z+$jets") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$t\\bar{t}W$") )  std::cout << "$t\\bar{t}V$";
    else if( table_labels[iSample].EqualTo("$W+$jets") ) std::cout << "$V+$jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=0; iCat2<NumCats-2; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$t\\bar{t}W$") || table_labels[iSample].EqualTo("$W+$jets") ){
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
    if( table_labels[iSample].Contains("t\\bar{t}H") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;





  std::cout << " *************** " << std::endl;
  std::cout << "    DIL yields   " << std::endl;
  std::cout << " *************** " << std::endl;

  std::cout << "    \\begin{tabular}{|l|c|c|} \\hline" << std::endl;
  std::cout << "& 2 jets + 2 tags & $\\geq$3 tags  \\\\ \\hline \\hline" << std::endl;
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$t\\bar{t}Z$") || table_labels[iSample].EqualTo("$Z+$jets") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$t\\bar{t}W$") )  std::cout << "$t\\bar{t}V$";
    else if( table_labels[iSample].EqualTo("$W+$jets") ) std::cout << "$V+$jets";
    else                                                std::cout << table_labels[iSample];

    for( int iCat2=NumCats-2; iCat2<NumCats; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$t\\bar{t}W$") || table_labels[iSample].EqualTo("$W+$jets") ){
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
    if( table_labels[iSample].Contains("t\\bar{t}H") || table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;
  }
  std::cout << "\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;


}
