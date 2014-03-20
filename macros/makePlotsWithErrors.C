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

root -b -q head.C generateToys.C+'("file.root","preFit",10,1,1)'

Do this each time you want plots

root -b -q head.C makePlotsWithErrors.C+'("file.root")'

*/

class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

void makePlotsWithErrors( TString dataFileName = "", TString dataFileName7TeV = "", bool blind=true, bool useLegend=true, bool debug=false, TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root", bool plotBDTFits_=true, bool plotSoverBCombined_=true, bool plotCategoryYieldCombined_=true ) {


  TString histoFilename = "HistoFiles/generateToys_histo.root";


  bool plotBDTFits = plotBDTFits_;
  bool plotSoverBCombined = plotSoverBCombined_;
  bool plotCategoryYieldCombined = plotCategoryYieldCombined_;

  TFile *file = new TFile(histoFilename);

  TH1D* h_numToys = (TH1D*)file->Get("h_numToys");

  TString imageDir = "Images/Images_2014_02_19_7and8TeV_ttH_allChan_postFitPlots";
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
  /*
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
*/

  color[0] = kBlack;
  //old
  color[bin_diboson] = kCyan;
  color[bin_zjets]   = kGreen+2;
  color[bin_wjets]   = kAzure+1;
  color[bin_ttW]     = kBlue-9;
  color[bin_ttZ]     = color[bin_ttW];

  //new
  color[bin_ttH]     = kBlue+2;

  color[bin_ttjets_bbbar]  = kRed+3;
  color[bin_ttjets_b]      = kRed-2;
  color[bin_ttjets_ccbar]  = kRed+1;
  color[bin_ttjets_other]  = kRed-7;

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

  label[bin_ttH] = "t#bar{t}H(125.6)";

  label[bin_ewk] = "EWK";
  label[bin_ttV] = "t#bar{t} + W,Z";

  label[bin_zjets] = "Z + jets";
  label[bin_wjets] = "W + jets";
  label[bin_diboson] = "WW, WZ, ZZ";



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


  std::vector<TString> axis_labels;
  axis_labels.push_back("#splitline{LJ}{6j2t}");
  axis_labels.push_back("#splitline{LJ}{4j3t}");
  axis_labels.push_back("#splitline{LJ}{5j3t}");
  axis_labels.push_back("#splitline{LJ}{6j3t}");
  axis_labels.push_back("#splitline{LJ}{4j4t}");
  axis_labels.push_back("#splitline{LJ}{5j4t}");
  axis_labels.push_back("#splitline{LJ}{6j4t}");
  axis_labels.push_back("#splitline{DIL}{2j2t}");
  axis_labels.push_back("#splitline{DIL}{3t}");

  axis_labels.push_back("#splitline{LJ}{6j2t}");
  axis_labels.push_back("#splitline{LJ}{4j3t}");
  axis_labels.push_back("#splitline{LJ}{5j3t}");
  axis_labels.push_back("#splitline{LJ}{6j3t}");
  axis_labels.push_back("#splitline{LJ}{4j4t}");
  axis_labels.push_back("#splitline{LJ}{5j4t}");
  axis_labels.push_back("#splitline{LJ}{6j4t}");
  axis_labels.push_back("#splitline{DIL}{3j2t}");
  axis_labels.push_back("#splitline{DIL}{4j2t}");
  axis_labels.push_back("#splitline{DIL}{3t}");
  axis_labels.push_back("#splitline{TAU}{4j1t}");
  axis_labels.push_back("#splitline{TAU}{5j1t}");
  axis_labels.push_back("#splitline{TAU}{6j1t}");
  axis_labels.push_back("#splitline{TAU}{4j2t}");
  axis_labels.push_back("#splitline{TAU}{5j2t}");
  axis_labels.push_back("#splitline{TAU}{6j2t}");



  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open( wsFileName );
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open( fitFileName );
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *dataFile7TeV = TFile::Open(dataFileName7TeV);


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
  fitResults.Add(new LabelInfo("fit_s","postFitS"));


  TString era = "8TeV";

  std::string cmsinfo =   "CMS                             #sqrt{s} = 8 TeV, L = 19.3 fb^{-1}"; //DPUIGH
  //std::string cmsinfo = "CMS Preliminary      #sqrt{s} = 8 TeV, L = 19.3 fb^{-1}"; //SBOUTLE
  //std::string cmsinfo = "CMS Preliminary,  #sqrt{s} = 8 TeV, L = 19.3 fb^{-1}";
  TLatex CMSInfoLatex(0.14, 0.91, cmsinfo.c_str()); //SBOUTLE
  //TLatex CMSInfoLatex(0.47, 0.91, cmsinfo.c_str());
  CMSInfoLatex.SetNDC();
  CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.055); //SBOUTLE


  std::string cmsinfo7TeV =   "CMS                             #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}"; //SBOUTLE
  //std::string cmsinfo7TeV = "CMS Preliminary      #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}"; //SBOUTLE
  //std::string cmsinfo = "CMS Preliminary,  #sqrt{s} = 8 TeV, L = 19.3 fb^{-1}";
  TLatex CMSInfoLatex7TeV(0.14, 0.91, cmsinfo7TeV.c_str()); //SBOUTLE
  //TLatex CMSInfoLatex(0.47, 0.91, cmsinfo.c_str());
  CMSInfoLatex7TeV.SetNDC();
  CMSInfoLatex7TeV.SetTextFont(42);
  CMSInfoLatex7TeV.SetTextSize(0.055); //SBOUTLE


  //CMSInfoLatex.SetTextSize(0.035);

  //Suppress the printout from making plots
  //gErrorIgnoreLevel = 2000;

  //RooFit::RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling) ;
 
  //gSystem->Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so");

  std::vector<std::vector<TString> > BDT_cat_names;
  std::vector<std::vector<double> >  BDT_chi2_value;


  double ratioMax = 2.3;
  double ratioMin = 0.0;

  TIter nextFit(&fitResults);
  LabelInfo *fitRes = 0;

  std::vector<TString> nuisance_names_noMCstat_postFitB;
  std::vector<double>  nuisance_pulls_noMCstat_postFitB;
  std::vector<double>  nuisance_constraints_noMCstat_postFitB;
  std::vector<double>  nuisance_AsymErrorLo_noMCstat_postFitB;
  std::vector<double>  nuisance_AsymErrorHi_noMCstat_postFitB;
  std::vector<double>  nuisance_getMax_noMCstat_postFitB;
  std::vector<double>  nuisance_getMin_noMCstat_postFitB;


  std::vector<TString> nuisance_names_noMCstat_postFitS;
  std::vector<double>  nuisance_pulls_noMCstat_postFitS;
  std::vector<double>  nuisance_constraints_noMCstat_postFitS;
  std::vector<double>  nuisance_AsymErrorLo_noMCstat_postFitS;
  std::vector<double>  nuisance_AsymErrorHi_noMCstat_postFitS;
  std::vector<double>  nuisance_getMax_noMCstat_postFitS;
  std::vector<double>  nuisance_getMin_noMCstat_postFitS;


  while ((fitRes = (LabelInfo *)nextFit())) {

    TString fitName = fitRes->name;
    TString fitLabel = fitRes->label;

    int nToys = 1;
    if( fitLabel.EqualTo("preFit") )        nToys = int(h_numToys->GetBinContent(1) + 0.0000001);
    else if( fitLabel.EqualTo("postFitB") ) nToys = int(h_numToys->GetBinContent(2) + 0.0000001);
    else if( fitLabel.EqualTo("postFitS") ) nToys = int(h_numToys->GetBinContent(3) + 0.0000001);

    if( nToys<1 ){
      std::cout << "\t nToys = " << nToys << ",\t Skipping " << fitLabel << std::endl;
      continue;
    }

    std::cout << " ===> Fit type = " << fitLabel << " (nToys = " << nToys << ")" << std::endl;

    RooFitResult *fitFR = (RooFitResult*)fitFile->Get( fitName );

    std::vector<TString> nuisance_names;
    nuisance_names.clear();

    TH1D* h_pulls = new TH1D("h_pulls","",80,-4,4);
    TH1D* h_pulls_noMCstats = new TH1D("h_pulls_noMCstats","",80,-4,4);
    
    std::vector<TString> nuisance_names_noMCstat;
    std::vector<double>  nuisance_pulls_noMCstat;
    std::vector<double>  nuisance_constraints_noMCstat;
    std::vector<double>  nuisance_AsymErrorLo_noMCstat;
    std::vector<double>  nuisance_AsymErrorHi_noMCstat;
    std::vector<double>  nuisance_getMax_noMCstat;
    std::vector<double>  nuisance_getMin_noMCstat;
    nuisance_names_noMCstat.clear();
    RooArgSet myArgs = fitFR->floatParsFinal();

    TIterator *nextArg(myArgs.createIterator());
    for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
      //RooRealVar *rrv = a;      
      //if( debug ) std::cout << " name = " << rrv->GetName() << ":\t fit value = " << rrv->getVal() << ",\t error = " << rrv->getError() << std::endl;
      TString nuis_name = rrv->GetName();
      nuisance_names.push_back( nuis_name );
      h_pulls->Fill(rrv->getVal());
      if( debug ) std::cout << " name = " << rrv->GetName() << ":\t fit value = " << rrv->getVal() << ",\t error = " << rrv->getError() << ",\t getAsymErrorHi = " << rrv->getAsymErrorHi() << ",\t getAsymErrorLo = " << rrv->getAsymErrorLo() << std::endl;
      if( !nuis_name.Contains("ANNbin") ){
	nuisance_names_noMCstat.push_back( nuis_name );
	nuisance_pulls_noMCstat.push_back( rrv->getVal() );
	nuisance_constraints_noMCstat.push_back( rrv->getError() );
	nuisance_AsymErrorLo_noMCstat.push_back( rrv->getAsymErrorLo() );
	nuisance_AsymErrorHi_noMCstat.push_back( rrv->getAsymErrorHi() );
	nuisance_getMax_noMCstat.push_back( rrv->getMax() );
	nuisance_getMin_noMCstat.push_back( rrv->getMin() );
	h_pulls_noMCstats->Fill(rrv->getVal());
      }
    }

    int NumPulls_noMCstat = int(nuisance_names_noMCstat.size());
    TH1D* h_pulls_all = new TH1D("h_pulls_all",";Index;Pull",NumPulls_noMCstat,0,NumPulls_noMCstat);
    TH2D* h_pulls_vert = new TH2D("h_pulls_vert",";Pull", 100, -2.5, 2.5, NumPulls_noMCstat,0,NumPulls_noMCstat);

    if( fitLabel.EqualTo("postFitB") ){
      for( int iPar=0; iPar<NumPulls_noMCstat; iPar++ ){
	nuisance_names_noMCstat_postFitB.push_back(nuisance_names_noMCstat[iPar]);
	nuisance_pulls_noMCstat_postFitB.push_back(nuisance_pulls_noMCstat[iPar]);
	nuisance_constraints_noMCstat_postFitB.push_back(nuisance_constraints_noMCstat[iPar]);
	nuisance_AsymErrorLo_noMCstat_postFitB.push_back(nuisance_AsymErrorLo_noMCstat[iPar]);
	nuisance_AsymErrorHi_noMCstat_postFitB.push_back(nuisance_AsymErrorHi_noMCstat[iPar]);
	nuisance_getMax_noMCstat_postFitB.push_back(nuisance_getMax_noMCstat[iPar]);
	nuisance_getMin_noMCstat_postFitB.push_back(nuisance_getMin_noMCstat[iPar]);
      }
    }
    else if( fitLabel.EqualTo("postFitS") ){
      for( int iPar=0; iPar<NumPulls_noMCstat; iPar++ ){
	nuisance_names_noMCstat_postFitS.push_back(nuisance_names_noMCstat[iPar]);
	nuisance_pulls_noMCstat_postFitS.push_back(nuisance_pulls_noMCstat[iPar]);
	nuisance_constraints_noMCstat_postFitS.push_back(nuisance_constraints_noMCstat[iPar]);
	nuisance_AsymErrorLo_noMCstat_postFitS.push_back(nuisance_AsymErrorLo_noMCstat[iPar]);
	nuisance_AsymErrorHi_noMCstat_postFitS.push_back(nuisance_AsymErrorHi_noMCstat[iPar]);
	nuisance_getMax_noMCstat_postFitS.push_back(nuisance_getMax_noMCstat[iPar]);
	nuisance_getMin_noMCstat_postFitS.push_back(nuisance_getMin_noMCstat[iPar]);
      }
    }

    double use_y[NumPulls_noMCstat];
    double use_yErr[NumPulls_noMCstat];
    double pulls_noMCstat[NumPulls_noMCstat];
    double constraints_noMCstat[NumPulls_noMCstat];

    for( int iPar=0; iPar<NumPulls_noMCstat; iPar++ ){
      use_y[iPar] = iPar + 0.5;
      use_yErr[iPar] = 0.;

      pulls_noMCstat[iPar] = nuisance_pulls_noMCstat[iPar];
      constraints_noMCstat[iPar] = nuisance_constraints_noMCstat[iPar];

      h_pulls_all->SetBinContent(iPar+1, nuisance_pulls_noMCstat[iPar]);
      h_pulls_all->SetBinError(iPar+1, nuisance_constraints_noMCstat[iPar]);

      h_pulls_vert->GetYaxis()->SetBinLabel(iPar+1,nuisance_names_noMCstat[iPar].Data());
    }

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(NumPulls_noMCstat,pulls_noMCstat,use_y,constraints_noMCstat,constraints_noMCstat,use_yErr,use_yErr); 
    gr->SetMarkerSize(1.2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);

    if( debug ){
      TCanvas *c1 = new TCanvas("c1","");
      h_pulls->Draw("hist");
      c1->Print(debugDir+"/pulls_"+fitLabel+".png");

      h_pulls_noMCstats->Draw("hist");
      c1->Print(debugDir+"/pulls_noMCstats_"+fitLabel+".png");

      h_pulls_all->SetStats(0);
      h_pulls_all->GetYaxis()->SetRangeUser(-2.5,2.5);
      h_pulls_all->SetMarkerStyle(20);
      h_pulls_all->SetLineWidth(2);
      h_pulls_all->Draw("pe1");

      TLine* myLine;
      myLine = new TLine(h_pulls_all->GetXaxis()->GetXmin(), 0., h_pulls_all->GetXaxis()->GetXmax(), 0.);
      myLine->SetLineStyle(7);
      myLine->Draw();

      c1->Print(debugDir+"/pulls_all_"+fitLabel+".png");


      TCanvas *c2 = new TCanvas("c2","",900,800);
      h_pulls_vert->SetStats(0);
      h_pulls_vert->Draw();
      gr->Draw("pe1same");

      TLine* myLine1;
      myLine1 = new TLine( 0., h_pulls_vert->GetYaxis()->GetXmin(), 0., h_pulls_vert->GetYaxis()->GetXmax());
      myLine1->SetLineStyle(7);
      myLine1->Draw();

      c2->GetPad(0)->SetLeftMargin(0.23);
      c2->GetPad(0)->SetTopMargin(0.05);
      c2->GetPad(0)->SetRightMargin(0.05);
      c2->Print(debugDir+"/pulls_vert_"+fitLabel+".png");

      delete c1;
      delete c2;
      delete myLine;
      delete myLine1;
    }
    delete gr;
    delete h_pulls_vert;
    delete h_pulls_all;
    delete h_pulls;
    delete h_pulls_noMCstats;


    std::cout << " Number of nuisances = " << int(nuisance_names.size()) << std::endl;
    std::cout << " Number of nuisances = " << int(nuisance_names_noMCstat.size()) << " (not including MCstats)" << std::endl;

    int NumPars = int(nuisance_names.size());

    std::string fitinfo = "";
    if( fitLabel.EqualTo("preFit") )        fitinfo = "Pre-Fit";
    else if( fitLabel.EqualTo("postFitB") ) fitinfo = "Post-Fit (B)";
    else if( fitLabel.EqualTo("postFitS") ) fitinfo = "Post-Fit (S+B)";


    TH1D** h_nuis = new TH1D*[NumPars];
    //TH1D* h_nuis[NumPars];
    for( int iPar=0; iPar<NumPars; iPar++ ){
      TString nuis_name = nuisance_names[iPar];
      h_nuis[iPar] = (TH1D*)file->Get("h_nuis_"+nuisance_names[iPar]+"_"+fitLabel);
      h_nuis[iPar]->GetXaxis()->SetTitle(nuis_name);
    }


    //TLatex FitinfoLatex(0.45, 0.97, fitinfo.c_str());
    TLatex FitinfoLatex(0.44, 0.97, fitinfo.c_str());//ccn
    FitinfoLatex.SetNDC();
    FitinfoLatex.SetTextFont(42);
    //FitinfoLatex.SetTextSize(0.04); //ccn: use for S/B plot all others use 0.045
    FitinfoLatex.SetTextSize(0.045);//ccn





    //Now get the fit central value
    w->saveSnapshot(fitLabel,RooArgSet(fitFR->floatParsFinal()),kTRUE);
    w->loadSnapshot(fitLabel);


    int numBins_most = 20;
    TH1D* fluctHist[numCats][numBins_most];
    RooAddition** fitTotal_cat = new RooAddition*[numCats];
    int numBins[numCats];

    TIter nextCat1(&categories);
    LabelInfo *category1 = 0;
    int iCat1=0;
    while ((category1 = (LabelInfo *)nextCat1())) {
      TString catLabel = category1->label;
      TString catName  = category1->name;

      TString dataCatName = catName;
      if( dataCatName.Contains("ch1_") ) dataCatName.ReplaceAll("ch1_","");
      if( dataCatName.Contains("hbb_7TeV_") ) dataCatName.ReplaceAll("hbb_7TeV_","");
      if( dataCatName.Contains("hbb_8TeV_") ) dataCatName.ReplaceAll("hbb_8TeV_","");
      if( dataCatName.Contains("htt_8TeV_") ) dataCatName.ReplaceAll("htt_8TeV_","");

      //Roofit plots histogram vs "bin number;" turn this into a histogram over ANN output
      //TH1 *dataHist = (TH1 *)dataFile->Get("data_obs_MVA_"+catName);
      TH1 *dataHist = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get("data_obs_CFMlpANN_"+dataCatName) : (TH1 *)dataFile->Get("data_obs_MVA_"+dataCatName);
      TH1 *fitHist = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_"+fitLabel+"_Clone_temp");
      numBins[iCat1] = fitHist->GetNbinsX();

      //Calculate the total in the fit...
      RooArgSet fitArgs;
      TIter nextProcess(&processes);
      LabelInfo *process = 0;
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
      for( int iBin=0; iBin<numBins[iCat1]; iBin++ ) fluctHist[iCat1][iBin] = (TH1D*)file->Get(Form("fluctHist_%s_%s_%d",fitLabel.Data(),catName.Data(),iBin));
      iCat1++;
    }


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
      c1->Print(debugDir+"/correlations_nuisances_2D_"+fitLabel+".pdf");


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


    if( plotSoverBCombined ){
      TH1D* h_data_SoverB = (TH1D*)file->Get(Form("data_SoverB_%s",fitLabel.Data()));
      TH1D* h_bkg_SoverB = (TH1D*)file->Get(Form("bkg_SoverB_%s",fitLabel.Data()));
      TH1D* h_sig_SoverB = (TH1D*)file->Get(Form("sig_SoverB_%s",fitLabel.Data()));

      h_data_SoverB->SetLineColor(1);
      h_data_SoverB->SetLineWidth(2);
      h_data_SoverB->SetMarkerStyle(20);
      h_data_SoverB->SetMarkerSize(0.75);


      h_sig_SoverB->SetFillColor(kRed);

      THStack *hs = new THStack("hs_SoverB_"+fitLabel,"");

      hs->Add(h_bkg_SoverB);
      hs->Add(h_sig_SoverB);


      TH1D* h_bkgToFit_SoverB = (TH1D*)h_bkg_SoverB->Clone(Form("bkgToFit_SoverB_%s",fitLabel.Data()));
      if( fitLabel.EqualTo("postFitS") ) h_bkgToFit_SoverB->Add(h_sig_SoverB);

      int numBinsSoverB = 20;
      TH1D* fluctHist_SoverB[numBinsSoverB];
      for( int iBin=0; iBin<numBinsSoverB; iBin++ ) fluctHist_SoverB[iBin] = (TH1D*)file->Get(Form("fluctHist_SoverB_%s_%d",fitLabel.Data(),iBin));



      TH1 *errHist_SoverB_1sig = (TH1 *)h_data_SoverB->Clone("errHist_SoverB_1sig_"+fitLabel);
      errHist_SoverB_1sig->Reset();
      errHist_SoverB_1sig->SetFillColor(1);
      errHist_SoverB_1sig->SetFillStyle(3654);
      errHist_SoverB_1sig->SetLineColor(1);
      errHist_SoverB_1sig->SetMarkerColor(0);
      errHist_SoverB_1sig->SetMarkerStyle(1);
      errHist_SoverB_1sig->SetMarkerSize(0.);

      for (int xBin=0; xBin<numBinsSoverB; xBin++) {

	//int yBinMin = fluctHist_SoverB[xBin]->FindBin(h_data_SoverB->GetBinContent(xBin+1));
	int yBinMin = fluctHist_SoverB[xBin]->FindBin(h_bkg_SoverB->GetBinContent(xBin+1));
	int yBinMax = yBinMin; //To start...

	double sum = fluctHist_SoverB[xBin]->GetBinContent(yBinMin);

	//printf("\t xBin=%d, data=%.0f, findBin=%d, sum=%.1f\t \n",xBin,h_data_SoverB->GetBinContent(xBin+1),yBinMin,sum);

	while (sum < 0.68269*nToys && yBinMin > 0 && yBinMax <= fluctHist_SoverB[xBin]->GetNbinsX()) {

	  double contMin = fluctHist_SoverB[xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist_SoverB[xBin]->GetBinContent(yBinMax+1);

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
	if (yBinMax == fluctHist_SoverB[xBin]->GetNbinsX()+1) std::cerr << "WARNING (1sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist_SoverB[xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist_SoverB[xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_SoverB_1sig->SetBinContent(xBin+1,yAve);
	errHist_SoverB_1sig->SetBinError(xBin+1,yMax-yAve);
      }



      TH1D* h_bkg_SoverB_blind = (TH1D*)h_bkgToFit_SoverB->Clone("h_bkgToFit_SoverB_blind");
      if( blind ){
	for( int iBin=0; iBin<numBinsSoverB; iBin++ ){
	  double lowEdge = h_data_SoverB->GetBinLowEdge(iBin+1);
	  if( lowEdge>-1.5 ){
	    h_data_SoverB->SetBinContent(iBin+1,0.);
	    h_bkg_SoverB_blind->SetBinContent(iBin+1,0.);
	  }
	}
      }

      //Now calculate the ratio
      TH1 *ratio_SoverB = (TH1 *)h_data_SoverB->Clone("ratio_SoverB_"+fitLabel);
      ratio_SoverB->SetMarkerStyle(1);

      TH1 *ratioErr_SoverB_1sig = (TH1 *)ratio_SoverB->Clone("ratioErr_SoverB_1sig_"+fitLabel);
      ratioErr_SoverB_1sig->Reset();
      ratioErr_SoverB_1sig->SetFillColor(kGreen);
      ratioErr_SoverB_1sig->SetLineColor(0);
      ratioErr_SoverB_1sig->SetLineWidth(0);
      ratioErr_SoverB_1sig->SetMarkerStyle(0);
      ratioErr_SoverB_1sig->SetMarkerColor(0);
      ratioErr_SoverB_1sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numBinsSoverB; iBin++) {

	double dataCentral = h_data_SoverB->GetBinContent(iBin);
	double mcCentral = h_bkg_SoverB->GetBinContent(iBin);
	double mcUp = errHist_SoverB_1sig->GetBinContent(iBin)+errHist_SoverB_1sig->GetBinError(iBin);
	double mcDown = errHist_SoverB_1sig->GetBinContent(iBin)-errHist_SoverB_1sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_SoverB_1sig->SetBinContent(iBin,ratioAve);
	ratioErr_SoverB_1sig->SetBinError(iBin, ratioUp-ratioAve);

	double mcCentral_maySig = h_bkgToFit_SoverB->GetBinContent(iBin);
	double myratio = dataCentral/mcCentral_maySig;
	double myratio_err = sqrt( dataCentral ) / mcCentral_maySig;

	ratio_SoverB->SetBinContent(iBin,myratio);
	ratio_SoverB->SetBinError(iBin,myratio_err);

	if( (myratio>ratioMax) && ((myratio - myratio_err)<ratioMax) ){
	  double minner = myratio - myratio_err;
	  ratio_SoverB->SetBinContent(iBin,ratioMax-0.0001);
	  ratio_SoverB->SetBinError(iBin,ratioMax-0.0001-minner);
	}
      }


      TH1 *errHist_SoverB_2sig = (TH1 *)h_data_SoverB->Clone("errHist_SoverB_2sig_"+fitLabel);
      errHist_SoverB_2sig->Reset();
      errHist_SoverB_2sig->SetFillColor(1);
      errHist_SoverB_2sig->SetFillStyle(3654);
      errHist_SoverB_2sig->SetLineColor(1);
      errHist_SoverB_2sig->SetMarkerColor(0);
      errHist_SoverB_2sig->SetMarkerStyle(1);
      errHist_SoverB_2sig->SetMarkerSize(0.);

      for (int xBin=0; xBin<numBinsSoverB; xBin++) {

	int yBinMin = fluctHist_SoverB[xBin]->FindBin(h_bkg_SoverB->GetBinContent(xBin+1));
	int yBinMax = yBinMin; //To start...

	double sum = fluctHist_SoverB[xBin]->GetBinContent(yBinMin);

	while (sum < 0.954499736*nToys && yBinMin > 0 && yBinMax <= fluctHist_SoverB[xBin]->GetNbinsX()) {

	  double contMin = fluctHist_SoverB[xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist_SoverB[xBin]->GetBinContent(yBinMax+1);

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

	if (yBinMin == 0  && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to underflow for bin " << xBin+1 << "!" << std::endl;
	if (yBinMax == (fluctHist_SoverB[xBin]->GetNbinsX()+1) && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist_SoverB[xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist_SoverB[xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_SoverB_2sig->SetBinContent(xBin+1,yAve);
	errHist_SoverB_2sig->SetBinError(xBin+1,yMax-yAve);
      }

      TH1 *ratioErr_SoverB_2sig = (TH1 *)ratio_SoverB->Clone("ratioErr_SoverB_2sig_"+fitLabel);
      ratioErr_SoverB_2sig->Reset();
      ratioErr_SoverB_2sig->SetFillColor(kYellow);
      ratioErr_SoverB_2sig->SetLineColor(0);
      ratioErr_SoverB_2sig->SetLineWidth(0);
      ratioErr_SoverB_2sig->SetMarkerStyle(0);
      ratioErr_SoverB_2sig->SetMarkerColor(0);
      ratioErr_SoverB_2sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numBinsSoverB; iBin++) {
	double mcCentral = h_bkg_SoverB->GetBinContent(iBin);
	double mcUp = errHist_SoverB_2sig->GetBinContent(iBin)+errHist_SoverB_2sig->GetBinError(iBin);
	double mcDown = errHist_SoverB_2sig->GetBinContent(iBin)-errHist_SoverB_2sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_SoverB_2sig->SetBinContent(iBin,ratioAve);
	ratioErr_SoverB_2sig->SetBinError(iBin, ratioUp-ratioAve);
      }


      if( debug ) {
	for( int xBin = 1; xBin <= numBinsSoverB; xBin++ ) {
	  std::cout << " xBin = " << xBin << ",\t dataHist = " << h_data_SoverB->GetBinContent(xBin) << ",\t fitHist = " << h_bkg_SoverB->GetBinContent(xBin) << ",\t errHist_1sig = " << errHist_SoverB_1sig->GetBinContent(xBin) << " +/- " << errHist_SoverB_1sig->GetBinError(xBin) << ",\t ratio = " << ratio_SoverB->GetBinContent(xBin) << ",\t ratioErr_1sig = " << ratioErr_SoverB_1sig->GetBinContent(xBin) << " +/- " << ratioErr_SoverB_1sig->GetBinError(xBin) << std::endl;
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

      ratioErr_SoverB_1sig->SetMinimum(ratioMin);
      ratioErr_SoverB_1sig->SetMaximum(ratioMax);
      ratioErr_SoverB_1sig->GetYaxis()->SetNdivisions(50000+404);
      ratioErr_SoverB_1sig->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratioErr_SoverB_1sig->GetXaxis()->SetLabelSize(0.1); //make y label bigger
      ratioErr_SoverB_1sig->GetXaxis()->SetTitleOffset(0.9); //SBOUTLE
      //ratioErr_SoverB_1sig->GetXaxis()->SetTitleOffset(1.1);
      ratioErr_SoverB_1sig->GetXaxis()->SetTitle(h_data_SoverB->GetXaxis()->GetTitle()); //make y label bigger
      ratioErr_SoverB_1sig->GetXaxis()->SetLabelSize(0.12); //SBOUTLE
      //ratioErr_SoverB_1sig->GetXaxis()->SetLabelSize(0.1);
      ratioErr_SoverB_1sig->GetXaxis()->SetLabelOffset(0.02);//SBOUTLE
      //ratioErr_SoverB_1sig->GetXaxis()->SetLabelOffset(0.04);
      ratioErr_SoverB_1sig->GetXaxis()->SetTitleSize(0.14); //SBOUTLE
      //ratioErr_SoverB_1sig->GetXaxis()->SetTitleSize(0.12);
      ratioErr_SoverB_1sig->GetYaxis()->SetTitle("Data/MC");
      ratioErr_SoverB_1sig->GetYaxis()->SetTitleSize(0.1);
      ratioErr_SoverB_1sig->GetYaxis()->CenterTitle();
      ratioErr_SoverB_1sig->GetYaxis()->SetTitleOffset(0.45);
      myC->cd(2);
      gPad->SetTopMargin(small);
      gPad->SetTickx();
      gPad->Modified();

      
      TLegend *legend = new TLegend(0.2,0.8,0.94,0.89);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);

      legend->SetNColumns(3);

      legend->AddEntry(h_data_SoverB,label[bin_data].Data(),"lpe");
      legend->AddEntry(errHist_SoverB_1sig,"Sum MC","f");
      legend->AddEntry(h_sig_SoverB,"t#bar{t}H125.6","f");



      TLegend *legend_ratio = new TLegend(0.14,0.87,0.48,0.97);

      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetShadowColor(kWhite);
      legend_ratio->SetTextFont(42);
      //legend_ratio->SetTextSize(0.055);
      legend_ratio->SetTextSize(0.10);//ccn

      if( !fitLabel.EqualTo("preFit") ){
	legend_ratio->SetNColumns(2);

	legend_ratio->AddEntry(ratioErr_SoverB_1sig,"Fit #pm1#sigma","f");
	legend_ratio->AddEntry(ratioErr_SoverB_2sig,"Fit #pm2#sigma","f");
      }
      else {
	legend_ratio->AddEntry(ratioErr_SoverB_1sig,"Stat+Syst Uncertainty","f");
      }

      ratioErr_SoverB_1sig->SetTitle(";Log_{10}(S/B);Data/MC");


      h_data_SoverB->SetTitle(";;Events");
      h_data_SoverB->GetYaxis()->SetTitleSize(0.05);
      h_data_SoverB->GetYaxis()->SetTitleOffset(1.1); //SBOUTLE
      //h_data_SoverB->GetYaxis()->SetTitleOffset(1.);

      double xmin = h_data_SoverB->GetBinLowEdge(1);
      double xmax = h_data_SoverB->GetBinLowEdge(numBinsSoverB) + h_data_SoverB->GetBinWidth(numBinsSoverB);



      myC->cd(1);
      myC->GetPad(1)->SetLogy(1);

      int maxBin_data = h_data_SoverB->GetMaximumBin();
      int maxBin_bkg  = errHist_SoverB_1sig->GetMaximumBin();

      double max_data = h_data_SoverB->GetBinContent(maxBin_data) + h_data_SoverB->GetBinError(maxBin_data);
      double max_bkg  = errHist_SoverB_1sig->GetBinContent(maxBin_bkg) + errHist_SoverB_1sig->GetBinError(maxBin_bkg);

      double histMax = TMath::Max( max_data, max_bkg );

      if( useLegend ) histMax = 10. * histMax;
      else            histMax = 10. * histMax;

      h_data_SoverB->SetMaximum(histMax);
      errHist_SoverB_1sig->SetMaximum(histMax);

      double chi2 = h_data_SoverB->Chi2Test(h_bkg_SoverB_blind,"UWP");

      TString str_chi2 = ( !fitLabel.EqualTo("postFitS") ) ? Form("B-only p-value (#chi^{2}) = %.3f",chi2) : Form("S+B p-value (#chi^{2}) = %.3f",chi2);
      TLatex Chi2Latex(0.59, 0.91, str_chi2);
      Chi2Latex.SetNDC();
      Chi2Latex.SetTextFont(42);
      //Chi2Latex.SetTextSize(0.055);
      Chi2Latex.SetTextSize(0.09);//ccn


      TString selectioninfo = "ttH(bb,#tau#tau): LJ + DIL + TAU";//"Lepton + #geq6 jets + 2 tags";
      //TLatex SELECTIONInfoLatex(0.12, 0.91, selectioninfo);
      TLatex SELECTIONInfoLatex(0.60, 0.84, selectioninfo);//ccn
      SELECTIONInfoLatex.SetNDC();
      SELECTIONInfoLatex.SetTextFont(42);
      //SELECTIONInfoLatex.SetTextSize(0.035);
      SELECTIONInfoLatex.SetTextSize(0.05);//ccn



      h_data_SoverB->SetMinimum(5);
      h_data_SoverB->Draw("e1");
      hs->Draw("histsame");
      errHist_SoverB_1sig->Draw("e2same");
      h_data_SoverB->Draw("e1SAME");

      if( useLegend ) legend->Draw();

      CMSInfoLatex.Draw();
      //SELECTIONInfoLatex.Draw();
      FitinfoLatex.Draw();


      myC->cd(2);
      ratioErr_SoverB_1sig->Draw("e2");
      if( !fitLabel.EqualTo("preFit") ){
	ratioErr_SoverB_2sig->Draw("e2same");
	ratioErr_SoverB_1sig->Draw("e2same");
      }
      ratio_SoverB->Draw("e1SAME");

      myC->GetPad(1)->RedrawAxis();
      myC->GetPad(2)->RedrawAxis();


      TLine* myLine;
      myLine = new TLine(xmin, 1.000, xmax, 1.000);
      myLine->Draw();
      legend_ratio->Draw();

      Chi2Latex.Draw();

      myC->SaveAs(imageDir+"/dataToMC_"+era+"_SoverB_"+fitLabel+".png");
      myC->SaveAs(imageDir+"/dataToMC_"+era+"_SoverB_"+fitLabel+".pdf");

      delete hs;
      delete myLine;
      delete myC;


      if( debug ){

	TLine* l_nom = new TLine( 0,0,1,1 );
	TLine* l_ave = new TLine( 0,0,1,1 );
	TLine* l_ave_p1 = new TLine( 0,0,1,1 );
	TLine* l_ave_m1 = new TLine( 0,0,1,1 );
	TLine* l_ave_p2 = new TLine( 0,0,1,1 );
	TLine* l_ave_m2 = new TLine( 0,0,1,1 );

	l_nom->SetLineColor(kBlack);
	l_ave->SetLineColor(kBlue);
	l_ave_p1->SetLineColor(kGreen+1);
	l_ave_m1->SetLineColor(kGreen+1);
	l_ave_p2->SetLineColor(kRed+1);
	l_ave_m2->SetLineColor(kRed+1);

	l_nom->SetLineWidth(2);
	l_ave->SetLineWidth(2);
	l_ave_p1->SetLineWidth(2);
	l_ave_m1->SetLineWidth(2);
	l_ave_p2->SetLineWidth(2);
	l_ave_m2->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1","",900,800);
	for( int iBin=0; iBin<numBinsSoverB; iBin++ ){

	  TString s_bin = Form("%d",iBin+1);
	  TString proj_name = "projection_bin"+s_bin+"_"+fitLabel;

	  TH1D* h_proj = (TH1D*)fluctHist_SoverB[iBin]->Clone(proj_name);
	  int rebin = 1000;
	  if( !fitLabel.EqualTo("preFit")) rebin = 100;
	  h_proj->Rebin(rebin);

	  int nProjBins = h_proj->GetNbinsX();
	  std::cout << "  h_proj nbins = " << nProjBins << std::endl;

	  double ave = errHist_SoverB_1sig->GetBinContent(iBin+1);
	  double aveErr = errHist_SoverB_1sig->GetBinError(iBin+1);
	  double nom = h_bkg_SoverB->GetBinContent(iBin+1);

	  double ave_2s = errHist_SoverB_2sig->GetBinContent(iBin+1);
	  double aveErr_2s = errHist_SoverB_2sig->GetBinError(iBin+1);

	  double scaleBand = 4;
	  double xminProj = std::max( h_proj->GetXaxis()->GetXmin(), ave-scaleBand*aveErr );
	  double xmaxProj = std::min( h_proj->GetXaxis()->GetXmax(), ave+scaleBand*aveErr );

	  double ymaxProj = 1.05 * h_proj->GetMaximum();

	  double integral_total = h_proj->Integral();
	  double integral_68 = h_proj->Integral( h_proj->FindBin(ave-aveErr), h_proj->FindBin(ave+aveErr) );
	  double integral_95 = h_proj->Integral( h_proj->FindBin(ave_2s-aveErr_2s), h_proj->FindBin(ave_2s+aveErr_2s) );

	  printf("\t bin = %d,\t integral_total = %.1f,\t integral_68 = %.1f,\t integral_95 = %.1f,\t 68/total = %.3f,\t 95/total = %.3f \n",iBin+1,integral_total, integral_68, integral_95,integral_68/integral_total,integral_95/integral_total);

	  h_proj->GetXaxis()->SetRangeUser(xminProj,xmaxProj);
	  h_proj->GetYaxis()->SetRangeUser(0,ymaxProj);

	  h_proj->SetLineColor(kBlack);
	  h_proj->Draw("hist");

	  l_nom->DrawLine(nom,0,nom,ymaxProj);
	  l_ave->DrawLine(ave,0,ave,ymaxProj);
	  l_ave_p1->DrawLine(ave+aveErr,0,ave+aveErr,ymaxProj);
	  l_ave_m1->DrawLine(ave-aveErr,0,ave-aveErr,ymaxProj);
	  if( !fitLabel.EqualTo("preFit") ){
	    l_ave_p2->DrawLine(ave_2s+aveErr_2s,0,ave_2s+aveErr_2s,ymaxProj);
	    l_ave_m2->DrawLine(ave_2s-aveErr_2s,0,ave_2s-aveErr_2s,ymaxProj);
	  }

	  printf("  bin = %d,\t nom = %.2f,\t ave = %.2f,\t ave - err = %.2f,\t ave+err = %.2f \n",iBin+1,nom,ave,ave-aveErr,ave+aveErr );

	  c1->Print(debugDir+"/projection_SoverB_"+fitLabel+"_bin"+s_bin+".png");
	}

	delete c1;
	delete l_nom;
	delete l_ave;
	delete l_ave_p1;
	delete l_ave_m1;
	delete l_ave_p2;
	delete l_ave_m2;
      }

    }




    if( plotCategoryYieldCombined ){
      TH1D* h_data_category_yield = (TH1D*)file->Get(Form("data_category_yield_%s",fitLabel.Data()));
      TH1D* h_bkg_category_yield = (TH1D*)file->Get(Form("bkg_category_yield_%s",fitLabel.Data()));
      TH1D* h_sig_category_yield = (TH1D*)file->Get(Form("sig_category_yield_%s",fitLabel.Data()));

      h_data_category_yield->SetLineColor(1);
      h_data_category_yield->SetLineWidth(4);
      //h_data_category_yield->SetLineWidth(2);
      h_data_category_yield->SetMarkerStyle(20);
      //h_data_category_yield->SetMarkerSize(0.75);
      h_data_category_yield->SetMarkerSize(1.4);//ccn


      TH1* h_bkg[NumBkgs];

      TIter nextProcess(&processes);
      LabelInfo *process = 0;
      int nProc = 0;
      //std::cout << "\t Filling individual sample components" << std::endl;
      while ((process = (LabelInfo *)nextProcess())) {
	TString procName = process->name;
	h_bkg[nProc] = (TH1 *)h_data_category_yield->Clone("pdf_bin_"+procName+"_"+fitLabel+"Clone_cat_yield");
	h_bkg[nProc]->Reset();

	TIter nextCat3(&categories);
	LabelInfo *category3 = 0;
	int iCat3=0;
	while ((category3 = (LabelInfo *)nextCat3())) {
	  TString catLabel = category3->label;
	  TString catName  = category3->name;

	  //Construct the name of the normalization variable.
	  TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;

	  RooAbsReal *fitN = w->function(normName);

	  if( fitN ){
	    h_bkg[nProc]->Fill(iCat3,fitN->getVal());
	  }
	  iCat3++;
	}
	nProc++;
      }


      TH1D* h_ttH = (TH1D*)h_sig_category_yield->Clone("h_ttH_category_yield"+fitLabel);
      h_ttH->SetLineColor(color[bin_ttH]);
      h_ttH->SetLineWidth(4);


      TIter nextCat3(&categories);
      LabelInfo *category3 = 0;
      double integral_ttH_original = 0;
      while ((category3 = (LabelInfo *)nextCat3())) {
	TString catLabel = category3->label;
	TString catName  = category3->name;

	TString dataCatName = catName;
	if( dataCatName.Contains("ch1_") ) dataCatName.ReplaceAll("ch1_","");
	if( dataCatName.Contains("hbb_7TeV_") ) dataCatName.ReplaceAll("hbb_7TeV_","");
	if( dataCatName.Contains("hbb_8TeV_") ) dataCatName.ReplaceAll("hbb_8TeV_","");
	if( dataCatName.Contains("htt_8TeV_") ) dataCatName.ReplaceAll("htt_8TeV_","");

	TH1 *sigHist_original = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get("ttH125.6_CFMlpANN_"+dataCatName)->Clone("ttH125.6_original_"+dataCatName+"_"+fitLabel) : (TH1 *)dataFile->Get("ttH125.6_MVA_"+dataCatName)->Clone("ttH125.6_original_"+dataCatName+"_"+fitLabel);
	integral_ttH_original += sigHist_original->Integral();
      }

      if( blind && fitLabel.EqualTo("postFitS") ){
	double integral_ttH_post = h_ttH->Integral();
	h_ttH->Scale( integral_ttH_original / integral_ttH_post );
      }

      TH1D* h_ewk_bkg = (TH1D*)h_bkg[bin_wjets]->Clone("h_ewk_category_yield_"+fitLabel);
      h_ewk_bkg->Add(h_bkg[bin_zjets]);
      h_ewk_bkg->Add(h_bkg[bin_diboson]);

      TH1D* h_ttbarV_bkg = (TH1D*)h_bkg[bin_ttW]->Clone("h_ttbarV_category_yield_"+fitLabel);
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


      THStack *hs = new THStack("hs_category_yield_"+fitLabel,"");

      hs->Add(h_ewk_bkg);
      hs->Add(h_ttbarV_bkg);
      hs->Add(h_bkg[bin_singlet]);
      hs->Add(h_bkg[bin_ttjets_bbbar]);
      hs->Add(h_bkg[bin_ttjets_b]);
      hs->Add(h_bkg[bin_ttjets_ccbar]);
      hs->Add(h_bkg[bin_ttjets_other]);






      TH1D* h_bkgToFit_category_yield = (TH1D*)h_bkg_category_yield->Clone(Form("bkgToFit_category_yield_%s",fitLabel.Data()));
      if( fitLabel.EqualTo("postFitS") ) h_bkgToFit_category_yield->Add(h_sig_category_yield);

      TH1D* fluctHist_category_yield[numCats];
      for( int iBin=0; iBin<numCats; iBin++ ) fluctHist_category_yield[iBin] = (TH1D*)file->Get(Form("fluctHist_category_yield_%s_%d",fitLabel.Data(),iBin));



      TH1 *errHist_category_yield_1sig = (TH1 *)h_data_category_yield->Clone("errHist_category_yield_1sig_"+fitLabel);
      errHist_category_yield_1sig->Reset();
      errHist_category_yield_1sig->SetFillColor(1);
      errHist_category_yield_1sig->SetFillStyle(3654);
      errHist_category_yield_1sig->SetLineColor(1);
      errHist_category_yield_1sig->SetMarkerColor(0);
      errHist_category_yield_1sig->SetMarkerStyle(1);
      errHist_category_yield_1sig->SetMarkerSize(0.);

      for (int xBin=0; xBin<numCats; xBin++) {

	//int yBinMin = fluctHist_category_yield[xBin]->FindBin(h_data_category_yield->GetBinContent(xBin+1));
	int yBinMin = fluctHist_category_yield[xBin]->FindBin(h_bkg_category_yield->GetBinContent(xBin+1));
	int yBinMax = yBinMin; //To start...

	double sum = fluctHist_category_yield[xBin]->GetBinContent(yBinMin);

	//printf("\t xBin=%d, data=%.0f, findBin=%d, sum=%.1f\t \n",xBin,h_data_category_yield->GetBinContent(xBin+1),yBinMin,sum);

	while (sum < 0.68269*nToys && yBinMin > 0 && yBinMax <= fluctHist_category_yield[xBin]->GetNbinsX()) {

	  double contMin = fluctHist_category_yield[xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist_category_yield[xBin]->GetBinContent(yBinMax+1);

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
	if (yBinMax == fluctHist_category_yield[xBin]->GetNbinsX()+1) std::cerr << "WARNING (1sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist_category_yield[xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist_category_yield[xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_category_yield_1sig->SetBinContent(xBin+1,yAve);
	errHist_category_yield_1sig->SetBinError(xBin+1,yMax-yAve);
      }



      TH1D* h_bkg_category_yield_blind = (TH1D*)h_bkgToFit_category_yield->Clone("h_bkgToFit_category_yield_blind");
//       if( blind ){
// 	for( int iBin=0; iBin<numCats; iBin++ ){
// 	  double lowEdge = h_data_category_yield->GetBinLowEdge(iBin+1);
// 	  if( lowEdge>-1.5 ){
// 	    h_data_category_yield->SetBinContent(iBin+1,0.);
// 	    h_bkg_category_yield_blind->SetBinContent(iBin+1,0.);
// 	  }
// 	}
//       }

      //Now calculate the ratio
      TH1 *ratio_category_yield = (TH1 *)h_data_category_yield->Clone("ratio_category_yield_"+fitLabel);
      ratio_category_yield->SetMarkerStyle(1);

      TH1 *ratioErr_category_yield_1sig = (TH1 *)ratio_category_yield->Clone("ratioErr_category_yield_1sig_"+fitLabel);
      ratioErr_category_yield_1sig->Reset();
      ratioErr_category_yield_1sig->SetFillColor(kGreen);
      ratioErr_category_yield_1sig->SetLineColor(0);
      ratioErr_category_yield_1sig->SetLineWidth(0);
      ratioErr_category_yield_1sig->SetMarkerStyle(0);
      ratioErr_category_yield_1sig->SetMarkerColor(0);
      ratioErr_category_yield_1sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numCats; iBin++) {

	double dataCentral = h_data_category_yield->GetBinContent(iBin);
	double mcCentral = h_bkg_category_yield->GetBinContent(iBin);
	double mcUp = errHist_category_yield_1sig->GetBinContent(iBin)+errHist_category_yield_1sig->GetBinError(iBin);
	double mcDown = errHist_category_yield_1sig->GetBinContent(iBin)-errHist_category_yield_1sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_category_yield_1sig->SetBinContent(iBin,ratioAve);
	ratioErr_category_yield_1sig->SetBinError(iBin, ratioUp-ratioAve);

	double mcCentral_maySig = h_bkgToFit_category_yield->GetBinContent(iBin);
	double myratio = dataCentral/mcCentral_maySig;
	double myratio_err = sqrt( dataCentral ) / mcCentral_maySig;

	ratio_category_yield->SetBinContent(iBin,myratio);
	ratio_category_yield->SetBinError(iBin,myratio_err);

	if( (myratio>ratioMax) && ((myratio - myratio_err)<ratioMax) ){
	  double minner = myratio - myratio_err;
	  ratio_category_yield->SetBinContent(iBin,ratioMax-0.0001);
	  ratio_category_yield->SetBinError(iBin,ratioMax-0.0001-minner);
	}
      }


      TH1 *errHist_category_yield_2sig = (TH1 *)h_data_category_yield->Clone("errHist_category_yield_2sig_"+fitLabel);
      errHist_category_yield_2sig->Reset();
      errHist_category_yield_2sig->SetFillColor(1);
      errHist_category_yield_2sig->SetFillStyle(3654);
      errHist_category_yield_2sig->SetLineColor(1);
      errHist_category_yield_2sig->SetMarkerColor(0);
      errHist_category_yield_2sig->SetMarkerStyle(1);
      errHist_category_yield_2sig->SetMarkerSize(0.);

      for (int xBin=0; xBin<numCats; xBin++) {

	int yBinMin = fluctHist_category_yield[xBin]->FindBin(h_bkg_category_yield->GetBinContent(xBin+1));
	int yBinMax = yBinMin; //To start...

	double sum = fluctHist_category_yield[xBin]->GetBinContent(yBinMin);

	while (sum < 0.954499736*nToys && yBinMin > 0 && yBinMax <= fluctHist_category_yield[xBin]->GetNbinsX()) {

	  double contMin = fluctHist_category_yield[xBin]->GetBinContent(yBinMin-1);
	  double contMax = fluctHist_category_yield[xBin]->GetBinContent(yBinMax+1);

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

	if (yBinMin == 0  && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to underflow for bin " << xBin+1 << "!" << std::endl;
	if (yBinMax == (fluctHist_category_yield[xBin]->GetNbinsX()+1) && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	double yMin = fluctHist_category_yield[xBin]->GetBinCenter(yBinMin);
	double yMax = fluctHist_category_yield[xBin]->GetBinCenter(yBinMax);
	double yAve = (yMin+yMax)/2;

	errHist_category_yield_2sig->SetBinContent(xBin+1,yAve);
	errHist_category_yield_2sig->SetBinError(xBin+1,yMax-yAve);
      }

      TH1 *ratioErr_category_yield_2sig = (TH1 *)ratio_category_yield->Clone("ratioErr_category_yield_2sig_"+fitLabel);
      ratioErr_category_yield_2sig->Reset();
      ratioErr_category_yield_2sig->SetFillColor(kYellow);
      ratioErr_category_yield_2sig->SetLineColor(0);
      ratioErr_category_yield_2sig->SetLineWidth(0);
      ratioErr_category_yield_2sig->SetMarkerStyle(0);
      ratioErr_category_yield_2sig->SetMarkerColor(0);
      ratioErr_category_yield_2sig->SetMarkerSize(0.);

      for (int iBin = 1; iBin <= numCats; iBin++) {
	double mcCentral = h_bkg_category_yield->GetBinContent(iBin);
	double mcUp = errHist_category_yield_2sig->GetBinContent(iBin)+errHist_category_yield_2sig->GetBinError(iBin);
	double mcDown = errHist_category_yield_2sig->GetBinContent(iBin)-errHist_category_yield_2sig->GetBinError(iBin);

	double ratioUp = mcUp/mcCentral;
	double ratioDown = mcDown/mcCentral;
	double ratioAve = (ratioUp+ratioDown)/2;

	ratioErr_category_yield_2sig->SetBinContent(iBin,ratioAve);
	ratioErr_category_yield_2sig->SetBinError(iBin, ratioUp-ratioAve);
      }


      if( debug ) {
	for( int xBin = 1; xBin <= numCats; xBin++ ) {
	  std::cout << " xBin = " << xBin << ",\t dataHist = " << h_data_category_yield->GetBinContent(xBin) << ",\t fitHist = " << h_bkg_category_yield->GetBinContent(xBin) << ",\t errHist_1sig = " << errHist_category_yield_1sig->GetBinContent(xBin) << " +/- " << errHist_category_yield_1sig->GetBinError(xBin) << ",\t ratio = " << ratio_category_yield->GetBinContent(xBin) << ",\t ratioErr_1sig = " << ratioErr_category_yield_1sig->GetBinContent(xBin) << " +/- " << ratioErr_category_yield_1sig->GetBinError(xBin) << std::endl;
	}
      }
 

      //Hack to get it plotted with ratio plot
      TCanvas* myC = new TCanvas("myC", "myC", 800,700);
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

      ratioErr_category_yield_1sig->SetMinimum(ratioMin);
      ratioErr_category_yield_1sig->SetMaximum(ratioMax);
      ratioErr_category_yield_1sig->GetYaxis()->SetNdivisions(50000+404);
      ratioErr_category_yield_1sig->GetYaxis()->CenterTitle();
      ratioErr_category_yield_1sig->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratioErr_category_yield_1sig->GetXaxis()->SetLabelSize(0.1); //make x label bigger
      ratioErr_category_yield_1sig->GetXaxis()->SetTitleOffset(1.1); //ccn
      //ratioErr_category_yield_1sig->GetXaxis()->SetTitleOffset(0.9); //SBOUTLE
      //ratioErr_category_yield_1sig->GetXaxis()->SetTitleOffset(1.1);
      ratioErr_category_yield_1sig->GetXaxis()->SetTitle(h_data_category_yield->GetXaxis()->GetTitle()); //make y label bigger
      ratioErr_category_yield_1sig->GetXaxis()->SetLabelSize(0.12);
      ratioErr_category_yield_1sig->GetXaxis()->SetLabelOffset(0.02); //SBOUTLE
      //ratioErr_category_yield_1sig->GetXaxis()->SetLabelOffset(0.04);
      ratioErr_category_yield_1sig->GetXaxis()->SetTitleSize(0.13); //SBOUTLE
      //ratioErr_category_yield_1sig->GetXaxis()->SetTitleSize(0.12);
      ratioErr_category_yield_1sig->GetYaxis()->SetTitle("Data/MC");
      ratioErr_category_yield_1sig->GetYaxis()->SetTitleSize(0.1);
      ratioErr_category_yield_1sig->GetYaxis()->SetTitleOffset(.45);
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

      //double scale_ttH = ( fitLabel.EqualTo("postFitS") && !blind ) ? 1. : 30.; // fitHist->Integral() / h_ttH->Integral()
      double scale_ttH = 30.; // fitHist->Integral() / h_ttH->Integral()

 
      legend->AddEntry(h_bkg[bin_ttjets_other],Form("%s (%.1f)",label[bin_ttjets_other].Data(),h_bkg[bin_ttjets_other]->Integral()),"f");
      legend->AddEntry(h_bkg[bin_ttjets_ccbar],Form("%s (%.1f)",label[bin_ttjets_ccbar].Data(),h_bkg[bin_ttjets_ccbar]->Integral()),"f");
      legend->AddEntry(h_bkg[bin_ttjets_b],Form("%s (%.1f)",label[bin_ttjets_b].Data(),h_bkg[bin_ttjets_b]->Integral()),"f");
      legend->AddEntry(h_bkg[bin_ttjets_bbbar],Form("%s (%.1f)",label[bin_ttjets_bbbar].Data(),h_bkg[bin_ttjets_bbbar]->Integral()),"f");
      legend->AddEntry(h_bkg[bin_singlet],Form("%s (%.1f)",label[bin_singlet].Data(),h_bkg[bin_singlet]->Integral()),"f");
      
      legend->AddEntry(h_ttbarV_bkg,Form("%s (%.1f)",label[bin_ttV].Data(),h_ttbarV_bkg->Integral()),"f");
      legend->AddEntry(h_ewk_bkg,Form("%s (%.1f)",label[bin_ewk].Data(),h_ewk_bkg->Integral()),"f");

      legend->AddEntry(h_data_category_yield,Form("%s (%.1f)",label[bin_data].Data(),h_data_category_yield->Integral()),"lpe");

      legend->AddEntry(errHist_category_yield_1sig,Form("Sum MC (%.1f)",h_data_category_yield->Integral()),"f");

      //if( fitLabel.EqualTo("postFitS") && !blind ) legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f)",h_ttH->Integral()),"l");
      //else                                         legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f#times%.0f)",h_ttH->Integral(),scale_ttH),"l");
      legend->AddEntry(h_ttH,Form("t#bar{t}H125.6 (%.1f#times%.0f)",h_ttH->Integral(),scale_ttH),"l");



      //TLegend *legend_ratio = new TLegend(0.14,0.87,0.34,0.97);
      TLegend *legend_ratio = new TLegend(0.14,0.77,0.48,0.97);//ccn

      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetShadowColor(kWhite);
      legend_ratio->SetTextFont(42);
      //legend_ratio->SetTextSize(0.055);
      legend_ratio->SetTextSize(0.10);//ccn

      if( !fitLabel.EqualTo("preFit") ){
	legend_ratio->SetNColumns(2);

	legend_ratio->AddEntry(ratioErr_category_yield_1sig,"Fit #pm1#sigma","f");
	legend_ratio->AddEntry(ratioErr_category_yield_2sig,"Fit #pm2#sigma","f");
      }
      else {
	legend_ratio->AddEntry(ratioErr_category_yield_1sig,"Stat+Syst Uncertainty","f");
      }

      ratioErr_category_yield_1sig->SetTitle(";channel and category");
      for( int iBin=0; iBin<ratioErr_category_yield_1sig->GetNbinsX(); iBin++ ) ratioErr_category_yield_1sig->GetXaxis()->SetBinLabel(iBin+1,axis_labels[iBin]);

      h_data_category_yield->SetTitle(";;Events");
      h_data_category_yield->GetYaxis()->SetTitleSize(0.05);
      h_data_category_yield->GetYaxis()->SetTitleOffset(1.1); //SBOUTLE
      //h_data_category_yield->GetYaxis()->SetTitleOffset(1.);

      double xmin = h_data_category_yield->GetBinLowEdge(1);
      double xmax = h_data_category_yield->GetBinLowEdge(numCats) + h_data_category_yield->GetBinWidth(numCats);



      myC->cd(1);
      myC->GetPad(1)->SetLogy(1);

      int maxBin_data = h_data_category_yield->GetMaximumBin();
      int maxBin_bkg  = errHist_category_yield_1sig->GetMaximumBin();

      double max_data = h_data_category_yield->GetBinContent(maxBin_data) + h_data_category_yield->GetBinError(maxBin_data);
      double max_bkg  = errHist_category_yield_1sig->GetBinContent(maxBin_bkg) + errHist_category_yield_1sig->GetBinError(maxBin_bkg);

      double histMax = TMath::Max( max_data, max_bkg );

      if( useLegend ) histMax = 10. * histMax;
      else            histMax = 10. * histMax;

      h_data_category_yield->SetMaximum(histMax);
      errHist_category_yield_1sig->SetMaximum(histMax);

      double chi2 = h_data_category_yield->Chi2Test(h_bkg_category_yield_blind,"UWP");

      TString str_chi2 = ( !fitLabel.EqualTo("postFitS") ) ? Form("B-only p-value (#chi^{2}) = %.3f",chi2) : Form("S+B p-value (#chi^{2}) = %.3f",chi2);
      TLatex Chi2Latex(0.67, 0.88, str_chi2);
      Chi2Latex.SetNDC();
      Chi2Latex.SetTextFont(42);
      //Chi2Latex.SetTextSize(0.055);
      Chi2Latex.SetTextSize(0.09);//ccn


      //ccn: Fig 28?
      TString selectioninfo = "ttH(bb,#tau#tau): LJ + DIL + TAU";//"Lepton + #geq6 jets + 2 tags";
      TLatex SELECTIONInfoLatex(0.60, 0.84, selectioninfo);
      SELECTIONInfoLatex.SetNDC();
      SELECTIONInfoLatex.SetTextFont(42);
      //SELECTIONInfoLatex.SetTextSize(0.035);
      SELECTIONInfoLatex.SetTextSize(0.05);//ccn


      //if( fitLabel.EqualTo("postFitS") && !blind ) hs->Add(h_ttH);
      //else                                         h_ttH->Scale( scale_ttH );
      h_ttH->Scale( scale_ttH );





      h_data_category_yield->SetMinimum(3);
      h_data_category_yield->GetYaxis()->SetLabelSize(0.05);//ccn - make y-axis labels bigger
      h_data_category_yield->Draw("e1");
      hs->Draw("histsame");
      errHist_category_yield_1sig->Draw("e2same");
      //if( !fitLabel.EqualTo("postFitS") || blind ) h_ttH->Draw("histsame");
      h_ttH->Draw("histsame");
      h_data_category_yield->Draw("e1SAME");

      if( useLegend ) legend->Draw();


      int myline_width = 3;
      int myline_style = 1;

      TLine* myLine_channel1;
      myLine_channel1 = new TLine(7, 0.000, 7, 2e4);
      myLine_channel1->SetLineWidth(myline_width);
      myLine_channel1->SetLineStyle(myline_style);
      myLine_channel1->Draw();

      TLine* myLine_channel2;
      myLine_channel2 = new TLine(10, 0.000, 10, 2e4);
      myLine_channel2->SetLineWidth(myline_width);
      myLine_channel2->SetLineStyle(myline_style);
      myLine_channel2->Draw();


      //ccn Fig 28?
      //TLatex CMSInfoLatex_wide(0.57, 0.91, cmsinfo.c_str());
      TLatex CMSInfoLatex_wide(0.14, 0.91, cmsinfo.c_str());//ccn
      CMSInfoLatex_wide.SetNDC();
      CMSInfoLatex_wide.SetTextFont(42);
      //CMSInfoLatex_wide.SetTextSize(0.035);
      CMSInfoLatex_wide.SetTextSize(0.055);//ccn


      CMSInfoLatex_wide.Draw();
      SELECTIONInfoLatex.Draw();
      FitinfoLatex.Draw();


      myC->cd(2);
      ratioErr_category_yield_1sig->Draw("e2");
      if( !fitLabel.EqualTo("preFit") ){
	ratioErr_category_yield_2sig->Draw("e2same");
	ratioErr_category_yield_1sig->Draw("e2same");
      }
      ratio_category_yield->Draw("e1SAME");


      double use_ratio[ratio_category_yield->GetNbinsX()];
      for( int iBin=0; iBin<ratio_category_yield->GetNbinsX(); iBin++ ) use_ratio[iBin] = ratio_category_yield->GetBinContent(iBin+1);

      TLatex *bin_ratio[ratio_category_yield->GetNbinsX()];
      for( int iBin=0; iBin<ratio_category_yield->GetNbinsX(); iBin++ ){
	//bin_ratio[iBin] = new TLatex(0.135+0.093*iBin,0.33,Form("%1.2f", use_ratio[iBin]) );
	bin_ratio[iBin] = new TLatex(0.123+0.0522*iBin,0.33,Form("%1.2f", use_ratio[iBin]) );
	bin_ratio[iBin]->SetNDC();
	bin_ratio[iBin]->SetTextFont(42);
	bin_ratio[iBin]->SetTextSize(0.07);
	//bin_ratio[iBin]->Draw();
      }

      myC->GetPad(1)->RedrawAxis();
      myC->GetPad(2)->RedrawAxis();


      myLine_channel1->Draw();
      myLine_channel2->Draw();

      TLine* myLine;
      myLine = new TLine(xmin, 1.000, xmax, 1.000);
      myLine->Draw();
      legend_ratio->Draw();

      Chi2Latex.Draw();

      myC->SaveAs(imageDir+"/dataToMC_"+era+"_category_yield_"+fitLabel+".png");
      myC->SaveAs(imageDir+"/dataToMC_"+era+"_category_yield_"+fitLabel+".pdf");

      delete hs;
      delete myLine;
      delete myC;


      if( debug ){

	TLine* l_nom = new TLine( 0,0,1,1 );
	TLine* l_ave = new TLine( 0,0,1,1 );
	TLine* l_ave_p1 = new TLine( 0,0,1,1 );
	TLine* l_ave_m1 = new TLine( 0,0,1,1 );
	TLine* l_ave_p2 = new TLine( 0,0,1,1 );
	TLine* l_ave_m2 = new TLine( 0,0,1,1 );

	l_nom->SetLineColor(kBlack);
	l_ave->SetLineColor(kBlue);
	l_ave_p1->SetLineColor(kGreen+1);
	l_ave_m1->SetLineColor(kGreen+1);
	l_ave_p2->SetLineColor(kRed+1);
	l_ave_m2->SetLineColor(kRed+1);

	l_nom->SetLineWidth(2);
	l_ave->SetLineWidth(2);
	l_ave_p1->SetLineWidth(2);
	l_ave_m1->SetLineWidth(2);
	l_ave_p2->SetLineWidth(2);
	l_ave_m2->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1","",900,800);
	for( int iBin=0; iBin<numCats; iBin++ ){

	  TString s_bin = Form("%d",iBin+1);
	  TString proj_name = "projection_bin"+s_bin+"_"+fitLabel;

	  TH1D* h_proj = (TH1D*)fluctHist_category_yield[iBin]->Clone(proj_name);
	  int rebin = 1000;
	  if( !fitLabel.EqualTo("preFit")) rebin = 100;
	  h_proj->Rebin(rebin);

	  int nProjBins = h_proj->GetNbinsX();
	  std::cout << "  h_proj nbins = " << nProjBins << std::endl;

	  double ave = errHist_category_yield_1sig->GetBinContent(iBin+1);
	  double aveErr = errHist_category_yield_1sig->GetBinError(iBin+1);
	  double nom = h_bkg_category_yield->GetBinContent(iBin+1);

	  double ave_2s = errHist_category_yield_2sig->GetBinContent(iBin+1);
	  double aveErr_2s = errHist_category_yield_2sig->GetBinError(iBin+1);

	  double scaleBand = 4;
	  double xminProj = std::max( h_proj->GetXaxis()->GetXmin(), ave-scaleBand*aveErr );
	  double xmaxProj = std::min( h_proj->GetXaxis()->GetXmax(), ave+scaleBand*aveErr );

	  double ymaxProj = 1.05 * h_proj->GetMaximum();

	  double integral_total = h_proj->Integral();
	  double integral_68 = h_proj->Integral( h_proj->FindBin(ave-aveErr), h_proj->FindBin(ave+aveErr) );
	  double integral_95 = h_proj->Integral( h_proj->FindBin(ave_2s-aveErr_2s), h_proj->FindBin(ave_2s+aveErr_2s) );

	  printf("\t bin = %d,\t integral_total = %.1f,\t integral_68 = %.1f,\t integral_95 = %.1f,\t 68/total = %.3f,\t 95/total = %.3f \n",iBin+1,integral_total, integral_68, integral_95,integral_68/integral_total,integral_95/integral_total);

	  h_proj->GetXaxis()->SetRangeUser(xminProj,xmaxProj);
	  h_proj->GetYaxis()->SetRangeUser(0,ymaxProj);

	  h_proj->SetLineColor(kBlack);
	  h_proj->Draw("hist");

	  l_nom->DrawLine(nom,0,nom,ymaxProj);
	  l_ave->DrawLine(ave,0,ave,ymaxProj);
	  l_ave_p1->DrawLine(ave+aveErr,0,ave+aveErr,ymaxProj);
	  l_ave_m1->DrawLine(ave-aveErr,0,ave-aveErr,ymaxProj);
	  if( !fitLabel.EqualTo("preFit") ){
	    l_ave_p2->DrawLine(ave_2s+aveErr_2s,0,ave_2s+aveErr_2s,ymaxProj);
	    l_ave_m2->DrawLine(ave_2s-aveErr_2s,0,ave_2s-aveErr_2s,ymaxProj);
	  }

	  printf("  bin = %d,\t nom = %.2f,\t ave = %.2f,\t ave - err = %.2f,\t ave+err = %.2f \n",iBin+1,nom,ave,ave-aveErr,ave+aveErr );

	  c1->Print(debugDir+"/projection_category_yield_"+fitLabel+"_bin"+s_bin+".png");
	}

	delete c1;
	delete l_nom;
	delete l_ave;
	delete l_ave_p1;
	delete l_ave_m1;
	delete l_ave_p2;
	delete l_ave_m2;
      }

    }





    if( plotBDTFits ){
      TIter nextCat3(&categories);
      LabelInfo *category3 = 0;
      int iCat3=0;

      std::vector<TString> temp_BDT_cat_names;
      std::vector<double>  temp_BDT_chi2_value;

      while ((category3 = (LabelInfo *)nextCat3())) {
	TString catLabel = category3->label;
	TString catName  = category3->name;

	bool isTauCat = ( catName.Contains("TTL_") );
	std::cout << "\n\t\t " << fitLabel << "\t category = " << catLabel << std::endl;
      
	TString selectioninfo = catLabel;//"Lepton + #geq6 jets + 2 tags";
	TLatex SELECTIONInfoLatex(0.45, 0.84, selectioninfo); //SBOUTLE
	//TLatex SELECTIONInfoLatex(0.12, 0.91, selectioninfo);
	SELECTIONInfoLatex.SetNDC();
	SELECTIONInfoLatex.SetTextFont(42);
	//SELECTIONInfoLatex.SetTextSize(0.035);
	SELECTIONInfoLatex.SetTextSize(0.05); //SBOUTLE

	TString dataCatName = catName;
	if( dataCatName.Contains("ch1_") ) dataCatName.ReplaceAll("ch1_","");
	if( dataCatName.Contains("hbb_7TeV_") ) dataCatName.ReplaceAll("hbb_7TeV_","");
	if( dataCatName.Contains("hbb_8TeV_") ) dataCatName.ReplaceAll("hbb_8TeV_","");
	if( dataCatName.Contains("htt_8TeV_") ) dataCatName.ReplaceAll("htt_8TeV_","");

	TH1 *dataHist = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get("data_obs_CFMlpANN_"+dataCatName) : (TH1 *)dataFile->Get("data_obs_MVA_"+dataCatName);
	dataHist->SetLineColor(1);
	dataHist->SetLineWidth(4); // SBOUTLE
	//dataHist->SetLineWidth(2);
	dataHist->SetMarkerStyle(20);
	dataHist->SetMarkerSize(1.2); //SBOUTLE
	//dataHist->SetMarkerSize(0.75);


	TH1 *sigHist_original = ( catName.BeginsWith("hbb_7TeV_") ) ? (TH1 *)dataFile7TeV->Get("ttH125.6_CFMlpANN_"+dataCatName)->Clone("ttH125.6_original_"+catName+"_"+fitLabel) : (TH1 *)dataFile->Get("ttH125.6_MVA_"+dataCatName)->Clone("ttH125.6_original_"+catName+"_"+fitLabel);

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


	//////

	TH1* sigTempHist = NULL;
	if( catName.BeginsWith("hbb_7TeV_") ){
	  TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";
	  RooAbsReal *fitTTH = w->function(normName_ttH);

	  sigTempHist = w->pdf("shapeSig_"+catName+"_ttH_morph")->createHistogram("CMS_th1x");
	  sigTempHist->Scale(fitTTH->getVal());
	}
	else {
	  // hbb
	  TString normName_ttH_bb = "n_exp_final_bin"+catName+"_proc_ttH_hbb";
	  RooAbsReal *fitTTH_bb = w->function(normName_ttH_bb);

	  //TH1* sigTempHist_bb = w->pdf("shapeSig_"+catName+"_ttH_hbb_morph")->createHistogram("CMS_th1x");
	  TH1* sigTempHist_bb = NULL;
	  if( w->pdf("shapeSig_"+catName+"_ttH_hbb_morph")!=NULL ){
	    sigTempHist_bb = w->pdf("shapeSig_"+catName+"_ttH_hbb_morph")->createHistogram("CMS_th1x");
	    sigTempHist_bb->Scale(fitTTH_bb->getVal());
	  }

	  // hcc
	  TString normName_ttH_cc = "n_exp_final_bin"+catName+"_proc_ttH_hcc";
	  RooAbsReal *fitTTH_cc = w->function(normName_ttH_cc);

	  TH1* sigTempHist_cc = NULL;
	  if( fitTTH_cc!=NULL ){
	    sigTempHist_cc = w->pdf("shapeSig_"+catName+"_ttH_hcc_morph")->createHistogram("CMS_th1x");
	    sigTempHist_cc->Scale(fitTTH_cc->getVal());
	  }

	  // hgg
	  TString normName_ttH_gg = "n_exp_final_bin"+catName+"_proc_ttH_hgg";
	  RooAbsReal *fitTTH_gg = w->function(normName_ttH_gg);

	  TH1* sigTempHist_gg = NULL;
	  if( fitTTH_gg!=NULL ){
	    sigTempHist_gg = w->pdf("shapeSig_"+catName+"_ttH_hgg_morph")->createHistogram("CMS_th1x");
	    sigTempHist_gg->Scale(fitTTH_gg->getVal());
	  }

	  // hgluglu
	  TString normName_ttH_gluglu = "n_exp_final_bin"+catName+"_proc_ttH_hgluglu";
	  RooAbsReal *fitTTH_gluglu = w->function(normName_ttH_gluglu);

	  TH1* sigTempHist_gluglu = NULL;
	  if( fitTTH_gluglu!=NULL ){
	    sigTempHist_gluglu = w->pdf("shapeSig_"+catName+"_ttH_hgluglu_morph")->createHistogram("CMS_th1x");
	    sigTempHist_gluglu->Scale(fitTTH_gluglu->getVal());
	  }

	  // htt
	  TString normName_ttH_tt = "n_exp_final_bin"+catName+"_proc_ttH_htt";
	  RooAbsReal *fitTTH_tt = w->function(normName_ttH_tt);

	  TH1* sigTempHist_tt = NULL;
	  if( fitTTH_tt !=NULL ){
	    sigTempHist_tt = w->pdf("shapeSig_"+catName+"_ttH_htt_morph")->createHistogram("CMS_th1x");
	    sigTempHist_tt->Scale(fitTTH_tt->getVal());
	  }

	  // hww
	  TString normName_ttH_ww = "n_exp_final_bin"+catName+"_proc_ttH_hww";
	  RooAbsReal *fitTTH_ww = w->function(normName_ttH_ww);

	  TH1* sigTempHist_ww = NULL;
	  if( fitTTH_ww!=NULL ){
	    sigTempHist_ww = w->pdf("shapeSig_"+catName+"_ttH_hww_morph")->createHistogram("CMS_th1x");
	    sigTempHist_ww->Scale(fitTTH_ww->getVal());
	  }

	  // hzg
	  TString normName_ttH_zg = "n_exp_final_bin"+catName+"_proc_ttH_hzg";
	  RooAbsReal *fitTTH_zg = w->function(normName_ttH_zg);

	  TH1* sigTempHist_zg = NULL;
	  if( fitTTH_zg!=NULL ){
	    sigTempHist_zg = w->pdf("shapeSig_"+catName+"_ttH_hzg_morph")->createHistogram("CMS_th1x");
	    sigTempHist_zg->Scale(fitTTH_zg->getVal());
	  }

	  // hzz
	  TString normName_ttH_zz = "n_exp_final_bin"+catName+"_proc_ttH_hzz";
	  RooAbsReal *fitTTH_zz = w->function(normName_ttH_zz);

	  TH1* sigTempHist_zz = NULL;
	  if( fitTTH_zz!=NULL ){
	    sigTempHist_zz = w->pdf("shapeSig_"+catName+"_ttH_hzz_morph")->createHistogram("CMS_th1x");
	    sigTempHist_zz->Scale(fitTTH_zz->getVal());
	  }

	  sigTempHist = (TH1*)sigTempHist_bb->Clone("CMS_th1x");
	  if( sigTempHist_cc!=NULL ) sigTempHist->Add(sigTempHist_cc);
	  if( sigTempHist_gg!=NULL ) sigTempHist->Add(sigTempHist_gg);
	  if( sigTempHist_gluglu!=NULL ) sigTempHist->Add(sigTempHist_gluglu);
	  if( sigTempHist_tt!=NULL ) sigTempHist->Add(sigTempHist_tt);
	  if( sigTempHist_ww!=NULL ) sigTempHist->Add(sigTempHist_ww);
	  if( sigTempHist_zg!=NULL ) sigTempHist->Add(sigTempHist_zg);
	  if( sigTempHist_zz!=NULL ) sigTempHist->Add(sigTempHist_zz);
	}
	//////

	//////
	//////


	TH1* h_ttH = (TH1 *)dataHist->Clone("pdf_bin"+catName+"_ttH_"+fitLabel+"Clone");
	h_ttH->Reset();
	for (int iBin = 0; iBin < numBins[iCat3]; iBin++) {
	  h_ttH->SetBinContent(iBin+1,sigTempHist->GetBinContent(iBin+1));
	}

	double blindSscale = sigHist_original->Integral() / h_ttH->Integral();

	if( blind ){
	  h_ttH->Scale( blindSscale );
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
	  hs->Add(h_ewk_bkg);
	  //hs->Add(h_bkg[bin_diboson]);
	  //hs->Add(h_bkg[bin_wjets]);
	  //hs->Add(h_bkg[bin_zjets]);
	  hs->Add(h_ttbarV_bkg);
	  hs->Add(h_bkg[bin_singlet]);
	  hs->Add(h_bkg[bin_ttjets_other]);
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

	  if (yBinMin == 0 && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to underflow for bin " << xBin+1 << "!" << std::endl;
	  if (yBinMax == (fluctHist[iCat3][xBin]->GetNbinsX()+1) && !fitLabel.EqualTo("preFit") ) std::cerr << "WARNING (2sig): error extends to overflow for bin " << xBin+1 << "!" << std::endl;

	  double yMin = fluctHist[iCat3][xBin]->GetBinCenter(yBinMin);
	  double yMax = fluctHist[iCat3][xBin]->GetBinCenter(yBinMax);
	  double yAve = (yMin+yMax)/2;

	  errHist_2sig->SetBinContent(xBin+1,yAve);
	  errHist_2sig->SetBinError(xBin+1,yMax-yAve);
	}



	//double chi2 = dataHistTemp->Chi2Test(fitHist,"UWP");
	//double chi2 = dataHist->Chi2Test(errHist_1sig,"WWP");

	TH1 *fitHist_blind = (TH1 *)fitHist->Clone("fitHist_blind"+catName+"_"+fitLabel);
	//std::cout << "\t Finished uncertainty band creation.  Making plots" << std::endl;
	if( blind ){
	  double threshold = 0.03;
	  bool aboveThreshold = false;
	  for (int iBin=0; iBin<numBins[iCat3]; iBin++) {
	    double bkgInBin = fitHist->GetBinContent(iBin+1);
	    double sigInBin = h_ttH->GetBinContent(iBin+1);
	    double SoverB = ( bkgInBin>0. ) ? sigInBin/bkgInBin : 0.;
	    aboveThreshold = ( aboveThreshold || (SoverB>threshold) );
	    if( aboveThreshold ){
	      dataHist->SetBinContent(iBin+1,0.);
	      fitHist_blind->SetBinContent(iBin+1,0.);
	    }
	  }
	}
	double chi2 = dataHist->Chi2Test(fitHist_blind,"UWP");


	temp_BDT_cat_names.push_back(catName);
	temp_BDT_chi2_value.push_back(chi2);


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
	ratioErr_1sig->GetYaxis()->CenterTitle();
	ratioErr_1sig->GetYaxis()->SetLabelSize(0.1); //make y label bigger
	ratioErr_1sig->GetXaxis()->SetLabelSize(0.1); //make y label bigger
	ratioErr_1sig->GetXaxis()->SetTitleOffset(0.9);
	//ratioErr_1sig->GetXaxis()->SetTitleOffset(1.1);
	ratioErr_1sig->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle()); //make y label bigger
	ratioErr_1sig->GetXaxis()->SetLabelSize(0.12);
	ratioErr_1sig->GetXaxis()->SetLabelOffset(0.02); //SBOUTLE
	//ratioErr_1sig->GetXaxis()->SetLabelOffset(0.04);
	ratioErr_1sig->GetXaxis()->SetTitleSize(0.14); //SBOUTLE
	//ratioErr_1sig->GetXaxis()->SetTitleSize(0.12);
	ratioErr_1sig->GetYaxis()->SetTitle("Data/MC");
	ratioErr_1sig->GetYaxis()->SetTitleSize(0.1);
	ratioErr_1sig->GetYaxis()->SetTitleOffset(.45);

	ratioErr_1sig->GetYaxis()->SetTitleFont(62);
	ratioErr_1sig->GetXaxis()->SetTitleFont(62);

	ratioErr_1sig->GetYaxis()->SetLabelFont(62);
	ratioErr_1sig->GetXaxis()->SetLabelFont(62);

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

	//double scale_ttH = ( fitLabel.EqualTo("postFitS") && !blind ) ? 1. : 30.; // fitHist->Integral() / h_ttH->Integral()
	//double scale_ttH = 30.; // fitHist->Integral() / h_ttH->Integral()
	double scale_ttH = ( fitLabel.EqualTo("postFitS") ) ? 10. : 30.; // fitHist->Integral() / h_ttH->Integral()


	if( isTauCat ){
	  //scale_ttH = ( fitLabel.EqualTo("postFitS") && !blind ) ? 1. : 100.;
	  //scale_ttH = 30.;

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

	//if( fitLabel.EqualTo("postFitS") && !blind ) legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f)",h_ttH->Integral()),"l");
	//else                                         legend->AddEntry(h_ttH,Form("t#bar{t}H125 (%.1f#times%.0f)",h_ttH->Integral(),scale_ttH),"l");
	legend->AddEntry(h_ttH,Form("t#bar{t}H125.6 (%.1f#times%.0f)",h_ttH->Integral(),scale_ttH),"l");



	//TLegend *legend_ratio = new TLegend(0.14,0.15,0.89,0.19);
	TLegend *legend_ratio = new TLegend(0.14,0.87,0.48,0.97);

	legend_ratio->SetFillColor(kWhite);
	legend_ratio->SetLineColor(kWhite);
	legend_ratio->SetShadowColor(kWhite);
	legend_ratio->SetTextFont(42);
	//legend_ratio->SetTextSize(0.055);
	legend_ratio->SetTextSize(0.10);//ccn

	if( !fitLabel.EqualTo("preFit") ){
	  legend_ratio->SetNColumns(2);

	  legend_ratio->AddEntry(ratioErr_1sig,"Fit #pm1#sigma","f");
	  legend_ratio->AddEntry(ratioErr_2sig,"Fit #pm2#sigma","f");
	}
	else {
	  legend_ratio->AddEntry(ratioErr_1sig,"Stat+Syst Uncertainty","f");
	}

	//ratioErr_1sig->SetTitle(";MVA output;Data/MC");
	ratioErr_1sig->SetTitle(";BDT output;Data/MC"); //SBOUTLE


	dataHist->SetTitle(";;Events");
	dataHist->GetYaxis()->SetTitleSize(0.05);
	dataHist->GetYaxis()->SetTitleOffset(1.1); //SBOUTLE
	//dataHist->GetYaxis()->SetTitleOffset(1.);

	double xmin = dataHist->GetBinLowEdge(1);
	double xmax = dataHist->GetBinLowEdge(numBins[iCat3]) + dataHist->GetBinWidth(numBins[iCat3]);

	dataHist->GetYaxis()->SetTitleFont(62);
	dataHist->GetXaxis()->SetTitleFont(62);

	dataHist->GetYaxis()->SetLabelFont(62);
	dataHist->GetXaxis()->SetLabelFont(62);


	//if( fitLabel.EqualTo("postFitS") && !blind ) hs->Add(h_ttH);
	//else                                         h_ttH->Scale( scale_ttH );
	h_ttH->Scale( scale_ttH );


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


	TString str_chi2 = Form("p-value (#chi^{2}) = %.3f",chi2);
	TLatex Chi2Latex(0.59, 0.91, str_chi2);
	Chi2Latex.SetNDC();
	Chi2Latex.SetTextFont(42);
	//Chi2Latex.SetTextSize(0.055);
	Chi2Latex.SetTextSize(0.09);//ccn



	dataHist->SetMinimum(1.5E-1);
	dataHist->Draw("e1");
	hs->Draw("histsame");
	errHist_1sig->Draw("e2same");
	//if( !fitLabel.EqualTo("postFitS") || blind ) h_ttH->Draw("histsame");
	h_ttH->Draw("histsame");
	dataHist->Draw("e1SAME");

	if( useLegend ) legend->Draw();

	if( catName.BeginsWith("hbb_7TeV_") ) CMSInfoLatex7TeV.Draw();
	else                             CMSInfoLatex.Draw();
	SELECTIONInfoLatex.Draw();
	//FitinfoLatex.Draw();


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

	Chi2Latex.Draw();

	TString saveName = catName;
	if( saveName.Contains("hbb_7TeV_") )      saveName.ReplaceAll("hbb_7TeV_","");
	else if( saveName.Contains("hbb_8TeV_") ) saveName.ReplaceAll("hbb_8TeV_","");
	else if( saveName.Contains("htt_8TeV_") ) saveName.ReplaceAll("htt_8TeV_","");

	if( catName.Contains("hbb_7TeV_") ) era = "7TeV";
	else                           era = "8TeV";

	myC->SaveAs(imageDir+"/dataToMC_"+era+"_finalMVA_"+fitLabel+"_"+saveName+".png");
	myC->SaveAs(imageDir+"/dataToMC_"+era+"_finalMVA_"+fitLabel+"_"+saveName+".pdf");



	if( debug ){
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

	  TCanvas *c1 = new TCanvas("c1","",900,800);
	  for( int iBin=0; iBin<numBins[iCat3]; iBin++ ){

	    TString s_bin = Form("%d",iBin+1);
	    TString proj_name = "projection_bin"+s_bin+"_"+fitLabel+"_"+catName;

	    TH1D* h_proj = (TH1D*)fluctHist[iCat3][iBin]->Clone(proj_name);
	    int rebin = 1000;
	    if( !fitLabel.EqualTo("preFit")) rebin = 100;
	    h_proj->Rebin(rebin);

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

	    c1->Print(debugDir+"/projection_"+fitLabel+"_"+catName+"_bin"+s_bin+".png");
	  }

	  delete c1;

	}



	iCat3++;

	// delete quantities defined for each category
	delete myC;
	delete legend;
	delete legend_ratio;
	delete myLine;
	delete hs;

      } // end loop over categories


      BDT_cat_names.push_back(temp_BDT_cat_names);
      BDT_chi2_value.push_back(temp_BDT_chi2_value);


    } // end conditional to loop over categories

    // delete quantities defined for each fit result
    for( int iPar=0; iPar<NumPars; iPar++ ) delete h_nuis[iPar];
    delete h_nuis;

  }// end loop over fitRes



  if( nuisance_names_noMCstat_postFitS.size()>0 && nuisance_names_noMCstat_postFitB.size()>0 ){
    int NumPulls_noMCstat_postFitS = int(nuisance_names_noMCstat_postFitS.size());
    int NumPulls_noMCstat_postFitB = int(nuisance_names_noMCstat_postFitB.size());
    TH2D* h_pulls_vert = new TH2D("h_pulls_vert",";Pull", 100, -2.5, 2.5, NumPulls_noMCstat_postFitS,0,NumPulls_noMCstat_postFitS);

    double use_y_postFitS[NumPulls_noMCstat_postFitS];
    double use_yErr_postFitS[NumPulls_noMCstat_postFitS];
    double pulls_noMCstat_postFitS[NumPulls_noMCstat_postFitS];
    //double constraints_noMCstat_postFitS[NumPulls_noMCstat_postFitS];
    double AsymErrorLo_noMCstat_postFitS[NumPulls_noMCstat_postFitS];
    double AsymErrorHi_noMCstat_postFitS[NumPulls_noMCstat_postFitS];

    double use_y_postFitB[NumPulls_noMCstat_postFitS];
    double use_yErr_postFitB[NumPulls_noMCstat_postFitS];
    double pulls_noMCstat_postFitB[NumPulls_noMCstat_postFitS];
    //double constraints_noMCstat_postFitB[NumPulls_noMCstat_postFitS];
    double AsymErrorLo_noMCstat_postFitB[NumPulls_noMCstat_postFitS];
    double AsymErrorHi_noMCstat_postFitB[NumPulls_noMCstat_postFitS];

    for( int iPar=0; iPar<NumPulls_noMCstat_postFitS; iPar++ ){
      use_y_postFitS[iPar] = iPar + 0.35;
      use_yErr_postFitS[iPar] = 0.;

      pulls_noMCstat_postFitS[iPar] = nuisance_pulls_noMCstat_postFitS[iPar];
      //constraints_noMCstat_postFitS[iPar] = nuisance_constraints_noMCstat_postFitS[iPar];

      double hiErr_postFitS = fabs(nuisance_AsymErrorHi_noMCstat_postFitS[iPar]);
      double loErr_postFitS = fabs(nuisance_AsymErrorLo_noMCstat_postFitS[iPar]);

      double maxError_postFitS = std::max<double>(std::max<double>(hiErr_postFitS, loErr_postFitS), nuisance_constraints_noMCstat_postFitS[iPar]);

      if( fabs(hiErr_postFitS) < 0.001*maxError_postFitS ) hiErr_postFitS = nuisance_constraints_noMCstat_postFitS[iPar];
      if( fabs(loErr_postFitS) < 0.001*maxError_postFitS ) loErr_postFitS = nuisance_constraints_noMCstat_postFitS[iPar];

      AsymErrorLo_noMCstat_postFitS[iPar] = fabs(loErr_postFitS);
      AsymErrorHi_noMCstat_postFitS[iPar] = fabs(hiErr_postFitS);

      h_pulls_vert->GetYaxis()->SetBinLabel(iPar+1,nuisance_names_noMCstat_postFitS[iPar].Data());


      for( int jPar=0; jPar<NumPulls_noMCstat_postFitB; jPar++ ){
	if( nuisance_names_noMCstat_postFitB[jPar].EqualTo(nuisance_names_noMCstat_postFitS[iPar]) ){
	  use_y_postFitB[jPar] = iPar + 0.65;
	  use_yErr_postFitB[jPar] = 0.;
	  pulls_noMCstat_postFitB[jPar] = nuisance_pulls_noMCstat_postFitB[jPar];
	  //constraints_noMCstat_postFitB[jPar] = nuisance_constraints_noMCstat_postFitB[jPar];

	  double hiErr_postFitB = fabs(nuisance_AsymErrorHi_noMCstat_postFitB[jPar]);
	  double loErr_postFitB = fabs(nuisance_AsymErrorLo_noMCstat_postFitB[jPar]);

	  double maxError_postFitB = std::max<double>(std::max<double>(hiErr_postFitB, loErr_postFitB), nuisance_constraints_noMCstat_postFitB[jPar]);

	  if( fabs(hiErr_postFitB) < 0.001*maxError_postFitB ) hiErr_postFitB = nuisance_constraints_noMCstat_postFitB[jPar];
	  if( fabs(loErr_postFitB) < 0.001*maxError_postFitB ) loErr_postFitB = nuisance_constraints_noMCstat_postFitB[jPar];

	  AsymErrorLo_noMCstat_postFitB[jPar] = fabs(loErr_postFitB);
	  AsymErrorHi_noMCstat_postFitB[jPar] = fabs(hiErr_postFitB);
	}
      }
    }

    //TGraphAsymmErrors *gr_postFitS = new TGraphAsymmErrors(NumPulls_noMCstat_postFitS,pulls_noMCstat_postFitS,use_y_postFitS,constraints_noMCstat_postFitS,constraints_noMCstat_postFitS,use_yErr_postFitS,use_yErr_postFitS); 
    TGraphAsymmErrors *gr_postFitS = new TGraphAsymmErrors(NumPulls_noMCstat_postFitS,pulls_noMCstat_postFitS,use_y_postFitS,AsymErrorLo_noMCstat_postFitS,AsymErrorHi_noMCstat_postFitS,use_yErr_postFitS,use_yErr_postFitS); 
    gr_postFitS->SetMarkerSize(1);
    gr_postFitS->SetMarkerStyle(20);
    gr_postFitS->SetMarkerColor(kBlack);
    gr_postFitS->SetLineColor(kBlack);


    //TGraphAsymmErrors *gr_postFitB = new TGraphAsymmErrors(NumPulls_noMCstat_postFitB,pulls_noMCstat_postFitB,use_y_postFitB,constraints_noMCstat_postFitB,constraints_noMCstat_postFitB,use_yErr_postFitB,use_yErr_postFitB); 
    TGraphAsymmErrors *gr_postFitB = new TGraphAsymmErrors(NumPulls_noMCstat_postFitB,pulls_noMCstat_postFitB,use_y_postFitB,AsymErrorLo_noMCstat_postFitB,AsymErrorHi_noMCstat_postFitB,use_yErr_postFitB,use_yErr_postFitB); 
    gr_postFitB->SetMarkerSize(1);
    gr_postFitB->SetMarkerStyle(21);
    gr_postFitB->SetMarkerColor(kRed);
    gr_postFitB->SetLineColor(kRed);

    if( debug ){
      TCanvas *c2 = new TCanvas("c2","",900,1400);
      h_pulls_vert->SetStats(0);
      h_pulls_vert->Draw();
      gr_postFitB->Draw("pe1same");
      gr_postFitS->Draw("pe1same");

      TLine* myLine1;
      myLine1 = new TLine( 0., h_pulls_vert->GetYaxis()->GetXmin(), 0., h_pulls_vert->GetYaxis()->GetXmax());
      myLine1->SetLineStyle(7);
      myLine1->Draw();

      TLegend *legend = new TLegend(0.3,0.955,0.94,0.99);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.03);

      legend->SetNColumns(2);

      legend->AddEntry(gr_postFitB,"Post-Fit B","pe1");
      legend->AddEntry(gr_postFitS,"Post-Fit S+B","pe1");

      legend->Draw();

      c2->GetPad(0)->SetLeftMargin(0.23);
      c2->GetPad(0)->SetTopMargin(0.05);
      c2->GetPad(0)->SetRightMargin(0.05);
      c2->Print(debugDir+"/pulls_vert_compare_postFitB_postFitS.png");

      delete c2;
      delete myLine1;
      delete legend;
    }
    delete gr_postFitS;
    delete gr_postFitB;
    delete h_pulls_vert;
  }


  int numUseFitRes = BDT_cat_names.size();

  if( numUseFitRes>0 && numUseFitRes<4 ){
    int numUseCats   = BDT_cat_names[0].size();

    if( numUseFitRes==1 ){
      std::cout << "\\begin{tabular}{|l|c|} \\hline" << std::endl;
      std::cout << "Category   & Pre-Fit \\\\  \\hline" << std::endl;
    }
    else if( numUseFitRes==2 ){
      std::cout << "\\begin{tabular}{|l|c|c|} \\hline" << std::endl;
      std::cout << "Category   & Pre-Fit & Post-Fit (B) \\\\  \\hline" << std::endl;
    }
    else if( numUseFitRes==3 ){
      std::cout << "\\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
      std::cout << "Category   & Pre-Fit & Post-Fit (B) & Post-Fit (S+B) \\\\  \\hline" << std::endl;
    }


    for( int iCat=0; iCat<numUseCats; iCat++ ){
      TString newname = BDT_cat_names[0][iCat];
      newname.ReplaceAll("_","\\_");
      std::cout << newname;

      for( int iRes=0; iRes<numUseFitRes; iRes++ ){
	printf("\t & \t %.3f", BDT_chi2_value[iRes][iCat]);
      }
      std::cout << " \\\\" << std::endl;
    }
    std::cout << "\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;

  } // end condition on at least one fit result


  file->Close();
  std::cout << "Done!" << std::endl;

}
