class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
};

void printNormsModSS (TString dataFileName = "", TString prefix_ttH = "ttH125") {

  gStyle->SetOptStat(0);

  //Suppress the printout from making plots
  gErrorIgnoreLevel = 2000;
 
  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libHiggsAnalysisCombinedLimit.so");
  
  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open("wsTest.root");
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the pdf with everything set to the postfit values
  TFile *postFitWSFile = TFile::Open("MaxLikelihoodFitResult.root");
  RooWorkspace *postFitW = (RooWorkspace *)postFitWSFile->Get("MaxLikelihoodFitResult");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open("mlfit.root");
  RooFitResult *preFitFR = (RooFitResult*)fitFile->Get("nuisances_prefit_res");  
  RooFitResult *postFitFR = (RooFitResult *)fitFile->Get("fit_b");
  RooFitResult *postFitFR_s = (RooFitResult *)fitFile->Get("fit_s");

  //Try making "snapshots" with preFit and postFit nuisance values
  w->saveSnapshot("prefit",RooArgSet(preFitFR->floatParsFinal()),true);
  w->saveSnapshot("postfit",RooArgSet(postFitFR->floatParsFinal()),true);
  w->saveSnapshot("postfit_s",RooArgSet(postFitFR_s->floatParsFinal()),true);

  //And if we want it, the data file
  TFile *dataFile = 0;
  if (dataFileName != "") {
    dataFile = TFile::Open(dataFileName);
  }


  //Now we need a list of categories and a list of processes
  TList processes;
  processes.Add(new LabelInfo("ttbar","$t\\bar{t}+$lf"));
  processes.Add(new LabelInfo("ttbarPlusBBbar","$t\\bar{t}+b\\bar{b}$"));
  processes.Add(new LabelInfo("ttbarPlusCCbar","$t\\bar{t}+c\\bar{c}$"));
  processes.Add(new LabelInfo("ttbarW","$t\\bar{t}W$"));
  processes.Add(new LabelInfo("ttbarZ","$t\\bar{t}Z$"));
  processes.Add(new LabelInfo("singlet","Single $t$"));
  processes.Add(new LabelInfo("wjets","$W+$jets"));
  processes.Add(new LabelInfo("zjets","$Z+$jets"));
  processes.Add(new LabelInfo("diboson","Diboson"));
  processes.Add(new LabelInfo("ttbarWW","$t\\bar{t}WW$"));
  processes.Add(new LabelInfo("ttbarttbar","$t\\bar{t}t\\bar{t}$"));
  processes.Add(new LabelInfo("WZZ","$WZZ$"));
  processes.Add(new LabelInfo("WWZ","$WWZ$"));
  processes.Add(new LabelInfo("WWW","$WWW$"));
  processes.Add(new LabelInfo("ZZZ","$WWW$"));

  int nRows = processes.GetEntries() + 1;
  if (dataFile) ++nRows;

  //These are the list of possible categories
  TList allCategories;
//   allCategories.Add(new LabelInfo("ljets_jge6_t2","LJ: $\\geq6$J,2T"));
//   allCategories.Add(new LabelInfo("ljets_j4_t3","LJ: 4J,3T"));
//   allCategories.Add(new LabelInfo("ljets_j5_t3","LJ: 5J,3T"));
//   allCategories.Add(new LabelInfo("ljets_jge6_t3","LJ: $\\geq6$J,3T"));
//   allCategories.Add(new LabelInfo("ljets_j4_t4","LJ: 4J,4T"));
//   allCategories.Add(new LabelInfo("ljets_j5_tge4","LJ: 5J,$\\geq3$T"));
//   allCategories.Add(new LabelInfo("ljets_jge6_tge4","LJ: $\\geq6$J,$\\geq4$T"));
//   allCategories.Add(new LabelInfo("e2je2t","DIL: 2J,2T"));
//   allCategories.Add(new LabelInfo("ge3t","DIL: $\\geq3$J,$\\geq3$T"));

  allCategories.Add(new LabelInfo("SS_e3je1t","DIL: $==3$J,$==1$T"));
  allCategories.Add(new LabelInfo("SS_ge4je1t","DIL: $\\geq4$J,$==1$T"));
  allCategories.Add(new LabelInfo("SS_e3jge2t","DIL: $==3$J,$\\geq2$T"));
  allCategories.Add(new LabelInfo("SS_ge4jge2t","DIL: $\\geq4$J,$\\geq2$T"));


  //Now look in the postfit workspace to see which ones were actually used
  TList categories;

  TIter nextCat(&allCategories);
  LabelInfo *c = 0;

  while (( c = (LabelInfo *)nextCat())) {
    TString name = c->name;
    if (postFitW->pdf("pdf_bin"+name+"_nuis")) {
      categories.Add(c);
      cout << "DEBUG: Including category " << name  << " in the output" << endl;
    } else {
      cout << "WARN: could not find fit result for category named " << name << endl;
    }
  }

  const int NumSamples = 9+3;
  const int NumCats = allCategories.GetEntries();

  cout << "DEBUG: NumCats = " << NumCats << endl;

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

    cout << "DEBUG: This is category " << catName << endl;

    TIter nextProcess(&processes);
    LabelInfo *process = 0;

    RooArgSet fitArgs;

    bool firstProc = true;

    double ttH_val_prefit = -1;
    double ttH_err_prefit = -1;
    double ttH_val_postfit = -1;
    double ttH_err_postfit = -1;

    double total_sumw2 = 0.;

    TString ttHHistName = prefix_ttH + "_CFMlpANN_" + catName;
    TH1 *ttHHist = (TH1 *)dataFile->Get(ttHHistName);
    double ttH_val_prefit = ttHHist->Integral();
    ttHHist->SetLineColor(kRed);
    ttHHist->SetLineWidth(2);

    int nbins_ttH = ttHHist->GetNbinsX();
    double sumw2_ttH = 0.;
    for( int b=0; b<nbins_ttH; b++ ) sumw2_ttH += ttHHist->GetBinError(b+1)*ttHHist->GetBinError(b+1);

    TString normName_ttH = "n_exp_final_bin"+catName+"_proc_ttH";

    RooAbsReal *fitTTH = w->function(normName_ttH);

    if( fitTTH ){
        w->loadSnapshot("prefit");        
        //ttH_val_prefit = fitTTH->getVal();
        ttH_err_prefit = fitTTH->getPropagatedError(*preFitFR);

        //Extract normalization and error on nomralization (after fit)
        w->loadSnapshot("postfit_s");        
        ttH_val_postfit = fitTTH->getVal();
        ttH_err_postfit = fitTTH->getPropagatedError(*postFitFR_s);
    }

    int iSample = 0;

    if( iCat==0 ) table_labels.push_back("$t\\bar{t}H(125)$");
    table_value[iSample][iCat] = ttH_val_prefit;
    table_value_syst_err[iSample][iCat] = ttH_err_prefit;
//     table_value[iSample][iCat] = ttH_val_postfit;
//     table_value_syst_err[iSample][iCat] = ttH_err_postfit;
    table_value_stat_err[iSample][iCat] = sqrt(sumw2_ttH);
    iSample++;

    while ((process = (LabelInfo *)nextProcess())) {

      TString procName = process->name;
      TString procLabel = process->label;

      cout << "DEBUG: Process " << procName
           << " cat = " << iCat 
           << endl;

      //Construct the name of the normalization variable.
      TString normName = "n_exp_final_bin"+catName+"_proc_"+procName;
      TString tempHistName = procName + "_CFMlpANN_"+catName;
      TH1 *procHist = (TH1 *)dataFile->Get(tempHistName);
    
      int nbins = procHist->GetNbinsX();
      double sumw2 = 0.;
      for( int b=0; b<nbins; b++ ) sumw2 += procHist->GetBinError(b+1)*procHist->GetBinError(b+1);

      //Extract normalization before fit
      double preFitNorm = -9e20, preFitErr = -9e10;
      double postFitNorm = -9e20, postFitErr = -9e10;
      double postFitNorm_s = -9e20, postFitErr_s = -9e10;

      RooAbsReal *fitN = w->function(normName);

      if (fitN) {

        fitArgs.add(*fitN);

        //Extract normalization and error on nomralization (before fit)
        w->loadSnapshot("prefit");
        preFitNorm = fitN->getVal();
        preFitErr = fitN->getPropagatedError(*preFitFR);

        cout << "    DEBUG: preFitNorm  for " << procName << " is " << preFitNorm
             << " error " << preFitErr 
             << " store in iSample = " << iSample
             << " iCat = " << iCat
             << endl;

        //Extract normalization and error on nomralization (after fit)
        w->loadSnapshot("postfit");        
        postFitNorm = fitN->getVal();
        postFitErr = fitN->getPropagatedError(*postFitFR);

        //Extract normalization and error on nomralization (after fit)
        w->loadSnapshot("postfit_s");        
        postFitNorm_s = fitN->getVal();
        postFitErr_s = fitN->getPropagatedError(*postFitFR_s);
      }

      total_sumw2 += sumw2;


      if (fitN) {

	if( iCat==0 ) table_labels.push_back(procLabel);
 	table_value[iSample][iCat] = preFitNorm;
 	table_value_syst_err[iSample][iCat] = preFitErr;
	// table_value[iSample][iCat] = postFitNorm;
	// table_value_syst_err[iSample][iCat] = postFitErr;
// 	table_value[iSample][iCat] = postFitNorm_s;
// 	table_value_syst_err[iSample][iCat] = postFitErr_s;
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

    w->loadSnapshot("prefit");
    double preFitTotalN = fitTotal.getVal();
    double preFitTotalErr = fitTotal.getPropagatedError(*preFitFR);
    
    w->loadSnapshot("postfit");
    double postFitTotalN = fitTotal.getVal();
    double postFitTotalErr = fitTotal.getPropagatedError(*postFitFR);

    w->loadSnapshot("postfit_s");
    double postFitTotalN_s = fitTotal.getVal();
    double postFitTotalErr_s = fitTotal.getPropagatedError(*postFitFR_s);


    if( iCat==0 ) table_labels.push_back("Total bkg");
    table_value[iSample][iCat] = preFitTotalN;
    table_value_syst_err[iSample][iCat] = preFitTotalErr;
    // table_value[iSample][iCat] = postFitTotalN;
    // table_value_syst_err[iSample][iCat] = postFitTotalErr;
//     table_value[iSample][iCat] = postFitTotalN_s;
//     table_value_syst_err[iSample][iCat] = postFitTotalErr_s;
    table_value_stat_err[iSample][iCat] = sqrt(total_sumw2);
    iSample++;


    if (dataFile) {
      //Get the data normalization
      TString dataHistName = "data_obs_CFMlpANN_" + catName;
      TH1 *dataHist = (TH1 *)dataFile->Get(dataHistName);
      double nData = dataHist->Integral();

      if( iCat==0 ) table_labels.push_back("Data");
      table_value[iSample][iCat] = nData;
      table_value_syst_err[iSample][iCat] = 0;
      table_value_stat_err[iSample][iCat] = 0;
      
    }


    ++iCat;

    
  } // while loop over categories

  //std::cout << " *************** " << std::endl;
  //std::cout << "    LJ yields    " << std::endl;
  // std::cout << " *************** " << std::endl;



  //std::cout << "    \\begin{tabular}{|l|c|c|c|c|c|c|c|} \\hline" << std::endl;
  //std::cout << "& $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets & 4 jets & 5 jets & $\\geq$6 jets \\\\" << std::endl;
  //std::cout << "& 2 tags & 3 tags & 3 tags & 3 tags & 4 tags & $\\geq$4 tags & $\\geq$4 tags \\\\ \\hline \\hline" << std::endl;

  std::cout << "\\documentclass{article}" << std::endl
            << "\\begin{document}" << std::endl;

  // this should be updated by numCats
  std::cout << "\\begin{tabular}{|l|c|c|c|c|}\n\\hline" << std::endl;

  
  TIter catLabels ( &allCategories);
  LabelInfo * eachCat =0;
  while (eachCat = (LabelInfo *) catLabels()) {
    std::cout << "& " << eachCat->label << "  ";
  }
  std::cout << " \\\\" << std::endl;
  
  for( int iSample=0; iSample<NumSamples; iSample++ ){
    if( table_labels[iSample].EqualTo("$t\\bar{t}Z$") || table_labels[iSample].EqualTo("$Z+$jets") ) continue;

    if( table_labels[iSample].Contains("Total") ) std::cout << " \\hline" << std::endl;

    if( table_labels[iSample].EqualTo("$t\\bar{t}W$") )  std::cout << "$t\\bar{t}V$";
    else if( table_labels[iSample].EqualTo("$W+$jets") ) std::cout << "$V+$jets";
    else                                                std::cout << table_labels[iSample];

    // only loop over the categories you have filled
    for( int iCat2=0; iCat2<NumCats; iCat2++ ){
      double val = table_value[iSample][iCat2];
      double syst_err2 = table_value_syst_err[iSample][iCat2] * table_value_syst_err[iSample][iCat2];
      double stat_err2 = table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2];
      if( table_labels[iSample].EqualTo("$t\\bar{t}W$")  ){
        val += table_value[iSample+1][iCat2];
        double sum_syst_err = table_value_syst_err[iSample][iCat2] + table_value_syst_err[iSample+1][iCat2];
        syst_err2 = sum_syst_err*sum_syst_err;
        
        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]) + (table_value_stat_err[iSample+1][iCat2] * table_value_stat_err[iSample+1][iCat2]);

        
      } else if (table_labels[iSample].EqualTo("$W+$jets")) {

        // do some extra work for W/Z+jets
        cout << "WJETS: --- This is Wjets ---" << endl
             << "WJETS: Norm = " << val << endl

        vector<int> otherSamples;
        otherSamples.push_back(iSample+1);
        otherSamples.push_back(iSample+3);
        otherSamples.push_back(iSample+4);
        otherSamples.push_back(iSample+5);
        otherSamples.push_back(iSample+6);
        otherSamples.push_back(iSample+7);
        otherSamples.push_back(iSample+8);
        double sum_syst_err = table_value_syst_err[iSample][iCat2];
        stat_err2 = (table_value_stat_err[iSample][iCat2] * table_value_stat_err[iSample][iCat2]);
        // and in Zjets
        for (int oSam = 0; oSam <otherSamples.size() ; oSam++){
          int nextSample = otherSamples[oSam];

          cout << "WJETS: Adding sample at " << nextSample
               << " with norm " << table_value[nextSample][iCat2] << endl;
          
          val += table_value[nextSample][iCat2];
          sum_syst_err += table_value_syst_err[nextSample][iCat2];
          
          stat_err2 = + (table_value_stat_err[nextSample][iCat2] * table_value_stat_err[nextSample][iCat2]);
          
        }
        syst_err2 = sum_syst_err*sum_syst_err;

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
  std::cout << "\\end{document}" << std::endl;


}
