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
#include "TObjArray.h"

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
#include <typeinfo>
#include <stdexcept>

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooStats/ModelConfig.h"
#include "RooAbsPdf.h"
#include "RooAbsCategoryLValue.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"

/*

Script to produce plots pre- and post-fit with uncertainty bands

to use, you need a datacard (e.g. datacard.dat), and a root file (e.g. file.root)

Do this once

combine -M MaxLikelihoodFit -m 125 --rMin -10 --rMax 10 --minos all datacard.dat

text2workspace.py -m 125 -D data_obs datacard.dat -b -o wsTest.root



Do this each time you want Asimov (i.e. fake data)

root -b -q head.C makeAsimovWithFitFromData.C+'("outfile.root")'

*/

class LabelInfo : public TObject {

public:
  TString name;
  TString label;

  LabelInfo(TString n, TString l): name(n), label(l) {}
  ClassDef(LabelInfo,1);
};

void copyAttributes(const RooAbsArg &from, RooAbsArg &to) {
  if (&from == &to) return;
  const std::set<std::string> attribs = from.attributes();
  if (!attribs.empty()) {
    for (std::set<std::string>::const_iterator it = attribs.begin(), ed = attribs.end(); it != ed; ++it) to.setAttribute(it->c_str());
  }
  const std::map<std::string, std::string> strattribs = from.stringAttributes();
  if (!strattribs.empty()) {
    for (std::map<std::string,std::string>::const_iterator it = strattribs.begin(), ed = strattribs.end(); it != ed; ++it) to.setStringAttribute(it->first.c_str(), it->second.c_str());
  }
}


//Stolen from HiggsAnalysis/CombinedLimit/src/utils.cc
RooAbsPdf *factorizePdf(const RooArgSet &observables, RooAbsPdf &pdf, RooArgList &constraints) {
  assert(&pdf);
  const std::type_info & id = typeid(pdf);
  if (id == typeid(RooProdPdf)) {
    //std::cout << " pdf is product pdf " << pdf.GetName() << std::endl;
    RooProdPdf *prod = dynamic_cast<RooProdPdf *>(&pdf);
    RooArgList newFactors; RooArgSet newOwned;
    RooArgList list(prod->pdfList());
    bool needNew = false;
    for (int i = 0, n = list.getSize(); i < n; ++i) {
      RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
      RooAbsPdf *newpdf = factorizePdf(observables, *pdfi, constraints);
      //std::cout << "    for " << pdfi->GetName() << "   newpdf  " << (newpdf == 0 ? "null" : (newpdf == pdfi ? "old" : "new"))  << std::endl;
      if (newpdf == 0) { needNew = true; continue; }
      if (newpdf != pdfi) { needNew = true; newOwned.add(*newpdf); }
      newFactors.add(*newpdf);
    }
    if (!needNew) { copyAttributes(pdf, *prod); return prod; }
    else if (newFactors.getSize() == 0) return 0;
    else if (newFactors.getSize() == 1) {
      RooAbsPdf *ret = (RooAbsPdf *) newFactors.first()->Clone(TString::Format("%s_obsOnly", pdf.GetName()));
      copyAttributes(pdf, *ret);
      return ret;
    }
    RooProdPdf *ret = new RooProdPdf(TString::Format("%s_obsOnly", pdf.GetName()), "", newFactors);
    ret->addOwnedComponents(newOwned);
    copyAttributes(pdf, *ret);
    return ret;
  } else if (id == typeid(RooSimultaneous)) {
    RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(&pdf);
    RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
    int nbins = cat->numBins((const char *)0);
    TObjArray factorizedPdfs(nbins); RooArgSet newOwned;
    bool needNew = false;
    for (int ic = 0, nc = nbins; ic < nc; ++ic) {
      cat->setBin(ic);
      RooAbsPdf *pdfi = sim->getPdf(cat->getLabel());
      RooAbsPdf *newpdf = factorizePdf(observables, *pdfi, constraints);
      factorizedPdfs[ic] = newpdf;
      if (newpdf == 0) {
        throw std::runtime_error(std::string("ERROR: channel ") + cat->getLabel() + " factorized to zero."); 
      }
      if (newpdf != pdfi) { needNew = true; newOwned.add(*newpdf); }
    }
    RooSimultaneous *ret = sim;
    if (needNew) {
      ret = new RooSimultaneous(TString::Format("%s_obsOnly", pdf.GetName()), "", (RooAbsCategoryLValue&) sim->indexCat());
      for (int ic = 0, nc = nbins; ic < nc; ++ic) {
        cat->setBin(ic);
        RooAbsPdf *newpdf = (RooAbsPdf *) factorizedPdfs[ic];
        if (newpdf) ret->addPdf(*newpdf, cat->getLabel());
      }
      ret->addOwnedComponents(newOwned);
    }
    delete cat;
    copyAttributes(pdf, *ret);
    return ret;
  } else if (pdf.dependsOn(observables)) {
    return &pdf;
  } else {
    if (!constraints.contains(pdf)) constraints.add(pdf);
    return 0;
  }
  
}



void makeAsimovWithFitFromData( TString asimovFileName = "", TString fitTypeName = "preFit", TString fitFileName = "mlfit.root", TString wsFileName = "wsTest.root", double scaleTTH=1.0 ) {

  //This file has the pdf with everything set to the values before the fit
  TFile *wsFile = TFile::Open( wsFileName );
  RooWorkspace *w = (RooWorkspace *)wsFile->Get("w");

  //This file has the actual fit results
  TFile *fitFile = TFile::Open( fitFileName );

  TFile *asimovFile = new TFile(asimovFileName, "RECREATE");
  TDirectory *outDir = asimovFile->mkdir("toys");

  TString fitLabel = fitTypeName;
  RooFitResult *fitRes;  
  
  if( fitTypeName.EqualTo("postFitB") )      fitRes = (RooFitResult*)fitFile->Get("fit_b");
  else if( fitTypeName.EqualTo("postFitS") ) fitRes = (RooFitResult*)fitFile->Get("fit_s");
  else                                       fitRes = (RooFitResult*)fitFile->Get("nuisances_prefit_res");

  RooArgSet fitargs = fitRes->floatParsFinal();
  if (!fitTypeName.EqualTo("postFitS")) {
    RooRealVar *r = w->var("r"); 
    fitargs.add(*r);  //Make sure we set the value of the signal strength
  }

  //r should be there now for every type
  RooRealVar *r = dynamic_cast<RooRealVar *>(fitargs.find("r"));
  if (r) {
    //Set the signal strength as desired.
    r->setVal(scaleTTH);
  } else {
    std::cerr << "Signal strength parameter r not found in fitargs." << std::endl;
    exit(2);
  }

  //Make snapshot!
  w->saveSnapshot("fitResult",fitargs,true);
  w->loadSnapshot("fitResult");

  asimovFile->cd("toys");
  
  //This is extracted from the combine code (various parts) for generating an Asimov dataset from a pdf.
  RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
  if (!mc) {
    std::cerr << "Error:  Could not load ModelConfig from workspace." << std::endl;
    exit(2);
  }

  RooAbsData *asimovData = 0;

  //Get the RooSimultaneous pdf
  RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(mc->GetPdf());
  RooArgSet *observables = simPdf->getObservables(*mc->GetObservables());

  RooAbsCategoryLValue *cat = const_cast<RooAbsCategoryLValue *>(&simPdf->indexCat());

  //Loop over each category and add data
  std::map<std::string,RooAbsData*> datasetPieces;

  RooArgList dummy;  
  RooRealVar *weightVar = new RooRealVar("_weight_","",1.0);
  for (int i = 0, n = cat->numBins((const char *)0); i < n; ++i) {  

    cat->setBin(i);

    RooAbsPdf *pdftemp = simPdf->getPdf(cat->getLabel());
    RooAbsPdf *pdf = factorizePdf(*observables, *pdftemp, dummy);    

    RooArgSet *obs = pdf->getObservables(*observables);

    RooAbsData *&data =  datasetPieces[cat->getLabel()]; delete data;
    
    RooArgList obsList(*obs);
    RooRealVar *x = (RooRealVar*)obsList.at(0);
    //Assume 1-D
    RooCmdArg ay = RooCmdArg::none();
    RooCmdArg az = RooCmdArg::none();

    TH1 *histoSpec = pdf->createHistogram("htemp", *x, ay, az); 
    histoSpec->SetDirectory(0);

    double expectedEvents = pdf->expectedEvents(*obs);
    histoSpec->Scale(expectedEvents/ histoSpec->Integral("width")); 
    RooArgSet obsPlusW(*obs); obsPlusW.add(*weightVar);
    data = new RooDataSet(TString::Format("%sData", pdf->GetName()), "", obsPlusW, weightVar->GetName());
    RooAbsArg::setDirtyInhibit(true); // don't propagate dirty flags while filling histograms 

    for (int ii = 1, nn = histoSpec->GetNbinsX(); ii <= nn; ++ii) {
      x->setVal(histoSpec->GetXaxis()->GetBinCenter(ii));
      double ww = histoSpec->GetXaxis()->GetBinWidth(ii);
      data->add(*obs,  ww*histoSpec->GetBinContent(ii));
    }
    RooAbsArg::setDirtyInhibit(false); // restore proper propagation of dirty flags
    delete histoSpec; histoSpec = 0;
  }

  RooArgSet vars(*observables), varsPlusWeight(*observables); varsPlusWeight.add(*weightVar);
  TString retName =  TString::Format("%sData", simPdf->GetName());
  asimovData = new RooDataSet(retName, "", varsPlusWeight, (weightVar ? weightVar->GetName() : 0));
  RooAbsArg::setDirtyInhibit(true); // don't propagate dirty flags while filling histograms 
  for (std::map<std::string,RooAbsData*>::iterator it = datasetPieces.begin(), ed = datasetPieces.end(); it != ed; ++it) {
    cat->setLabel(it->first.c_str());
    for (unsigned int i = 0, n = it->second->numEntries(); i < n; ++i) {
      vars = *it->second->get(i);
      asimovData->add(vars, it->second->weight());
    }
  }
  RooAbsArg::setDirtyInhibit(false); // restore proper propagation of dirty flags

  asimovData->Print("v");

  outDir->WriteTObject(asimovData, "toy_asimov");

  std::cout << "Done!" << std::endl;

}
