#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1.h"
#include "TF1.h"
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
#include "TKey.h"
#include "TGraph.h"
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
#include "TSystem.h"


const int LIM7TEV = 0x1;
const int LIM8TEV = 0x2;

void makePlots_limit_v4(TString limitFileName, std::string  pname, std::string  ptitle, int ana, float maxHeight = -1, bool addObs = true, bool add_inj=true) {
  
  TCanvas *c1 = new TCanvas("c1");
  

  const int nMax = 100;
  double x[nMax], y[nMax], y2[nMax], y_obs[nMax], y_inj[nMax];
  double ymin_1sig[nMax], ymax_1sig[nMax], ymin_2sig[nMax], ymax_2sig[nMax];
    
 
  ifstream infile;
  infile.open(limitFileName.Data());
  if(!infile) { // file couldn't be opened
    cerr << "Error: file " << limitFileName.Data() << " could not be opened" << endl;
    return;
  }

  //Table header
  std::cout << "{\\small" << std::endl
            << "\\begin{tabular}{|c|c|ccc|} \\hline " << std::endl
            << "           &          & \\multicolumn{3}{c|}{Expected} \\\\" << std::endl
            << "Higgs Mass & Observed & Median & 68\\% C.L. Range &  95\\% C.L. Range  \\\\ \\hline" << std::endl;

  TString inLine = "";

  int i = 0;
  double xMin = 9e20, xMax = -9e20;
  while (inLine.ReadLine(infile)) {

    if (inLine.BeginsWith("#")) continue;

    if(nMax < i+1){ // Too many mass bins?
      cerr << "Error: nMax (" << nMax << ") < i+1 (" << i+1 << ")" << endl;
      return;
    }

    //Read mass
    x[i] = atof(inLine.Data());

    if (x[i] > xMax) xMax = x[i];
    if (x[i] < xMin) xMin = x[i];

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);
    
    if (readSuccess) {
      TString obsStr = inLine.Data();
      if (obsStr == "NO OBS") {
        y_obs[i] = -9e20;
        addObs = false; //If we have points without observed limits, don't bother plotting
      } else {
        y_obs[i] = atof(inLine.Data());
      }

    } else {
      std::cerr << "Missing observed limit for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);

    if (readSuccess) {
      ymin_2sig[i] = atof(inLine.Data());
    } else {
      std::cerr << "Missing -2 sigma for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);

    if (readSuccess) {
      ymin_1sig[i] = atof(inLine.Data());
    } else {
      std::cerr << "Missing -1 sigma for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);

    if (readSuccess) {
      y[i] = atof(inLine.Data());
    } else {
      std::cerr << "Missing +2 sigma for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);

    if (readSuccess) {
      ymax_1sig[i] = atof(inLine.Data());
    } else {
      std::cerr << "Missing +1 sigma for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }

    bool readSuccess = inLine.ReadLine(infile);
    while (readSuccess && inLine.BeginsWith("#")) readSuccess = inLine.ReadLine(infile);
    
    if (readSuccess) {
      ymax_2sig[i] = atof(inLine.Data());
    } else {
      std::cerr << "Missing +2 sigma for mass point " << i << ", m = " << x[i] << std::endl;
      return;
    }
    
    cout << "$" << x[i] << "~\\GeVcc$ & ";
    cout << (addObs ?  Form("%.1f",y_obs[i]) : " X ");
    cout << " & ";
    cout << Form("%.1f",y[i]) << " & ["
         << Form("%.1f",ymin_1sig[i]) << ","
         << Form("%.1f",ymax_1sig[i]) << "] & ["
         << Form("%.1f",ymin_2sig[i]) << ","
         << Form("%.1f",ymax_2sig[i]) << "] \\\\" 
         << endl;

    //Increment point
    ++i;

  }

  ifstream injfile;
  if (add_inj) {
    injfile.open(limitFileName.ReplaceAll(".", "_inj.").Data());
    if (!injfile) {
      cerr << "Error: file " << limitFileName.ReplaceAll(".", "_inj.").Data() << " could not be opened" << endl;
      return;
    }

    int j = 0;
    TString injline = "";
    while (injline.ReadLine(injfile)) {
      if (injline.BeginsWith("#"))
        continue;

      float x_tmp = atof(injline.Data());
      if (x[j] != x_tmp) {
        cerr << "Error: injected and expected mass points are different!" << endl;
        cerr << "x[" << j << "] = " << x[j] << " != " << x_tmp << endl;
        return;
      }

      bool succ = injline.ReadLine(injfile);
      while (succ && injline.BeginsWith("#")) succ = injline.ReadLine(injfile);
      y_inj[j++] = atof(injline.Data());
    }
    injfile.close();
  }

  std::cout << "\\hline" << std::endl
            << "\\end{tabular}}" << std::endl;

  int n = i;

  std::cout << "i = " << i << std::endl;
	
  if (maxHeight < 0) maxHeight = ymax_2sig[n-1] ;



  TLatex *CMSInfoLatex = new TLatex(0.11, 0.91, ("CMS preliminary          " + ptitle).c_str());
  CMSInfoLatex->SetNDC();
  CMSInfoLatex->SetTextFont(42);

  TString lumiinfo = " ";
  if (ana & LIM7TEV) {
    lumiinfo += "#sqrt{s} = 7 TeV, L = 5.0 fb^{-1}";
  }

  double textSize = 0.04;
  double offset = 0.0;
  textSize = 0.03;
  if (ana == (LIM7TEV|LIM8TEV)) {
    lumiinfo += "; ";
    offset = 0.07;
  }
  offset = -0.07;

  if (ana & LIM8TEV) {
    lumiinfo += "#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}";
  }


  TLatex *LUMIInfoLatex = new TLatex(0.54-offset, 0.91, lumiinfo);
  LUMIInfoLatex->SetNDC();
  LUMIInfoLatex->SetTextFont(42);

  //Set to same size
  CMSInfoLatex->SetTextSize(textSize);
  LUMIInfoLatex->SetTextSize(textSize);


  TH1D* h_dummy = new TH1D("h_dummy","",100, xMin, xMax );

  TGraph *gr     = new TGraph(n,x,y);
  TGraph *gr_obs = new TGraph(n,x,y_obs);
  TGraph *grshade_1sig = new TGraph(2*n);
  TGraph *grshade_2sig = new TGraph(2*n);

  for( int i=0;i<n;i++) {
    grshade_1sig->SetPoint(i,x[i],ymax_1sig[i]);
    grshade_1sig->SetPoint(n+i,x[n-i-1],ymin_1sig[n-i-1]);

    grshade_2sig->SetPoint(i,x[i],ymax_2sig[i]);
    grshade_2sig->SetPoint(n+i,x[n-i-1],ymin_2sig[n-i-1]);
  }

  grshade_1sig->SetFillColor(kGreen);
  grshade_1sig->SetLineColor(1);
  grshade_1sig->SetLineStyle(2);
  grshade_2sig->SetFillColor(kYellow);
  grshade_2sig->SetLineColor(1);
  grshade_2sig->SetLineStyle(2);

  TGraph *gr_inj = new TGraph(n, x, y_inj);
  gr_inj->SetLineWidth(2);
  gr_inj->SetLineStyle(2);
  gr_inj->SetLineColor(kRed);

  gr->SetLineWidth(2);
  gr->SetLineStyle(2);

  gr_obs->SetLineWidth(2);
  gr_obs->SetLineStyle(1);
  gr_obs->SetLineColor(1);
  gr_obs->SetMarkerStyle(20);
  gr_obs->SetMarkerColor(1);
  gr_obs->SetMarkerSize(0.9);

  h_dummy->SetTitle(";m_{H} (GeV);95% CL limit on #sigma/#sigma_{SM}");


  TLegend *legend     = new TLegend(0.17,0.65,0.42,0.87);
  TLegend *legend_obs = new TLegend(0.17,0.59,0.45,0.87);

  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);

  legend->AddEntry(grshade_1sig,"Expected #pm 1#sigma","fl");
  legend->AddEntry(grshade_2sig,"Expected #pm 2#sigma","fl");

  legend_obs->SetFillColor(kWhite);
  legend_obs->SetLineColor(kWhite);
  legend_obs->SetShadowColor(kWhite);
  legend_obs->SetTextFont(42);
  legend_obs->SetTextSize(0.04);

  if (addObs) legend_obs->AddEntry(gr_obs, "Observed ","lp");
  legend_obs->AddEntry(gr_inj, "ttH(125) injected", "l");
  legend_obs->AddEntry(grshade_1sig,"Expected #pm 1#sigma","fl");
  legend_obs->AddEntry(grshade_2sig,"Expected #pm 2#sigma","fl");


  h_dummy->SetStats(0);    
  h_dummy->GetYaxis()->SetRangeUser(0.,20.);



  h_dummy->SetMaximum(1.05*maxHeight);
  cout << "maxHeight is " << maxHeight << endl;
  h_dummy->Draw();

  c1->SetGridx(1);
  c1->SetGridy(1);

  grshade_2sig->Draw("f");
  grshade_1sig->Draw("f");
  gr->Draw("l");
  if (add_inj)
    gr_inj->Draw("l");
  if (addObs) gr_obs->Draw("pl");


  TLine* line = new TLine( xMin, 1., xMax, 1. );
  line->SetLineColor(kRed);
  line->SetLineWidth(4);
  line->Draw();

  gPad->RedrawAxis();
  gPad->RedrawAxis("G");

  if (addObs) legend_obs->Draw();
  else        legend->Draw();

  CMSInfoLatex->Draw();
  LUMIInfoLatex->Draw();

  TString dirprefix = "plots";

  int dirOK = 0;
  if (gSystem->AccessPathName(dirprefix)) { //Backwards: False means path is there, true means not there!
    dirOK = gSystem->MakeDirectory(dirprefix);
  }

  if (dirOK == 0) {

    TString plotname = dirprefix;
    plotname += "/";
    plotname += pname;
    plotname += "_ttH_mH_limit_Exp";
    if (addObs || add_inj) {
      if (addObs)
        plotname += "AndObs";
      if (add_inj)
        plotname += "AndInj";
    } else {
      plotname += "Only";
    }
    plotname += "_shape";
    c1->Print(plotname+".pdf");
    c1->Print(plotname+".png");
  } else {
    std::cout << "Couldn't write plots to directory " << dirprefix << std::endl;
  }

}


