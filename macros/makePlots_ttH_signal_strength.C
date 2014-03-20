#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
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


void makePlots_ttH_signal_strength( TString observedFileName = "", TString expectedFileName = "", TString prefix = "", double rMin=0, double rMax=6) {

  TFile* f_obs = new TFile( observedFileName );
  TFile* f_exp = new TFile( expectedFileName );

  TGraph* g_obs = (TGraph*)f_obs->Get("gr_r_tt_deltaNLL")->Clone("g_obs");
  TGraph* g_exp = (TGraph*)f_exp->Get("gr_r_tt_deltaNLL")->Clone("g_exp");

  g_obs->SetLineWidth(2.5); g_obs->SetLineStyle(1); 
  g_exp->SetLineWidth(2.5); g_exp->SetLineStyle(7); 

  g_obs->SetMarkerStyle(20);
  g_exp->SetMarkerStyle(20);


  double xMin = rMin, xMax = rMax-0.2;
  TH1D* hdummy = new TH1D("hdummy",";#mu_{t#bar{t}H};-2 #Delta ln L",100,xMin,xMax);



  // find minimum of obs
  int n_obs = g_obs->GetN();
  double minNLL_obs = 999;
  double best_obs = 999;
  for( int i=0; i<n_obs; i++ ){
    if( g_obs->GetY()[i]<minNLL_obs ){
      minNLL_obs = g_obs->GetY()[i];
      best_obs = g_obs->GetX()[i];
    }
  }

  // find minimum of exp
  int n_exp = g_exp->GetN();
  double minNLL_exp = 999;
  double best_exp = 999;
  for( int i=0; i<n_exp; i++ ){
    if( g_exp->GetY()[i]<minNLL_exp ){
      minNLL_exp = g_exp->GetY()[i];
      best_exp = g_exp->GetX()[i];
    }
  }



  // find 1 sig and 2 sig crossing
  double yThres_1sig = 1.;
  double yThres_2sig = 3.84;

  double obs_p1sig = -99;
  double obs_m1sig = -99;
  double obs_p2sig = -99;
  double obs_m2sig = -99;
  double last_y_obs = 0;
  double last_y_exp = 0;

  double exp_p1sig = -99;
  double exp_m1sig = -99;

    
  // fix obs
  for( int i=0; i<n_obs; i++ ){
    double newX = g_obs->GetX()[i];
    double newY = 2*(g_obs->GetY()[i] - minNLL_obs);
    g_obs->SetPoint(i,newX,newY);

    if( last_y_obs>yThres_1sig && newY<yThres_1sig ) obs_m1sig = newX;
    if( last_y_obs<yThres_1sig && newY>yThres_1sig ) obs_p1sig = newX;

    if( last_y_obs>yThres_2sig && newY<yThres_2sig ) obs_m2sig = newX;
    if( last_y_obs<yThres_2sig && newY>yThres_2sig ) obs_p2sig = newX;
    last_y_obs = newY;
  }

  // fix exp
  for( int i=0; i<n_exp; i++ ){
    double newX = g_exp->GetX()[i];
    double newY = 2*(g_exp->GetY()[i] - minNLL_exp);
    g_exp->SetPoint(i,newX,newY);

    if( last_y_exp>yThres_1sig && newY<yThres_1sig ) exp_m1sig = newX;
    if( last_y_exp<yThres_1sig && newY>yThres_1sig ) exp_p1sig = newX;
    last_y_exp = newY;
  }


  printf("\t Observed: best-fit mu = %.2f  +%.2f -%.2f \n", best_obs, obs_p1sig-best_obs, best_obs-obs_m1sig );
  printf("\t Expected: best-fit mu = %.2f  +%.2f -%.2f \n", best_exp, exp_p1sig-best_exp, best_exp-exp_m1sig );



  //std::string cmsinfo =    "CMS Preliminary, ttH, b#bar{b}, #tau#tau, #gamma#gamma, WW, ZZ  #sqrt{s} = 8 TeV, L = 19.7 fb^{-1} #sqrt{s} = 7 TeV, L = 19.7 fb^{-1}";
  TString cmsinfo =   "CMS                            #sqrt{s} = 8 TeV, L = 19.7 fb^{-1};  #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}";
  //TString cmsinfo = "CMS Preliminary,    #sqrt{s} = 8 TeV, L = 19.7 fb^{-1};  #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}";
  TLatex CMSInfoLatex(0.14, 0.92, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.03);

  TLegend *legend_fake = new TLegend(0.17,0.64,0.55,0.897);
  legend_fake->SetFillColor(kWhite);
  legend_fake->SetLineColor(kWhite);
  legend_fake->SetShadowColor(kWhite);

  //TString decayinfo = "t#bar{t}H, H #rightarrow b#bar{b},#tau#tau,#gamma#gamma,WW,ZZ";
  TString decayinfo = "t#bar{t}H, H #rightarrow b#bar{b}, #tau#tau, #gamma#gamma, WW, ZZ";
  TLatex DECAYInfoLatex(0.18, 0.84, decayinfo);
  DECAYInfoLatex.SetNDC(); DECAYInfoLatex.SetTextFont(42);
  DECAYInfoLatex.SetTextSize(0.04);
  TString massinfo = "m_{H} = 125.6 GeV/c^{2}";
  TLatex MASSInfoLatex(0.18, 0.76, massinfo);
  MASSInfoLatex.SetNDC(); MASSInfoLatex.SetTextFont(42);
  MASSInfoLatex.SetTextSize(0.04);
  TString resultinfo = Form("#mu_{t#bar{t}H} = %.2f ^{+%.2f}_{ -%.2f}",best_obs,obs_p1sig-best_obs,best_obs-obs_m1sig);
  TLatex RESULTInfoLatex(0.18, 0.68, resultinfo);
  RESULTInfoLatex.SetNDC(); RESULTInfoLatex.SetTextFont(42);
  RESULTInfoLatex.SetTextSize(0.04);



  //TLegend *legend = new TLegend(0.15,0.8,0.89,0.89);
  TLegend *legend = new TLegend(0.65,0.7,0.89,0.89);

  legend->SetFillColor(kWhite);
  //legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);

  legend->AddEntry(g_obs,"Observed","l");
  legend->AddEntry(g_exp,"Expected","l");


  TCanvas* c1 = new TCanvas("c1","c1",600,600);

  c1->SetLeftMargin(.12);
  c1->SetBottomMargin(.12);
  c1->SetRightMargin(.05);

  hdummy->SetStats(0);
  hdummy->GetYaxis()->SetTitleSize(0.05);
  hdummy->GetXaxis()->SetTitleSize(0.05);


  TLine* line_obs_1sig = new TLine(xMin,yThres_1sig,xMax,yThres_1sig);
  line_obs_1sig->SetLineStyle(7);
  line_obs_1sig->SetLineColor(kRed);

  TLine* line_obs_m1sig = new TLine(obs_m1sig,0,obs_m1sig,yThres_1sig);
  line_obs_m1sig->SetLineStyle(1);
  line_obs_m1sig->SetLineColor(kRed);

  TLine* line_obs_p1sig = new TLine(obs_p1sig,0,obs_p1sig,yThres_1sig);
  line_obs_p1sig->SetLineStyle(1);
  line_obs_p1sig->SetLineColor(kRed);



  TLine* line_obs_2sig = new TLine(xMin,yThres_2sig,xMax,yThres_2sig);
  line_obs_2sig->SetLineStyle(7);
  line_obs_2sig->SetLineColor(kRed);

  TLine* line_obs_m2sig = new TLine(obs_m2sig,0,obs_m2sig,yThres_2sig);
  line_obs_m2sig->SetLineStyle(1);
  line_obs_m2sig->SetLineColor(kRed);

  TLine* line_obs_p2sig = new TLine(obs_p2sig,0,obs_p2sig,yThres_2sig);
  line_obs_p2sig->SetLineStyle(1);
  line_obs_p2sig->SetLineColor(kRed);


  hdummy->GetYaxis()->SetRangeUser(0.,8.5);

  hdummy->Draw("axis");
  g_obs->Draw("Lsame");
  g_exp->Draw("Lsame");

  line_obs_1sig->Draw();
  line_obs_m1sig->Draw();
  line_obs_p1sig->Draw();

  line_obs_2sig->Draw();
  line_obs_m2sig->Draw();
  line_obs_p2sig->Draw();

  legend->Draw();
  legend_fake->Draw();
  CMSInfoLatex.Draw();
  DECAYInfoLatex.Draw();
  MASSInfoLatex.Draw();
  RESULTInfoLatex.Draw();

  //c1->RedrawAxis();

  c1->Print("ttH_signal_strength_scan_"+prefix+".png");
  c1->Print("ttH_signal_strength_scan_"+prefix+".pdf");



  f_obs->Close();
  f_exp->Close();
  std::cout << " Done!" << std::endl;

}
