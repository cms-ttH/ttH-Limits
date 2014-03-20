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


void makePlots_ttH_coupling( TString inputFileName = "", TString prefix = "") {

  TFile* file = new TFile( inputFileName );

  TH2D* h_CVCF_deltaNLL = (TH2D*)file->Get("h_CVCF_deltaNLL");
  TGraph* gr68_p_old = (TGraph*)file->Get("gr68_p")->Clone("gr68_p_old");
  TGraph* gr95_p_old = (TGraph*)file->Get("gr95_p")->Clone("gr95_p_old");
  TGraph* gr95_n_old = (TGraph*)file->Get("gr95_n")->Clone("gr95_n_old");
  TGraph* gr_SM1 = (TGraph*)file->Get("gr_SM")->Clone("gr_SM1");
  TGraph* gr_SM2 = (TGraph*)file->Get("gr_SM")->Clone("gr_SM2");
  TGraph* gr_best = (TGraph*)file->Get("gr_best");



  printf("\t best-fit k_{V} = %.2f, k_{f} = %.2f \n", gr_best->GetX()[0], gr_best->GetY()[0] );



  TString cmsinfo =   "CMS";
  //TString cmsinfo =   "CMS                                        #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
  //TString cmsinfo = "CMS                             #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
  //TString cmsinfo = "CMS Preliminary                 #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
  //TString cmsinfo =    "CMS Preliminary,    #sqrt{s} = 8 TeV, L = 19.7 fb^{-1};  #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}";
  TLatex CMSInfoLatex(0.12, 0.94, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.04);

  TString lumiinfo = "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";

  TLatex LumiInfoLatex(0.58, 0.94, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.035);


  TLegend *legend_fake = new TLegend(0.15,0.62,0.55,0.897);
  legend_fake->SetFillColor(kWhite);
  legend_fake->SetLineColor(kWhite);
  legend_fake->SetShadowColor(kWhite);

  TString decayinfo = "t#bar{t}H, H #rightarrow b#bar{b}, #tau#tau, #gamma#gamma, WW, ZZ";
  TLatex DECAYInfoLatex(0.44, 0.87, decayinfo);
  DECAYInfoLatex.SetNDC(); DECAYInfoLatex.SetTextFont(42);
  DECAYInfoLatex.SetTextSize(0.04);
  TString massinfo = "m_{H} = 125.6 GeV/c^{2}";
  TLatex MASSInfoLatex(0.57, 0.81, massinfo);
  MASSInfoLatex.SetNDC(); MASSInfoLatex.SetTextFont(42);
  MASSInfoLatex.SetTextSize(0.04);


  TLatex DECAYInfoLatex1(0.48, 0.87, decayinfo);
  DECAYInfoLatex1.SetNDC(); DECAYInfoLatex1.SetTextFont(42);
  DECAYInfoLatex1.SetTextSize(0.04);
  TLatex MASSInfoLatex1(0.62, 0.81, massinfo);
  MASSInfoLatex1.SetNDC(); MASSInfoLatex1.SetTextFont(42);
  MASSInfoLatex1.SetTextSize(0.04);
//   TString resultinfo = Form("#mu_{t#bar{t}H} = %.2f ^{+%.2f}_{ -%.2f}",best_obs,obs_p1sig-best_obs,best_obs-obs_m1sig);
//   TLatex RESULTInfoLatex(0.16, 0.66, resultinfo);
//   RESULTInfoLatex.SetNDC(); RESULTInfoLatex.SetTextFont(42);
//   RESULTInfoLatex.SetTextSize(0.04);

  TString deltaNLLinfo = "- 2 #Delta ln L";
  TLatex DELTANLLInfoLatex(0.98, 0.72, deltaNLLinfo);
  DELTANLLInfoLatex.SetNDC(); DELTANLLInfoLatex.SetTextFont(42);
  DELTANLLInfoLatex.SetTextSize(0.04);
  DELTANLLInfoLatex.SetTextAngle(90);



//   TLegend *legend_lines = new TLegend(0.65,0.15,0.85,0.25);
//   //legend_lines->SetFillColor(kWhite);
//   //legend_lines->SetLineColor(kWhite);
//   //legend_lines->SetShadowColor(kWhite);
//   legend_lines->SetFillStyle(4000);
//   legend_lines->SetBorderSize(0);
//   legend_lines->SetTextFont(42);
//   legend_lines->SetTextSize(0.04);


  TLegend *legend_lines_nocontour = new TLegend(0.68,0.13,0.90,0.31);
  //legend_lines->SetFillColor(kWhite);
  legend_lines_nocontour->SetLineColor(kBlack);
  //legend_lines->SetShadowColor(kWhite);
  legend_lines_nocontour->SetFillStyle(4000);
  //legend_lines_nocontour->SetBorderSize(0);
  legend_lines_nocontour->SetTextFont(42);
  legend_lines_nocontour->SetTextSize(0.04);

  int npoints68_p = gr68_p_old->GetN();
  int npoints95_p = gr95_p_old->GetN();
  int npoints95_n = gr95_n_old->GetN();
  Double_t *xi68_p = gr68_p_old->GetX(), *yi68_p = gr68_p_old->GetY();
  Double_t *xi95_p = gr95_p_old->GetX(), *yi95_p = gr95_p_old->GetY();
  Double_t *xi95_n = gr95_n_old->GetX(), *yi95_n = gr95_n_old->GetY();

  int reduction68_p = 15;//20 //25
  int reduction95_p = 10;
  int reduction95_n = 10;
  int npoints_new68_p = ( npoints68_p%reduction68_p==0 ) ? npoints68_p/reduction68_p : (npoints68_p+npoints68_p%reduction68_p)/reduction68_p;
  int npoints_new95_p = ( npoints95_p%reduction95_p==0 ) ? npoints95_p/reduction95_p : (npoints95_p+npoints95_p%reduction95_p)/reduction95_p;
  int npoints_new95_n = ( npoints95_n%reduction95_n==0 ) ? npoints95_n/reduction95_n : (npoints95_n+npoints95_n%reduction95_n)/reduction95_n;

  double x68_p_new[npoints_new68_p+1], y68_p_new[npoints_new68_p+1];
  double x95_p_new[npoints_new95_p+1], y95_p_new[npoints_new95_p+1];
  double x95_n_new[npoints_new95_n+1], y95_n_new[npoints_new95_n+1];
  int iPoint68_p_new=0;
  int iPoint95_p_new=0;
  int iPoint95_n_new=0;
  for( int iPoint=0; iPoint<npoints68_p; iPoint++ ){
    if( iPoint%reduction68_p!=0 ) continue;
    std::cout << iPoint << "\t" << xi68_p[iPoint] << "\t" << yi68_p[iPoint] << std::endl;
    x68_p_new[iPoint68_p_new] = xi68_p[iPoint];
    y68_p_new[iPoint68_p_new] = yi68_p[iPoint];
    iPoint68_p_new++;
  }

  x68_p_new[iPoint68_p_new] = xi68_p[0];
  y68_p_new[iPoint68_p_new] = yi68_p[0];
  iPoint68_p_new++;


  for( int iPoint=0; iPoint<npoints95_p; iPoint++ ){
    if( iPoint%reduction95_p!=0 ) continue;
    x95_p_new[iPoint95_p_new] = xi95_p[iPoint];
    y95_p_new[iPoint95_p_new] = yi95_p[iPoint];
    iPoint95_p_new++;
  }

  x95_p_new[iPoint95_p_new] = xi95_p[0];
  y95_p_new[iPoint95_p_new] = yi95_p[0];
  iPoint95_p_new++;


  for( int iPoint=0; iPoint<npoints95_n; iPoint++ ){
    if( iPoint%reduction95_n!=0 ) continue;
    //std::cout << iPoint << "\t" << xi95_n[iPoint] << "\t" << yi95_n[iPoint] << std::endl;
    x95_n_new[iPoint95_n_new] = xi95_n[iPoint];
    y95_n_new[iPoint95_n_new] = yi95_n[iPoint];
    iPoint95_n_new++;
  }

  x95_n_new[iPoint95_n_new] = xi95_n[0];
  y95_n_new[iPoint95_n_new] = yi95_n[0];
  iPoint95_n_new++;



  TGraph *gr68_p = new TGraph(iPoint68_p_new,x68_p_new,y68_p_new);
  TGraph *gr95_p = new TGraph(iPoint95_p_new,x95_p_new,y95_p_new);
  TGraph *gr95_n = new TGraph(iPoint95_n_new,x95_n_new,y95_n_new);


  gr68_p->SetLineStyle(1);
  gr95_p->SetLineStyle(1);
  gr95_n->SetLineStyle(1);
  gr68_p->SetLineColor(82);
  gr95_p->SetLineColor(89);
  gr95_n->SetLineColor(89);

  h_CVCF_deltaNLL->GetXaxis()->SetTitleSize(0.06);
  h_CVCF_deltaNLL->GetYaxis()->SetTitleSize(0.06);
  h_CVCF_deltaNLL->GetXaxis()->SetTitleOffset(0.7);
  h_CVCF_deltaNLL->GetYaxis()->SetTitleOffset(0.7);


  for( int yBin=0; yBin<h_CVCF_deltaNLL->GetNbinsY(); yBin++ ) 
    h_CVCF_deltaNLL->SetBinContent(h_CVCF_deltaNLL->GetNbinsX()-1,yBin+1,h_CVCF_deltaNLL->GetBinContent(h_CVCF_deltaNLL->GetNbinsX()-2,yBin+1));
    for( int yBin=0; yBin<h_CVCF_deltaNLL->GetNbinsY(); yBin++ ) 
    h_CVCF_deltaNLL->SetBinContent(h_CVCF_deltaNLL->GetNbinsX(),yBin+1,h_CVCF_deltaNLL->GetBinContent(h_CVCF_deltaNLL->GetNbinsX()-1,yBin+1));
    
  //h_CVCF_deltaNLL->GetXaxis()->SetRangeUser(0.,3.0);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  h_CVCF_deltaNLL->SetStats(0);

  c1->SetTopMargin(0.08);
  c1->SetRightMargin(0.08);


  h_CVCF_deltaNLL->Draw("axis");

  gr68_p->SetLineColor(kBlack);
  gr95_p->SetLineColor(kBlack);
  gr95_n->SetLineColor(kBlack);

  gr68_p->SetLineStyle(1);
  gr95_p->SetLineStyle(7);
  gr95_n->SetLineStyle(7);

  gr68_p->SetLineWidth(2);
  gr95_p->SetLineWidth(2);
  gr95_n->SetLineWidth(2);


  gr68_p->Draw("C SAME");
  gr95_p->Draw("L SAME");
  gr95_n->Draw("L SAME");

  gr_SM1->SetMarkerStyle(33);
  gr_SM1->SetMarkerColor(kRed);
  gr_SM1->SetMarkerSize(3);

  TGraph* gr_SM1_dummy = (TGraph*)gr_SM1->Clone("gr_SM1_dummy");

  gr_SM2->SetMarkerStyle(33);
  gr_SM2->SetMarkerColor(kOrange);
  gr_SM2->SetLineColor(kRed);
  gr_SM2->SetLineWidth(2);

  TGraph* gr_SM2_dummy = (TGraph*)gr_SM2->Clone("gr_SM1_dummy");
  gr_SM2_dummy->SetPoint(0,4.444,-2.61);

  gr_SM1->Draw("P SAME");
  gr_SM2->Draw("P SAME");
  gr_best->Draw("P SAME");

  legend_lines_nocontour->AddEntry(gr_best,"Best fit","p");
  legend_lines_nocontour->AddEntry(gr68_p,"68% CL","l");
  legend_lines_nocontour->AddEntry(gr95_p,"95% CL","l");
  legend_lines_nocontour->AddEntry(gr_SM1_dummy,"SM Higgs","p");

  legend_lines_nocontour->Draw();
  gr_SM2_dummy->Draw("P SAME");

  CMSInfoLatex.Draw();
  DECAYInfoLatex1.Draw();
  MASSInfoLatex1.Draw();

  LumiInfoLatex.Draw();

  c1->Print("ttH_cVcF_coupling_nocontour_"+prefix+".png");
  c1->Print("ttH_cVcF_coupling_nocontour_"+prefix+".pdf");



  TCanvas *c2 = new TCanvas("c2","c2",650,600);

  c2->SetTopMargin(0.08);
  c2->SetRightMargin(0.08);
  c2->SetRightMargin(0.15);
  gStyle->SetNumberContours(50);

  //h_CVCF_deltaNLL->GetYaxis()->SetRangeUser(-2.6,3.0);

  h_CVCF_deltaNLL->Draw("colz");

  gr68_p->SetLineColor(kBlack);
  gr95_p->SetLineColor(kBlack);

  gr68_p->SetLineStyle(1);
  gr95_p->SetLineStyle(7);

  gr68_p->SetLineWidth(2);
  gr95_p->SetLineWidth(2);


  TLegend *legend_lines = new TLegend(0.67,0.13,0.83,0.31);
  //legend_lines->SetFillColor(kWhite);
  legend_lines->SetLineColor(kBlack);
  //legend_lines->SetShadowColor(kWhite);
  legend_lines->SetFillStyle(4000);
  //legend_lines->SetBorderSize(0);
  legend_lines->SetTextFont(42);
  legend_lines->SetTextSize(0.03);


  legend_lines->AddEntry(gr_best,"Best fit","p");
  legend_lines->AddEntry(gr68_p,"68% CL","l");
  legend_lines->AddEntry(gr95_p,"95% CL","l");
  legend_lines->AddEntry(gr_SM1_dummy,"SM Higgs","p");

//   legend_lines->AddEntry(gr68_p," 68% CL","l");
//   legend_lines->AddEntry(gr95_p," 95% CL","l");

  gr68_p->Draw("C SAME");
  gr95_p->Draw("L SAME");
  gr95_n->Draw("L SAME");


  gr_SM1->Draw("P SAME");
  gr_SM2->Draw("P SAME");
  gr_best->Draw("P SAME");

  legend_lines->Draw();
  gr_SM2_dummy->SetPoint(0,4.72,-2.61);
  gr_SM2_dummy->Draw("P SAME");

  CMSInfoLatex.Draw();
  DECAYInfoLatex.Draw();
  MASSInfoLatex.Draw();
  DELTANLLInfoLatex.Draw();

  LumiInfoLatex.Draw();

  c2->Print("ttH_cVcF_coupling_contour_"+prefix+".png");
  c2->Print("ttH_cVcF_coupling_contour_"+prefix+".pdf");


  file->Close();
  std::cout << " Done!" << std::endl;

}
