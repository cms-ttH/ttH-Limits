#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TList.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include <sstream>


TGraph* bestFit(double xbest, double ybest);
TGraph* contourPlot(TH2D* hist, TGraph *bestFit, int side=1);


void runOnHiggsCombine_CV_CF( int maxNentries=-1, TString inputFile="higgsCombine_CVCF_LJ_OSDIL_grid_1000.MultiDimFit.mH125.root", TString postfix = "", int nBinsX_=120, int nBinsY_=120, double xMax_=3.0, double yMax_=3.0 ){

  TFile *treefile = new TFile( inputFile );

  TTree *tree = (TTree*)treefile->Get("limit");


  Float_t CV, CF, quantileExpected, deltaNLL;

  tree->SetBranchAddress("CV", &CV );
  tree->SetBranchAddress("CF", &CF );
  tree->SetBranchAddress("quantileExpected", &quantileExpected );
  tree->SetBranchAddress("deltaNLL", &deltaNLL );


  TString histofilename = "HistoFiles/runOnHiggsCombine"+postfix+".root";

  TFile histofile(histofilename,"recreate");
  histofile.cd();


  TH1D* h_quantileExpected = new TH1D("h_quantileExpected",";quantileExpected",10000,-1,1);
  TH1D* h_deltaNLL = new TH1D("h_deltaNLL",";2*deltaNLL",10000,-3,3);

  int nbinsX = nBinsX_;
  int nbinsY = nBinsY_;
  double xMax = xMax_;
  double yMax = yMax_;

  TH2D* h_CVCF_quantileExpected = new TH2D("h_CVCF_quantileExpected",";#kappa_{V};#kappa_{f}",nbinsX,0,xMax,nbinsY,-yMax,yMax);
  TH2D* h_CVCF_deltaNLL_full = new TH2D("h_CVCF_deltaNLL_full",";#kappa_{V};#kappa_{f}",nbinsX,0,xMax,nbinsY,-yMax,yMax);
  TH2D* h_CVCF_deltaNLL = new TH2D("h_CVCF_deltaNLL",";#kappa_{V};#kappa_{f}",nbinsX,0,xMax,nbinsY,-yMax,yMax);

  TH2D* h_CVCF_CL68 = new TH2D("h_CVCF_CL68",";#kappa_{V};#kappa_{f}",nbinsX,0,xMax,nbinsY,-yMax,yMax);
  TH2D* h_CVCF_CL95 = new TH2D("h_CVCF_CL95",";#kappa_{V};#kappa_{f}",nbinsX,0,xMax,nbinsY,-yMax,yMax);

  int nentries = tree->GetEntries();

  double best_CV = -99, best_CF = -99;

  for (Long64_t ievt=0; ievt<tree->GetEntries();ievt++) {    //Long64_t

    if (ievt%1000 == 0){
      std::cout << "--- ... Processing event: " << ievt << std::endl;
    }

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%100000==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					    << int(double(ievt)/double(nentries)*100) << "% done" << std::endl;

    if( ievt==maxNentries ) break;

    tree->GetEntry(ievt);

    h_quantileExpected->Fill( quantileExpected );

    if( quantileExpected<0 ) continue;

    h_CVCF_deltaNLL_full->SetBinContent( h_CVCF_deltaNLL->FindBin(CV,CF) , 2*deltaNLL );
    h_CVCF_quantileExpected->SetBinContent( h_CVCF_quantileExpected->FindBin(CV,CF) , quantileExpected );

    //if( deltaNLL<=1 ){
      h_CVCF_deltaNLL->SetBinContent( h_CVCF_deltaNLL->FindBin(CV,CF) , 2*deltaNLL );
      //}
    h_deltaNLL->Fill( 2*deltaNLL );
 

    if( quantileExpected==1 ){
      best_CV = CV;
      best_CF = CF;
    }

  }


  int NbinsX = h_CVCF_quantileExpected->GetNbinsX();
  int NbinsY = h_CVCF_quantileExpected->GetNbinsY();

  for( int ybin=0; ybin<NbinsY; ybin++ ){
    for( int xbin=0; xbin<NbinsX; xbin++ ){
      if( ybin==0 || ybin==NbinsY-1 || xbin==0 || xbin==NbinsX-1 ) h_CVCF_quantileExpected->SetBinContent(xbin+1,ybin+1,0);
    }
  }

  double use68 = 0.32;
  double use95 = 0.05;


  for( int xbin=0; xbin<NbinsX; xbin++ ){
    double previous_bin_content = 0;
    for( int ybin=1; ybin<NbinsY; ybin++ ){
      double content = h_CVCF_quantileExpected->GetBinContent(xbin+1,ybin+1);
      if( (previous_bin_content<use68 && content>=use68) ||
	  (previous_bin_content>use68 && content<=use68) ) h_CVCF_CL68->SetBinContent(xbin+1,ybin+1,1);
      if( (previous_bin_content<use95 && content>=use95) ||
	  (previous_bin_content>use95 && content<=use95) ) h_CVCF_CL95->SetBinContent(xbin+1,ybin+1,1);
      previous_bin_content = content;
    }
  }

  for( int ybin=0; ybin<NbinsY; ybin++ ){
    double previous_bin_content = 0;
    for( int xbin=1; xbin<NbinsX; xbin++ ){
      double content = h_CVCF_quantileExpected->GetBinContent(xbin+1,ybin+1);
      if( (previous_bin_content<use68 && content>=use68) ||
	  (previous_bin_content>use68 && content<=use68) ) h_CVCF_CL68->SetBinContent(xbin+1,ybin+1,1);
      if( (previous_bin_content<use95 && content>=use95) ||
	  (previous_bin_content>use95 && content<=use95) ) h_CVCF_CL95->SetBinContent(xbin+1,ybin+1,1);
      previous_bin_content = content;
    }
  }


  std::cout << "Best-fit value at #kappa_{f} = " << best_CF << ", #kappa_{V} = " << best_CV << std::endl;

  int bin_best_CV = h_CVCF_deltaNLL_full->GetXaxis()->FindBin(best_CV);
  int bin_sm_CV = h_CVCF_deltaNLL_full->GetXaxis()->FindBin(1.0);

  int bin_best_CF = h_CVCF_deltaNLL_full->GetYaxis()->FindBin(best_CF);
  int bin_sm_CF = h_CVCF_deltaNLL_full->GetYaxis()->FindBin(1.0);

  TH1D* h_deltaNLL_CF_bestCV = (TH1D*)h_CVCF_deltaNLL_full->ProjectionY("h_deltaNLL_CF_bestCV",bin_best_CV,bin_best_CV);
  TH1D* h_deltaNLL_CF_smCV = (TH1D*)h_CVCF_deltaNLL_full->ProjectionY("h_deltaNLL_CF_smCV",bin_sm_CV,bin_sm_CV);

  TH1D* h_deltaNLL_CV_bestCF = (TH1D*)h_CVCF_deltaNLL_full->ProjectionX("h_deltaNLL_CV_bestCF",bin_best_CF,bin_best_CF);
  TH1D* h_deltaNLL_CV_smCF = (TH1D*)h_CVCF_deltaNLL_full->ProjectionX("h_deltaNLL_CV_smCF",bin_sm_CF,bin_sm_CF);

  TH1D* h_quantileExpected_CF_bestCV = (TH1D*)h_CVCF_quantileExpected->ProjectionY("h_quantileExpected_CF_bestCV",bin_best_CV,bin_best_CV);
  TH1D* h_quantileExpected_CF_smCV = (TH1D*)h_CVCF_quantileExpected->ProjectionY("h_quantileExpected_CF_smCV",bin_sm_CV,bin_sm_CV);

  TH1D* h_quantileExpected_CV_bestCF = (TH1D*)h_CVCF_quantileExpected->ProjectionX("h_quantileExpected_CV_bestCF",bin_best_CF,bin_best_CF);
  TH1D* h_quantileExpected_CV_smCF = (TH1D*)h_CVCF_quantileExpected->ProjectionX("h_quantileExpected_CV_smCF",bin_sm_CF,bin_sm_CF);

  TGraph *gr_best = bestFit(best_CV,best_CF);
  TGraph *gr_SM = bestFit(1.,1.);
  TGraph *gr68_p = contourPlot(h_CVCF_CL68,gr_best,1);
  TGraph *gr95_p = contourPlot(h_CVCF_CL95,gr_best,1);
  TGraph *gr68_n = contourPlot(h_CVCF_CL68,gr_best,-1);
  TGraph *gr95_n = contourPlot(h_CVCF_CL95,gr_best,-1);

  if( gr68_p!=NULL ){
    gr68_p->SetLineWidth(1); gr68_p->SetLineStyle(1); gr68_p->SetLineColor(1); gr68_p->SetFillStyle(1001); gr68_p->SetFillColor(82); 
  }

  if( gr95_p!=NULL ){
    gr95_p->SetLineWidth(1); gr95_p->SetLineStyle(7); gr95_p->SetLineColor(1); gr95_p->SetFillStyle(1001); gr95_p->SetFillColor(89);
  }

  if( gr68_n!=NULL ){
    gr68_n->SetLineWidth(1); gr68_n->SetLineStyle(1); gr68_n->SetLineColor(1); gr68_n->SetFillStyle(1001); gr68_n->SetFillColor(82);  
  }

  if( gr95_n!=NULL ){
    gr95_n->SetLineWidth(1); gr95_n->SetLineStyle(7); gr95_n->SetLineColor(1); gr95_n->SetFillStyle(1001); gr95_n->SetFillColor(89);
  }

  gr_SM->SetMarkerStyle(33);
  gr_SM->SetMarkerColor(kBlue);

  gr_SM->Write("gr_SM");
  gr_best->Write("gr_best");
  if( gr68_p!=NULL ) gr68_p->Write("gr68_p");
  if( gr95_p!=NULL ) gr95_p->Write("gr95_p");
  if( gr68_n!=NULL ) gr68_n->Write("gr68_n");
  if( gr95_n!=NULL ) gr95_n->Write("gr95_n");

  treefile->Close();

  histofile.Write();
  histofile.Close();

  std::cout << " Done! " << std::endl;

}


TGraph* bestFit(double xbest, double ybest) {

  Double_t x[1], y[1];
  x[0] = xbest;
  y[0] = ybest;

  TGraph *gr0 = new TGraph(1,x,y);
  gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
  return gr0;
}



TGraph* contourPlot(TH2D* hist, TGraph *bestFit, int side) {

    std::vector<double> use_xi,use_yi;
    for( int xbin=0; xbin<hist->GetNbinsX(); xbin++ ){
      for( int ybin=0; ybin<hist->GetNbinsY(); ybin++ ){
	if( side>0 && hist->GetYaxis()->GetBinCenter(ybin+1)<=0 ) continue;
	if( side<0 && hist->GetYaxis()->GetBinCenter(ybin+1)>=0 ) continue;

	if( hist->GetBinContent(xbin+1,ybin+1)>0 ){
	  use_xi.push_back( hist->GetXaxis()->GetBinCenter(xbin+1) );
	  use_yi.push_back( hist->GetYaxis()->GetBinCenter(ybin+1) );
	}
      }
    }

    int size = int( use_xi.size() );
    std::vector<int> usedPoint;
    usedPoint.clear();
    
    if( size==0 ){
      std::cout << " No points in range " << std::endl;
      TGraph *grNULL = NULL;
      return grNULL;
    }

    double lastX = use_xi[0], lastY = use_yi[0];
    usedPoint.push_back(0);

    std::vector<double> xi,yi;
    xi.push_back(lastX);
    yi.push_back(lastY);
    for( int iPoint=1; iPoint<size; iPoint++ ){
      double minDist = 999;
      int minDistPoint = -1;
      for( int jPoint=0; jPoint<size; jPoint++ ){
	bool alreadyUsed = false;
	for( int kPoint=0; kPoint<int(usedPoint.size()); kPoint++ ){
	  if( jPoint==usedPoint[kPoint] ){
	    alreadyUsed = true;
	    break;
	  }
	}
	if( alreadyUsed ) continue;

	double curX = use_xi[jPoint];
	double curY = use_yi[jPoint];

	double dist2 = (curX-lastX)*(curX-lastX) + (curY-lastY)*(curY-lastY);

	if( dist2<minDist ){
	  minDist = dist2;
	  minDistPoint = jPoint;
	}
      }

      usedPoint.push_back(minDistPoint);
      lastX = use_xi[minDistPoint];
      lastY = use_yi[minDistPoint];
      xi.push_back(lastX);
      yi.push_back(lastY);
    }


    xi.push_back(use_xi[0]);
    yi.push_back(use_yi[0]);

    size = int( xi.size() );
    double myxi[size];
    double myyi[size];
    for( int i=0; i<size; i++ ){
      myxi[i] = xi[i];
      myyi[i] = yi[i];
    }

    TGraph *gr = new TGraph( size, myxi, myyi );
    return gr;
}


/*


TCanvas *c1 = new TCanvas("c1","c1");
h_deltaNLL_CF_bestCV->SetStats(0);
h_deltaNLL_CF_bestCV->SetMaximum(3);
h_deltaNLL_CF_bestCV->GetXaxis()->SetRangeUser(-3.2,3.2);
h_deltaNLL_CF_bestCV->Draw();
h_deltaNLL_CF_smCV->SetLineColor(kRed);
h_deltaNLL_CF_smCV->Draw("same");


gr68->SetLineStyle(1);
gr95->SetLineStyle(1);
gr68->SetLineColor(82);
gr95->SetLineColor(89);

TLegend *legend = new TLegend(0.6,0.67,0.89,0.89);

legend->SetFillColor(kWhite);
legend->SetLineColor(kWhite);
legend->SetShadowColor(kWhite);
legend->SetTextFont(42);
legend->SetTextSize(0.04);
//legend->SetNColumns(1);

legend->AddEntry(gr_SM,"SM Higgs","p");
legend->AddEntry(gr_best, "Best-fit","p");
legend->AddEntry(gr68,"68% CL","f");
legend->AddEntry(gr95,"95% CL","f");

h_CVCF_deltaNLL->GetXaxis()->SetTitleSize(0.06);
h_CVCF_deltaNLL->GetYaxis()->SetTitleSize(0.06);
h_CVCF_deltaNLL->GetXaxis()->SetTitleOffset(0.7);
h_CVCF_deltaNLL->GetYaxis()->SetTitleOffset(0.7);

h_CVCF_deltaNLL->GetXaxis()->SetRangeUser(0.,3.0);

TCanvas *c2 = new TCanvas("c2","c2",600,600);
h_CVCF_deltaNLL->SetStats(0);
h_CVCF_deltaNLL->Draw("axis");
gr68->Draw("LF SAME");
gr95->Draw("LF SAME");
gr68->Draw("LF SAME");
gr_SM->Draw("P SAME");
gr_best->Draw("P SAME");

legend->Draw();

c2->SetTopMargin(0.08);
c2->SetRightMargin(0.08);

c2->Print("Images/Images_2013_08_27_ttH_couplings/kV_kf_bb_tt_gg.png");
c2->Print("Images/Images_2013_08_27_ttH_couplings/kV_kf_bb_tt_gg.pdf");

c2->Print("Images/Images_2013_08_27_ttH_couplings/kV_kf_bb_tt.png");
c2->Print("Images/Images_2013_08_27_ttH_couplings/kV_kf_bb_tt.pdf");

TCanvas *c2 = new TCanvas("c2","c2");
h_CVCF_deltaNLL->SetStats(0);
h_CVCF_deltaNLL->Draw("axis");
gr68->Draw("LF SAME");
gr95->Draw("LF SAME");
gr68->Draw("LF SAME");
gr_best->Draw("P SAME");


TCanvas *c3 = new TCanvas("c3","c3");
h_CVCF_deltaNLL_full->SetMaximum(4);
h_CVCF_deltaNLL_full->Draw("colz");


TCanvas *c4 = new TCanvas("c4","c4");
h_CVCF_quantileExpected->Draw("colz");


TCanvas *c5 = new TCanvas("c5","c5");
h_deltaNLL_CV_bestCF->SetStats(0);
h_deltaNLL_CV_bestCF->SetMaximum(3);
h_deltaNLL_CV_bestCF->GetXaxis()->SetRangeUser(-3.2,3.2);
h_deltaNLL_CV_bestCF->Draw();
h_deltaNLL_CV_smCF->SetLineColor(kRed);
h_deltaNLL_CV_smCF->Draw("same");



*/
