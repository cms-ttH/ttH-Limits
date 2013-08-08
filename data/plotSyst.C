class CategoryInfo : public TObject {

public:
  TString name;
  int nRebin;
  TString label;

  CategoryInfo(TString n, TString l, int nb) :
    name(n), label(l), nRebin(nb) {}

};

class SampleInfo : public TObject {

public:
  TString name;
  
  TString label;

  SampleInfo(TString n, TString l) :
    name(n), label(l) {}

};


void plotSyst(TString systName, TString systLabel, TString fileName = "AN_v1_LJ_DIL.root") {

  gStyle->SetOptStat(0);

  TFile *histFile = TFile::Open(fileName);

  TList samples;
//   samples.Add(new SampleInfo("ttbarPlusB","t#bar{t}+b"));
   samples.Add(new SampleInfo("ttbarPlusBBbar","t#bar{t}+b#bar{b}"));
//  samples.Add(new SampleInfo("ttbar","t#bar{t}+LF"));
//   samples.Add(new SampleInfo("ttH125","t#bar{t}H"));


  TList categories;
//   categories.Add(new CategoryInfo("ljets_j4_t2", "Lepton+Jets, 4 jet, 2 tag",1));
//  categories.Add(new CategoryInfo("ljets_j4_t3", "Lepton+Jets, 4 jet, 3 tag",1));
//  categories.Add(new CategoryInfo("ljets_j4_t4", "Lepton+Jets, 4 jet, 4 tag",1));
//   categories.Add(new CategoryInfo("ljets_j5_t2", "Lepton+Jets, 5 jet, 2 tag",1));
//  categories.Add(new CategoryInfo("ljets_j5_t3", "Lepton+Jets, 5 jet, 3 tag",1));
//  categories.Add(new CategoryInfo("ljets_j5_tge4", "Lepton+Jets, 5 jet, #geq 4 tag",1));
//  categories.Add(new CategoryInfo("ljets_jge6_t2", "Lepton+Jets, #geq 6 jet, 2 tag",1));
//  categories.Add(new CategoryInfo("ljets_jge6_t3", "Lepton+Jets, #geq 6 jet, 3 tag",1));
  categories.Add(new CategoryInfo("ljets_jge6_tge4", "Lepton+Jets, #geq 6 jet, #geq 4 tag",1));
  categories.Add(new CategoryInfo("ge3t", "Dilepton, #geq 3 tag",1));
//   categories.Add(new CategoryInfo("e3je2t", "Dilepton, 3 jet, 2 tag",1));
//   categories.Add(new CategoryInfo("ge4je2t", "Dilepton, 3 jet, 2 tag",1));

//   categories.Add(new TObjString("e3je2t", "Dilepton, 3 jet, 2 tag",1));
//   categories.Add(new TObjString("ge4je2t", "Dilepton, #geq 4 jet, 2 tag",1));


  TIter nextCat(&categories);
  CategoryInfo *cat = 0;

  TString summaryText = "";
  summaryText += "Sample/Category                           Up        Down\n";
  summaryText += "----------------------------------------  ----------  ----------\n";
  
  while ((cat = (CategoryInfo *)nextCat())) {

    std::cout << "==>Doing cat " << cat->name << std::endl;

    //Loop over samples, make plots and add to legend
    TIter nextSample(&samples);
    SampleInfo *sample = 0;

    while ((sample = (SampleInfo *)nextSample())) {

      std::cout << "Doing sample " << sample->name << std::endl;

      //Central sample
      TString histName = sample->name + "_MVA_" + cat->name;
      std::cout << "  ==>Extracting histogram " << histName << std::endl;
      TH1 *hist = (TH1 *)histFile->Get(histName);
      std::cout << "    ==> nBins = " << hist->GetNbinsX() << " , GetName() = " << hist->GetName() << std::endl;
      hist->Rebin(cat->nRebin);
      double norm = hist->Integral();
      hist->Scale(1./hist->Integral());
      hist->SetTitle(systLabel+" "+sample->label+" "+cat->label);
      hist->GetYaxis()->SetTitle("Unit Normalized");
      hist->SetTitle(cat->label);
      hist->GetXaxis()->SetTitle("ANN Output");

      hist->SetLineColor(kBlack);
      hist->SetFillColor(0);
      hist->SetLineWidth(2);

      //Syst Up sample
      TString histNameUp = histName + "_";
      histNameUp += systName;
      histNameUp += "Up";
      std::cout << "  ==>Extracting histogram " << histNameUp << std::endl;
      TH1 *histUp = (TH1 *)histFile->Get(histNameUp);
      std::cout << "    ==> nBins = " << histUp->GetNbinsX() << " , GetName() = " << histUp->GetName() << std::endl;
      histUp->Rebin(cat->nRebin);
      double normUp = histUp->Integral();
      histUp->Scale(1./histUp->Integral());

      histUp->SetTitle(systLabel+" "+sample->label+" "+cat->label);
      histUp->SetLineColor(kRed);
      histUp->SetFillColor(0);
      histUp->SetLineWidth(2);

      //Syst Down sample
      TString histNameDown = histName + "_";
      histNameDown += systName;
      histNameDown += "Down";
      std::cout << "  ==>Extracting histogram " << histNameDown << std::endl;
      TH1 *histDown = (TH1 *)histFile->Get(histNameDown);
      std::cout << "    ==> nBins = " << histDown->GetNbinsX() << " , GetName() = " << histDown->GetName() << std::endl;
      histDown->Rebin(cat->nRebin);
      double normDown = histDown->Integral();
      histDown->Scale(1./histDown->Integral());

      histDown->SetTitle(systLabel+" "+sample->label+" "+cat->label);
      histDown->SetLineColor(kBlue);
      histDown->SetFillColor(0);
      histDown->SetLineWidth(2);


      //Calculate rate shift
      double shiftUp = normUp/norm - 1;
      double shiftDown = normDown/norm -1;

      std::cout << "    --> norm = " << norm << std::endl
                << "    --> normUp = " << normUp << std::endl
                << "    --> normDown = " << normDown << std::endl
                << "    --> shiftUp = " << shiftUp << std::endl
                << "    --> shiftDown = " << shiftDown << std::endl;

      TString summaryLabel = sample->name + " " + cat->name;
      summaryText += Form("%40s  %9.2f%%  %9.2f%%\n",summaryLabel.Data(),100.*shiftUp,100.*shiftDown);


      double histMax = TMath::Max(hist->GetMaximum(),
                                  histUp->GetMaximum());

      histMax = TMath::Max(histMax,
                           histDown->GetMaximum());

      hist->SetMaximum(1.3*histMax);
      histUp->SetMaximum(1.3*histMax);
      histDown->SetMaximum(1.3*histMax);

      //Make a canvas for plotting
      TCanvas *canv = new TCanvas("SystPlot_"+systName+"_"+
                                  sample->name+"_"+cat->name,
                                  "SBPlot_"+systName+"_"
                                  +sample->name+"_"+cat->name,500,500);
      canv->cd();
    
      TLegend *leg = new TLegend(0.11,0.8,0.89,0.89);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetNColumns(3);
      leg->AddEntry(hist,"Central","l");
      leg->AddEntry(histUp,"Up","l");
      leg->AddEntry(histDown,"Down","l");
            
      histUp->Draw("HIST");
      histDown->Draw("SAMEHIST");
      hist->Draw("SAMEHIST");
      leg->Draw();
    
      canv->Print(".pdf");



    }
  }

  std::cout << summaryText << std::endl;

}
        



    
