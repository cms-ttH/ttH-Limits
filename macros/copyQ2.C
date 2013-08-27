

{



  TFile *myFile = new TFile ("combined_lj_dil_split3.root", "UPDATE");


  myFile->cd();

  TH1* tempHist;
  ////----------------------------- 0 parton
  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j4_t3_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_j4_t3_Q2scale_ttH_ttbar0pUp");
  tempHist->SetDirectory(myFile);


  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j4_t3_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_j4_t3_Q2scale_ttH_ttbar0pDown");
  tempHist->SetDirectory(myFile);

  
  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j4_t4_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_j4_t4_Q2scale_ttH_ttbar0pUp");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j4_t4_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_j4_t4_Q2scale_ttH_ttbar0pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_e2je2t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_e2je2t_Q2scale_ttH_ttbar0pDown");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_e2je2t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_e2je2t_Q2scale_ttH_ttbar0pUp");
  tempHist->SetDirectory(myFile);
  //// DIL SS
  tempHist = (TH1*) ttbar_CFMlpANN_SS_e3je1t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_SS_e3je1t_Q2scale_ttH_ttbar0pDown");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_SS_e3je1t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_SS_e3je1t_Q2scale_ttH_ttbar0pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_SS_e3jge2t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_SS_e3jge2t_Q2scale_ttH_ttbar0pDown");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_SS_e3jge2t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_SS_e3jge2t_Q2scale_ttH_ttbar0pUp");
  tempHist->SetDirectory(myFile);
  ////----------------------------- 1 parton
  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j5_t3_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_j5_t3_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j5_t3_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_j5_t3_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j5_tge4_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_j5_tge4_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_j5_tge4_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_j5_tge4_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);


  tempHist = (TH1*) ttbar_CFMlpANN_ge3t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ge3t_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ge3t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ge3t_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_e3je2t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_e3je2t_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_e3je2t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_e3je2t_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);
  //// DIL SS
  tempHist = (TH1*) ttbar_CFMlpANN_SS_ge4je1t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_SS_ge4je1t_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_SS_ge4je1t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_SS_ge4je1t_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_SS_ge4jge2t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_SS_ge4jge2t_Q2scale_ttH_ttbar1pDown");
  tempHist->SetDirectory(myFile);
  
  tempHist = (TH1*) ttbar_CFMlpANN_SS_ge4jge2t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_SS_ge4jge2t_Q2scale_ttH_ttbar1pUp");
  tempHist->SetDirectory(myFile);

  
  ////----------------------------- 2 parton
  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_t2_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_jge6_t2_Q2scale_ttH_ttbar2pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_t2_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_jge6_t2_Q2scale_ttH_ttbar2pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_t3_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_jge6_t3_Q2scale_ttH_ttbar2pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_t3_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_jge6_t3_Q2scale_ttH_ttbar2pDown");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_tge4_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ljets_jge6_tge4_Q2scale_ttH_ttbar2pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ljets_jge6_tge4_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ljets_jge6_tge4_Q2scale_ttH_ttbar2pDown");
  tempHist->SetDirectory(myFile);


  tempHist = (TH1*) ttbar_CFMlpANN_ge4je2t_Q2scale_ttH_ttbarUp->Clone("ttbar_CFMlpANN_ge4je2t_Q2scale_ttH_ttbar2pUp");
  tempHist->SetDirectory(myFile);

  tempHist = (TH1*) ttbar_CFMlpANN_ge4je2t_Q2scale_ttH_ttbarDown->Clone("ttbar_CFMlpANN_ge4je2t_Q2scale_ttH_ttbar2pDown");
  tempHist->SetDirectory(myFile);

  myFile->Write();
  myFile->Close();
  
  



}
