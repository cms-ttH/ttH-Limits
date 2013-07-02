


void replaceWithFakeData (TString newFileName, double scaleB = 1.0, double scaleBB = 1.0, double scaleTTH = 1.0) {

  TFile * myFile = new TFile(newFileName, "UPDATE");
  
  TString discName = "MVA";
  const int numChannels = 16;
  TString channelNames[] = {

    "ljets_j4_t3",
    "ljets_j4_t4",
	
    "ljets_j5_t3",
    "ljets_j5_tge4",

    "ljets_jge6_t2",
    "ljets_jge6_t3",
    "ljets_jge6_tge4",

    "ge3t",		
    "e3je2t",
    "ge4je2t",

    "TTL_1b_1nb",
    "TTL_1b_2nb",
    "TTL_1b_3+nb",
    "TTL_2b_0nb",
    "TTL_2b_1nb",
    "TTL_2b_2+nb"
  };

  int numSamples = 11;
  TString sampleNames[] = {"ttH125","ttbar","ttbarPlusB","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson"};

  
  //  first, get the data_obs histogram then delete it
  //  then create a new data_obs histogram and replace it with   
  for ( int iChan = 0; iChan < numChannels; iChan++ ){
    TString dataName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan] + TString(";*");
    TString dataRealName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan];
    TH1 * checkData = (TH1*) myFile->Get(dataRealName);
    double theNorm = checkData->Integral();
    cout << "The data  " << channelNames[iChan]  << "  norm is  " << theNorm  << endl;
    myFile->Delete(dataName);
    //myFile->Write();
    TH1 * sumHisto;
    double channelSumTotal = 0;
    
    for (int iSample = 0; iSample < numSamples; iSample++) {
      TString histoName = sampleNames[iSample] + "_" + discName + "_" +  channelNames[iChan];

      TH1 * tempHisto = (TH1*) myFile->Get(histoName);

      if (! tempHisto) {
        std::cerr << "WARNING: no histogram found for " << histoName << std::endl
                  << "I'm assuming that's OK" << std::endl;
        continue;
      }

      if (histoName.Contains("PlusBBbar")) {
        tempHisto->Scale(scaleBB);
      } else if (histoName.Contains("PlusB")) {
        tempHisto->Scale(scaleB);
      } else if (histoName.BeginsWith("ttH")) {
        tempHisto->Scale(scaleTTH);
      }

      if ( iSample == 0 ){
        TString cloneName =  TString("data_obs_") + discName + TString("_") +  channelNames[iChan];
        sumHisto = (TH1*) tempHisto->Clone(cloneName);
        channelSumTotal += tempHisto->Integral();
        cout << histoName << "   norm =  " << tempHisto->Integral() << ", now sum is " << channelSumTotal << endl;
      } else {
        sumHisto->Add(tempHisto);
        channelSumTotal += tempHisto->Integral();
        cout << histoName << "   norm =  " << tempHisto->Integral() << ", now sum is " << channelSumTotal << endl;
      }
      
    }

    //Now, fix the errors so they look like data
    for (int iBin = 1; iBin <= sumHisto->GetNbinsX(); ++iBin) {

      double n = sumHisto->GetBinContent(iBin);
      sumHisto->SetBinError(iBin,sqrt(n));
      
    }

    cout << "Channel " << channelNames[iChan] << " norm =  " <<  sumHisto->Integral() << endl;
    sumHisto->SetDirectory(myFile);
    sumHisto->Write();

    
  }
  
}
