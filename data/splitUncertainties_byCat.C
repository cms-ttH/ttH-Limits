
void splitUncertainties_byCat(TString inFileName, TString outFileName) {

  int numChannels = 15;
  TString channelNames[] = {
    //"ljets_j4_t2", 		//	 "j4_2t",
    "ljets_j4_t3", 		//	 "j4_3t",
    "ljets_j4_t4", 		//	 "j4_3t",
    //"ljets_j5_t2", 		//	 "j5_2t",
    "ljets_j5_t3", 		//	 "j5_3t",
    "ljets_j5_tge4", 		//	 "j5_3t",
    "ljets_jge6_t2", 		//	 "j6_2t",
    "ljets_jge6_t3", 		//	 "j6_3t",
    "ljets_jge6_tge4",    	//	 "j6_3t",
    "ge3t",	       		// dilep >= 3tags
    "e2je2t",   		// dilep 2 jets 2 tags
    "e3je2t",   		// dilep 3 jets 2 tags
    "ge4je2t",   		// dilep >=4 jets 2 tags
    "SS_e3je1t",   		// dilep SS 3 jets 1 tags
    "SS_ge4je1t",   		// dilep SS >=4 jets 1 tags
    "SS_e3jge2t",   		// dilep SS 3 jets >=2 tags
    "SS_ge4jge2t",   		// dilep SS >=4 jets >=2 tags
  };

  //Open the root file
  TFile *inFile = TFile::Open(inFileName);
  TFile *outFile = TFile::Open(outFileName,"RECREATE");

  //Define the list of uncertainties we're going to split
  TList nuisancePars;
  nuisancePars.Add(new TObjString("CMS_eff_b"));
  nuisancePars.Add(new TObjString("CMS_fake_b"));

  //Now loop over all the histograms in the file and either copy or 
  //Split as needed
  TList *keys = inFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();
    TString newName = name;

//     std::cout << "Copying histogram " << name << std::endl;

    TH1 *hist = (TH1 *)inFile->Get(name);

    if( name.Contains("bShape") ){
      for (int j=0; j<numChannels; j++) {
	if( name.Contains(channelNames[j]) ) newName.ReplaceAll("CMS",channelNames[j]);
      }
    }

    //Automatically copy the original histogram over
    outFile->cd();
    hist->Write(newName);

  }

}
