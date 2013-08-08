
void rename7TeV(TString fileName) {

  TFile *file = TFile::Open(fileName,"UPDATE");

  //Go looking for histograms with the right names to replace
  TList *keys = file->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();
    TH1 *hist = (TH1 *)file->Get(name);

    //Replace if it matches one of the names we're changing.
    TString replaceString = "Q2scale_ttH_ttbar_bb";
    if (name.Contains(replaceString) && !name.Contains("_7TeV")) {

      TString newName = name;
      newName.Replace(newName.Index(replaceString),replaceString.Length(),
                      replaceString+"_7TeV");


      if (file->Get(newName)) {
        std::cout << "Histogram " << newName 
                  << " already exists.  Doing nothing." << std::endl;
      } else {

        std::cout << "Replacing " << name << " by " << newName << std::endl;
        TH1 *histClone = (TH1 *)hist->Clone(newName);
        histClone->Write();
        file->Delete(name);
      }
    }

    replaceString = "Q2scale_ttH_ttbar_cc";
    if (name.Contains(replaceString) && !name.Contains("_7TeV")) {

      TString newName = name;
      newName.Replace(newName.Index(replaceString),replaceString.Length(),
                      replaceString+"_7TeV");


      if (file->Get(newName)) {
        std::cout << "Histogram " << newName 
                  << " already exists.  Doing nothing." << std::endl;
      } else {

        std::cout << "Replacing " << name << " by " << newName << std::endl;
        TH1 *histClone = (TH1 *)hist->Clone(newName);
        histClone->Write();
        file->Delete(name);
      }
    }

    replaceString = "CMS_scale_j";
    if (name.Contains(replaceString) && !name.Contains("_7TeV")) {

      TString newName = name;
      newName.Replace(newName.Index(replaceString),replaceString.Length(),
                      replaceString+"_7TeV");


      if (file->Get(newName)) {
        std::cout << "Histogram " << newName 
                  << " already exists.  Doing nothing." << std::endl;
      } else {

        std::cout << "Replacing " << name << " by " << newName << std::endl;
        TH1 *histClone = (TH1 *)hist->Clone(newName);
        histClone->Write();
        file->Delete(name);
      }
    }

  }

}
