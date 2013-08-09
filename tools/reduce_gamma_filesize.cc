// Script to reduce the filesize of the gamma gamma limit file:
//
// Usage: compile and run in appropriate directory
//
// Contains hard-coded paths to the input and output files;  Change these
// accordingly.  Also has hard-coded category and mass-point filters.
//
// Compile with: c++ $(root-config --cflags --ldflags --libs) -I/pscratch/osg/app/cmssoft/cms/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1/external/slc5_amd64_gcc472/bin/../../../../../../lcg/roofit/5.34.03-cms4/include -L/pscratch/osg/app/cmssoft/cms/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1/external/slc5_amd64_gcc472/bin/../../../../../../lcg/roofit/5.34.03-cms4/lib -lRooFit -lRooStats -lRooFitCore reduce_gamma_filesize.cc -o reduce_gamma_filesize
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"

int
main()
{
   TFile infile("orig.root");
   TFile outfile("CMS-HGG_sigdata_interpolated_interpolated.root", "RECREATE");

   RooWorkspace *inws, *ws;

   infile.GetObject("cms_hgg_workspace", inws);
   ws = new RooWorkspace("cms_hgg_workspace");

   for (const auto& data: inws->allData()) {
      std::string name = data->GetName();
      if (name.find("cat6") != std::string::npos || name.find("cat7") != std::string::npos) {
         auto pos = name.find("mass_m1");
         // Filter out additional mass points (keep 5 GeV steps) - but only do
         // so where needed (first check)
         if (pos == std::string::npos || (name[pos + 9] != '.' && (name[pos + 8] == '5' || name[pos + 8] == '0')))
            ws->import(*data);
      }
   }

   auto set = inws->allPdfs();
   auto it = set.createIterator();
   while (TObject *o = it->Next()) {
      auto pdf = dynamic_cast<RooAbsPdf*>(o);
      if (pdf) {
         std::string name = pdf->GetName();
         if (name.find("cat6") != std::string::npos || name.find("cat7") != std::string::npos)
            ws->import(*pdf);
      }
   }

   ws->Write();
   infile.Close();
}
