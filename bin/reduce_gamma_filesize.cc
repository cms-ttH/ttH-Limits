// Script to reduce the filesize of the gamma gamma limit file:
//
// Usage: compile and run with appropriate input/output filenames
//
// Caveat: hard-coded category and mass-point filters.
#include <iostream>

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"

int
main(int argc, char *argv[])
{
   if (argc != 3) {
      std::cerr << "usage: reduce_gamma_filesize infile outfile\n" << std::endl;
      return 1;
   }

   TFile infile(argv[1]);
   TFile outfile(argv[2], "RECREATE");

   if (!infile.IsOpen()) {
      std::cerr << "could not open input file" << std::endl;
      return 1;
   } else if (!outfile.IsOpen()) {
      std::cerr << "could not open output file" << std::endl;
      return 1;
   }

   RooWorkspace *inws, *ws;

   infile.GetObject("cms_hgg_workspace", inws);
   ws = new RooWorkspace("cms_hgg_workspace");

   for (const auto& data: inws->allData()) {
      std::string name = data->GetName();
      if (name.find("cat6") != std::string::npos || name.find("cat7") != std::string::npos) {
         auto pos = name.find("mass_m1");
         // Let all masses starting with 125 slip through
         auto finegrain = name.find("125");
         // Filter out additional mass points (keep 5 GeV steps) - but only do
         // so where needed (first check)
         if (pos == std::string::npos || finegrain != std::string::npos ||
               (name[pos + 9] != '.' && (name[pos + 8] == '5' || name[pos + 8] == '0')))
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
