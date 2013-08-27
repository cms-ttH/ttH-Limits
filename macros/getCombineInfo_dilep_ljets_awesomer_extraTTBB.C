 //Root includes                                   
 #include "TROOT.h"
 #include "Riostream.h"
 #include "TFile.h"
 #include "TH1.h"
 #include "TF1.h"
 #include "TH1F.h"
 #include "TSystem.h"
 #include "TStyle.h"
 #include "TTree.h"
 #include "TString.h"
 #include "TMath.h"
 #include "TAxis.h"
 #include "TKey.h"
 #include "TList.h"

//
//
// To use:
// root [0] .L getCombineInfo_dilep_ljets_new.C

//
// To generate the "combination" datacard with all channels enabled, call like this:
//     root [1] getCombineInfo_dilep_ljets("ttH120","combined_ljets_dilep","result_combined.root",7,7,7,0,0,0,1,1);
//     ===> the outout data card is then: combined_ljets_dilep_mu777_el777_dilep3t1_dilep2t1.dat
//
// To generate a datacard with all channels enabled, call like this:
//     root [1] getCombineInfo_dilep_ljets("ttH120","combined_ljets_dilep","result_combined.root",7,7,7,7,7,7,1,1);
//     ===> the outout data card is then: combined_ljets_dilep_mu777_el777_dilep3t1_dilep2t1.dat
//
// To generate a datacard with just l+jets channels enabled, call like this:
//     root [1] getCombineInfo_dilep_ljets("ttH120","combined_ljets_dilep","combined_ljets_dilep.root",7,7,7,7,7,7,0,0);
//     ===> the outout data card is then: combined_ljets_dilep_mu777_el777_dilep3t0_dilep2t0.dat
//
// To generate a datacard with just 6 jet channels enabled, call like this:
//     root [1] getCombineInfo_dilep_ljets("ttH120","combined_ljets_dilep","combined_ljets_dilep.root",0,0,7,0,0,7,0,0);
//     ===> the outout data card is then: combined_ljets_dilep_mu007_el007_dilep3t0_dilep2t0.dat
//
//
//
void getCombineInfo_dilep_ljets_awesomer_extraTTBB(TString hmass="ttH125",TString outputTextFileName="dataCardOutput_ttH_ljets_dilep_120",
						   TString histName = "result_combined.root",
						   int ljets4cat=7,int ljets5cat=7,int ljets6cat=7,
                                                   int dilepCat1 = 1, int dilepCat2 = 7,int SS1tCat = 3, int SS2tCat = 3, bool is8TeV = true, TString ttbb_unc="1.5",
                                                   TString NPSFValueString  = "1.5",
                                                   bool noNorm = true,
						   int systFlagEnable = 0x3FFFFF)
{
//
//	TString discName = "MLP";
	TString discName = "CFMlpANN";
	const int numChannels = 17;
	TString channelNames[] = {
										"ljets_j4_t2", 		//	0 "j4_2t",
										"ljets_j4_t3", 		//	1 "j4_3t",
										"ljets_j4_t4", 		//	2 "j4_3t",
										"ljets_j5_t2", 		//	3 "j5_2t",
										"ljets_j5_t3", 		//	4 "j5_3t",
										"ljets_j5_tge4", 		//	5 "j5_3t",
										"ljets_jge6_t2", 		//	6 "j6_2t",
										"ljets_jge6_t3", 		//	7 "j6_3t",
										"ljets_jge6_tge4", 	//	8 "j6_3t",
										"ge3t",				// 9 dilep >= 3tags  cat1 bit0
										"e2je2t",				// 10 dilep 2 jets 2 tags cat2 bit0
										"e3je2t",           // 11 dilep 3 jets, 2 tags cat2 bit1
										"ge4je2t",          // 12 dilep >=4 jets, 2 tags cat2 bit2
										"SS_e3je1t",   		// 13 dilep SS 3 jets 1 tags SScat1 bit0
										"SS_ge4je1t",   	// 14 dilep SS >=4 jets 1 tags SScat1 bit1
										"SS_e3jge2t",   	// 15 dilep SS 3 jets >=2 tags SScat2 bit0
										"SS_ge4jge2t",   	// 16 dilep SS >=4 jets >=2 tags SScat2 bit1
										};
	TString channelNamesOutput[] = {
										"ljets_j4_t2", 		//	 "j4_2t",
										"ljets_j4_t3", 		//	 "j4_3t",
										"ljets_j4_t4", 		//	 "j4_3t",
										"ljets_j5_t2", 		//	 "j5_2t",
										"ljets_j5_t3", 		//	 "j5_3t",
										"ljets_j5_tge4", 		//	 "j5_3t",
										"ljets_jge6_t2", 		//	 "j6_2t",
										"ljets_jge6_t3", 		//	 "j6_3t",
										"ljets_jge6_tge4", 	//	 "j6_3t",
										"ge3t",			// dilep >= 3tags
										"e2je2t",		// dilep 2 jets 2 tags
										"e3je2t",           // dilep 3 jets, 2 tags
										"ge4je2t",          // dilep >=4 jets, 2 tags
										"SS_e3je1t",   		// dilep SS 3 jets 1 tags 
										"SS_ge4je1t",   	// dilep SS >=4 jets 1 tags
										"SS_e3jge2t",   	// dilep SS 3 jets >=2 tags
										"SS_ge4jge2t",   	// dilep SS >=4 jets >=2 tags
										};
	int fileNumber[] = {
							0,0,0,0,0,0,0,0,0,
							1,1,1,1,2,2,2,2};				// last 4 are from dilepton
	
//
// Show which categories are included
	bool channelEnable[numChannels];
	for (int i=0; i<numChannels; i++) {
		channelEnable[i] = false;
	}
	for (int i=0; i<3; i++){
		if (ljets4cat & (1<< i)) channelEnable[i] = true;
		if (ljets5cat & (1<< i)) channelEnable[i+3] = true;
		if (ljets6cat & (1<< i)) channelEnable[i+6] = true;
	}
	if (dilepCat1 & 1) channelEnable[9] = true;
	if (dilepCat2 & 1) channelEnable[10] = true;
	if (dilepCat2 & 2) channelEnable[11] = true;
	if (dilepCat2 & 4) channelEnable[12] = true;

	///// DIL SS
	if (SS1tCat & 1) channelEnable[13] = true;
	if (SS1tCat & 2) channelEnable[14] = true;
	if (SS2tCat & 1) channelEnable[15] = true;
	if (SS2tCat & 2) channelEnable[16] = true;

	for (int i=0; i<numChannels; i++) {
		cout << channelNames[i] << ": enabled? " << channelEnable[i] << endl;
	}
	
	char catVals[50];     // "_muXXX_elXXX.dat"
	int nnn = sprintf(catVals,"_ljet%d%d%d_dilep3t%d_dilep2t%d_SS1t%d_SS2t%d.dat",ljets4cat,ljets5cat,ljets6cat,dilepCat1,dilepCat2,SS1tCat,SS2tCat);
	outputTextFileName.Append(catVals);
    outputTextFileName = "dataCardOutput/" + outputTextFileName;
	cout << "File name = " << outputTextFileName << endl;

    //exit (10);
//
// Setup the output file
	ofstream of;
	of.open(outputTextFileName.Data());
// 
	int numSamples = 10;
	TString sampleNames[] = {"ttH125","ttbar","ttbarPlusBBbar","ttbarPlusCCbar","singlet","wjets","zjets","ttbarW","ttbarZ","diboson"};
	sampleNames[0] = hmass;
	TString signalName = "ttH";

//	TString histName = "multinet_reh.root";
	
	TFile* inputFile = new TFile(histName);
	//TFile* inputfile[2];
//	inputFile[0] = new TFile("ttH_ljets_scale_7TeV.root");
//	inputFile[1] = new TFile("result_combined.root");
	//inputFile[0] = new TFile("ttH_7TeV.root");
	//inputFile[1] = new TFile("ttH_7TeV.root");
	
	//inputFile->ls();
//
// Print out header info
//	of << "imax " << numChannels << " # number of channels" << endl;
//	of << "jmax " << numSamples-1 << " # number of backgrounds" << endl;
//	of << "kmax " << "*" << " # number of nuisance parameters" << endl;
//	of << "---------------" << endl;

	of << "imax * # number of channels" << endl;
	of << "jmax * # number of backgrounds" << endl;
	of << "kmax * # number of nuisance parameters" << endl;
	of << "---------------" << endl;

//
// Do observed data
//
// Loop over output lines
	for (int n=0; n<2; n++) {
		if (n==0) {
			of << "bin                ";
		} else {
			of << "observation       ";
		}
		
//
// get the data
		TString hhname = "data_obs";
		hhname.Append("_");
		hhname.Append(discName);
		hhname.Append("_");
//
// Loop over the channels for each sample
		for (int j=0; j<numChannels; j++) {
			TString hname = hhname;
			hname.Append(channelNames[j]);
			if (n==0) {
				if (channelEnable[j]) of << " " << channelNamesOutput[j];
//				of << setw(10) << channelNames[j];
			} else {
              if (channelEnable[j]) {
				cout << "Getting hist " << hname << " " << fileNumber[j] << endl;
				TH1F* th = (TH1F*)inputFile->Get(hname);
//				TH1F* th = (TH1F*)inputFile[fileNumber[j]]->Get(hname);
				float num = th->Integral();
				cout << "data " << channelNames[j] << " " << num << endl;
//				of.setf(ios::fixed, ios::floatfield);
//				of << setprecision(0) << setw(12) << num;
				of << "   " << Form("%f",num);
                              }
			}
		}
		of << endl;
	}  	
	of << "---------------" << endl;
//
// 
	of << "shapes * * " << histName << " $PROCESS_CFMlpANN_$CHANNEL $PROCESS_CFMlpANN_$CHANNEL_$SYSTEMATIC" << endl;
	of << "shapes ttH * " << histName << " $PROCESS$MASS_CFMlpANN_$CHANNEL $PROCESS$MASS_CFMlpANN_$CHANNEL_$SYSTEMATIC" << endl;
	of << "---------------" << endl;
//
// NOw do backgrounds
//
// Loop over output lines
	for (int n=0; n<4; n++) {
		if (n==0) {
			of << "bin        ";
		} else if ((n==1) || (n==2)) {
			of << "process    ";
		} else if (n==3) {
			of << "rate       ";
		}
		
//
// Loop over the samples
		for (int i=0; i<numSamples; i++) {
			TString hhname = sampleNames[i];
			hhname.Append("_");
			hhname.Append(discName);
			hhname.Append("_");
//
// Loop over the channels for each sample
			for (int j=0; j<numChannels; j++) {
				TString hname = hhname;
				hname.Append(channelNames[j]);
				if (n==0) {
//					of << setw(12) << channelNames[j];
					if (channelEnable[j]) of << " " << channelNamesOutput[j];
				} else if (n==1) {
//					of << setw(12) << sampleNames[i];
					if (i == 0) {
						if (channelEnable[j]) of << " " << signalName;
					} else {
						if (channelEnable[j]) of << " " << sampleNames[i];
					}
				} else if (n==2) {
//					of << setw(12) << i;
					if (channelEnable[j]) of << " " << i;
				} else if (n==3) {
                                  if (channelEnable[j]) {
                                    cout << "hanme " << hname << endl;
                                    TH1F* th = (TH1F*)inputFile->Get(hname);
                                    //TH1F* th = (TH1F*)inputFile[fileNumber[j]]->Get(hname);
                                   float num = th->Integral();
                                   //of.setf(ios::fixed, ios::floatfield);
                                    //					of << setprecision(2) << setw(12) << num;
                                   if (noNorm && i == 0) {
                                     of << "  -1";
                                   } else {
                                     of << "  " << num;
                                   }
                                  }
				}
			}
		}
		of << endl;
	}  	
	of << "---------------" << endl;
//
// Loop over the samples
	for (int j=0; j<numChannels; j++) {
		if (channelEnable[j]) {
			cout << "channel " << channelNames[j] << " ";
			float sum = 0.0;
			for (int i=0; i<numSamples; i++) {
				TString hhname = sampleNames[i];
				hhname.Append("_");
				hhname.Append(discName);
				hhname.Append("_");
				TString hname = hhname;
				hname.Append(channelNames[j]);
				TH1F* th = (TH1F*)inputFile->Get(hname);
//				TH1F* th = (TH1F*)inputFile[fileNumber[j]]->Get(hname);
				float num = th->Integral();
				if (i == 0) cout << "  ttH: " << num;
				if (i == 1) cout << "; summed background: ";
				if (i >= 1) sum = sum + num;
			}
			cout << sum << endl;
		}
	}
//
// Loop over the channels for each sample

	int numSysts = 55;//21;
	TString systNames[] = {
	  "lumi",
	  "QCDscale_ttH",
	  "QCDscale_ttbar",
	  "pdf_gg",
	  "pdf_qqbar",
	  "pdf_qg",
	  "QCDscale_V",
	  "QCDscale_VV",
	  "CMS_ttH_pu",
	  "CMS_res_j",
	  "CMS_ttH_eff_lep",
	  "Q2scale_ttH_ttbar0p",
	  "Q2scale_ttH_ttbar_bb",
	  "Q2scale_ttH_ttbar_cc",
	  "Q2scale_ttH_V",
	  "CMS_fake_bRate",
	  "CMS_eff_bRate",
	  "CMS_scale_j",
	  "CMS_ttH_topPtcorr",
	  "ljets_j4_t3_fake_bShape",
	  "ljets_j4_t4_fake_bShape",
	  "ljets_j5_t3_fake_bShape",
	  "ljets_j5_tge4_fake_bShape",
	  "ljets_jge6_t2_fake_bShape",
	  "ljets_jge6_t3_fake_bShape",
	  "ljets_jge6_tge4_fake_bShape",
	  "ge3t_fake_bShape",
	  "e2je2t_fake_bShape",
	  "e3je2t_fake_bShape",
	  "ge4je2t_fake_bShape",
	  "SS_e3je1t_fake_bShape",
	  "SS_ge4je1t_fake_bShape",
	  "SS_e3jge2t_fake_bShape",
	  "SS_ge4jge2t_fake_bShape",
	  "ljets_j4_t3_eff_bShape",
	  "ljets_j4_t4_eff_bShape",
	  "ljets_j5_t3_eff_bShape",
	  "ljets_j5_tge4_eff_bShape",
	  "ljets_jge6_t2_eff_bShape",
	  "ljets_jge6_t3_eff_bShape",
	  "ljets_jge6_tge4_eff_bShape",
	  "ge3t_eff_bShape",
	  "e2je2t_eff_bShape",
	  "e3je2t_eff_bShape",
	  "ge4je2t_eff_bShape",
	  "SS_e3je1t_eff_bShape",
	  "SS_ge4je1t_eff_bShape",
	  "SS_e3jge2t_eff_bShape",
	  "SS_ge4jge2t_eff_bShape",
	  "Q2scale_ttH_ttbar1p",
	  "Q2scale_ttH_ttbar2p",
	  "CMS_ttH_QCDscale_ttbb",
      "SS_NPSF_4j1t",
      "SS_NPSF_3j2t",
      "SS_NPSF_4j2t"
	};

	TString sysTypes[] = {
										"lnN",			 //	"lumi",
										"lnN",			 //	"QCDscale_ttH",
										"lnN",			 //	"QCDscale_ttbar",
										"lnN",			 //	"pdf_gg",
										"lnN",			 //	"pdf_qqbar",
										"lnN",			 //	"pdf_qg",
										"lnN",			 //	"QCDscale_V",
										"lnN",			 //	"QCDscale_VV",
										"lnN",			 //	"CMS_ttH_pu",
										"lnN",			 //   "CMS_res_j",
										"lnN",			 //   "CMS_ttH_eff_lep",
										"shape",			 //   "Q2scale_ttH_ttbar",
										"shape",    	 //   "Q2scale_ttH_ttbar_bb",
										"shape", 		 //	"Q2scale_ttH_ttbar_cc",
										"lnN", 		 //	"Q2scale_ttH_V",
										"shape", 		 //	"CMS_fake_bRate",
										"shape",			 //	"CMS_eff_bRate",
										"shape",			 //	"CMS_scale_j",
										"shape",              // "CMS_ttH_PUcorr"
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape", 		 //	"CMS_fake_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",			 //	"CMS_eff_bShape",
										"shape",              // 1p
										"shape",              // 2p
										"lnN",			 //	"CMS_ttH_QCDscale_ttbb"
                                        "lnN",          // NPSF 4j1t
                                        "lnN",          // NPSF 3j2t
                                        "lnN"          //  NPSF 4j2t
        };									


    
	TString lumi_unc = ( is8TeV ) ? "1.044" : "1.022";
	TString eff_lep_unc = ( is8TeV ) ? "1.04" : "1.018";

//        "ttH120", "ttbar","ttbarPlusBBbar","ttbarPlusCCbar", "singlet", "wjets", "zjets", "ttbarW", "ttbarZ", "diboson"};
	TString valSyst[] = {  
	  lumi_unc, lumi_unc,    lumi_unc,  	 lumi_unc, 	 	  lumi_unc,       lumi_unc,     lumi_unc,  lumi_unc,       lumi_unc,      lumi_unc, // "lumi",
	  "1.125",  "-",   	 "-", 		 "-", 	 	 	  "-",		  "-",		"-",  	   "-",  	   "-",  	  "-",     // "QCDscale_ttH",
	  "-",      "1.12",      "1.12",         "1.12",      	          "1.02",	  "-",		"-",  	   "1.15",  	   "1.15",  	  "-",     // "QCDscale_ttbar",
	  "1.08",   "1.09",      "1.09",   	 "1.09", 	    	  "-",		  "-",		"-",  	   "-",  	   "1.09",  	  "-",     // "pdf_gg",
	  "-",      "-",         "-", 		 "-",  	    	          "-",	          "1.048",	"1.042",   "1.07",         "-",		  "-",     // "pdf_qqbar",
	  "-",      "-",         "-", 		 "-",  		 	  "1.046",	  "-",   	"-",	   "-",	           "-",		  "-",     // "pdf_qg",
	  "-",      "-",   	 "-", 		 "-",     	 	  "-",	          "1.013",	"1.012",   "-",            "-", 	  "-",     // "QCDscale_V",
	  "-",      "-",   	 "-", 		 "-", 	 	 	  "-",		  "-",		"-", 	   "-",  	   "-",  	  "1.035", // "QCDscale_VV",
	  "1.01",   "1.01",	 "1.01", 	 "1.01",     	          "1.01",         "1.01",       "1.01",    "1.01",         "1.01",	  "1.01",  // "CMS_ttH_pu",
	  "1.015",  "1.015",	 "1.015",	 "1.015",    	          "1.015",        "1.015",      "1.015",   "1.015",        "1.015",	  "1.015", // "CMS_res_j",
	  eff_lep_unc,   eff_lep_unc,	 eff_lep_unc,	 eff_lep_unc,     eff_lep_unc,    eff_lep_unc,  eff_lep_unc, eff_lep_unc, eff_lep_unc,	  eff_lep_unc, // "CMS_ttH_eff_lep",
	  "-",      "1",   	 "-", 		 "-", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "Q2scale_ttH_ttbar",
	  "-",      "-",   	 "1", 		 "-", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "Q2scale_ttH_ttbar_bb",
	  "-",      "-",   	 "-", 		 "1", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "Q2scale_ttH_ttbar_cc",
	  "-",      "-",   	 "-", 		 "-", 	 	 	  "-",		  "-",	        "-", 	   "-",  	   "-",  	  "-",     // "Q2scale_ttH_V" (value set below),
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bRate",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bRate",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "-",   	"1", 	   "1", 	   "1",  	  "1",	     // "CMS_scale_j",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",	     // CMS_ttH_PUcorr,
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_fake_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "1",      "1",   	 "1", 		 "1", 	 	 	  "1",		  "1",   	"1", 	   "1", 	   "1",  	  "1",     // "CMS_eff_bShape",
	  "-",      "1",   	 "-", 		 "-", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "Q2scale_ttH_ttbar1p",
	  "-",      "1",   	 "-", 		 "-", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "Q2scale_ttH_ttbar2p",
	  "-",      "-",   	 ttbb_unc, 	 "-", 	 	 	  "-",		  "-",   	"-", 	   "-", 	   "-",  	  "-",     // "CMS_ttH_QCDscale_ttbb",
      "-",      NPSFValueString,   	 NPSFValueString, 	 NPSFValueString,   NPSFValueString,		  NPSFValueString,   	NPSFValueString, 	   "-", 	   "-",  	  "-",     // NPSF 4j1t
      "-",      NPSFValueString,   	 NPSFValueString, 	 NPSFValueString, 	NPSFValueString,		  NPSFValueString,   	NPSFValueString, 	   "-", 	   "-",  	  "-",     // NPSF 3j2t
      "-",      NPSFValueString,   	 NPSFValueString, 	 NPSFValueString, 	NPSFValueString,		  NPSFValueString,   	NPSFValueString, 	   "-", 	   "-",  	  "-"     //  NPSF 4j2t
        };
	
        
	bool systEnabled[] = {
								true,	  //  "lumi",
								true,	  //  "QCDscale_ttH",
								true,	  //  "QCDscale_ttbar",
								true,	  //  "pdf_gg", 												
								true,	  //  "pdf_qqbar", 												
								true,	  //  "pdf_qg", 												
								true,	  //  "QCDscale_V",  										
								true,	  //  "QCDscale_VV",
								true,	  //  "CMS_ttH_pu",
								true,	  //  "CMS_res_j",
								true,	  //  "CMS_ttH_eff_lep",
								true,   //  "Q2scale_ttH_ttbar",
								true,   //  "Q2scale_ttH_ttbar_bb",
								true,	  //  "Q2scale_ttH_ttbar_cc",
								true,	  //  "Q2scale_ttH_V",
								true,	  //  "CMS_fake_bRate",
								true,   //  "CMS_eff_bRate",
								true,	  //  "CMS_scale_j",
								true,      // CMS_ttH_PUcorr
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,	  //  "CMS_fake_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   //  "CMS_eff_bShape",
								true,   // 1p
								true,   // 2p
								true,     // "CMS_ttH_QCDscale_ttbb"
                                true,   // NPSF 4j1t
                                true,   // NPSF 3j2t
                                true   // NPSF 4j2t 
								}; 			
								

//
// Loop over systematica
	for (int n=0; n<numSysts; n++) {
          // if (!(systFlagEnable & (1<<n)))
          //   systEnabled[n] = false;
          printf("setting enables: %d %x %x; %d\n",n,systFlagEnable,(1<<n),systEnabled[n]);
	}							
    //
    // Loop over systematica
	for (int n=0; n<numSysts; n++) {
      if (systEnabled[n]) {

	TString eraSuffix = ( is8TeV ) ? "_8TeV" : "_7TeV";
	TString systName = systNames[n];

	if( systName.Contains("lumi") ) systName += eraSuffix;
	
        of <<	systName << "  " << 	sysTypes[n] << "  ";

        //
        // Loop over the samples
        for (int i=0; i<numSamples; i++) {
          //
          // Loop over the channels for each sample
          for (int j=0; j<numChannels; j++) {
		
            if (channelEnable[j]) {
                  
              const char *syst = valSyst[n*numSamples+i];
                  
              if( systNames[n].Contains("bShape") ){
                if( systNames[n].Contains(channelNames[j]) ) syst = Form("1");
                else                                         syst = Form("-");
              }


                                              
              // Q2_ttH_ttbar0p only for ljets_j4_t3, ljets_j4_t4, e2je2t
              if (systNames[n] == "Q2scale_ttH_ttbar0p"){
                // if you're not any of these, then skip
                if (!channelNames[j].Contains("e2je2t")
                    && !channelNames[j].Contains("SS_e3je1t")
                    && !channelNames[j].Contains("SS_e3jge2t")
                    && !channelNames[j].Contains("ljets_j4_t2")
                    && !channelNames[j].Contains("ljets_j4_t3")
                    && !channelNames[j].Contains("ljets_j4_t4")) {                      
                  syst = Form("-");
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }
                  
              if (systNames[n] == "Q2scale_ttH_ttbar1p"){
                // if you're not any of these, then skip
                if (!channelNames[j].Contains("ge3t") && channelNames[j]!="e3je2t"
                    && !channelNames[j].Contains("SS_ge4je1t")
                    && !channelNames[j].Contains("SS_ge4jge2t")
                    && !channelNames[j].Contains("ljets_j5_t2")
                    && !channelNames[j].Contains("ljets_j5_t3")
                    && !channelNames[j].Contains("ljets_j5_tge4")) {                      
                  syst = Form("-"); // off
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }

              if (systNames[n] == "Q2scale_ttH_ttbar2p"){
                // if you're not any of these, then skip
                if (   !channelNames[j].Contains("ljets_jge6_t2") && channelNames[j]!="ge4je2t"
                       && !channelNames[j].Contains("ljets_jge6_t3")
                       && !channelNames[j].Contains("ljets_jge6_tge4")) {                      
                  syst = Form("-");
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }

              if (systNames[n] == "SS_NPSF_4j1t"){
                // if you're not any of these, then skip
                if (   !channelNames[j].Contains("SS_ge4je1t") ) {                      
                  syst = Form("-");
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }

          
              if (systNames[n] == "SS_NPSF_3j2t"){
                // if you're not any of these, then skip
                if (   !channelNames[j].Contains("SS_e3jge2t") ) {                      
                  syst = Form("-");
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }

              if (systNames[n] == "SS_NPSF_4j2t"){
                // if you're not any of these, then skip
                if (   !channelNames[j].Contains("SS_ge4jge2t") ) {                      
                  syst = Form("-");
                } else {
                  cout << systNames[n] << ": using channel " << channelNames[j]  << "  sample name  " << sampleNames[i] << " value " << syst  <<  endl;
                }
              }



          
              //Make this a little more complicated so that some of the systematics
              //get bigger as the jet bin increases
              if (systNames[n] == "Q2scale_ttH_V") {
                if (sampleNames[i] == "wjets" ||
                    sampleNames[i] == "zjets") {
                  int iJet = 0;
                  if (channelNames[j].Contains("e2je2t")) {
                    iJet = 2;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("e3je2t")) {
                    iJet = 3;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("ge4je2t")) {
                    iJet = 4;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("ge3t")) {
                    iJet = 3;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("e3je1t")) {
                    iJet = 3;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("ge4je1t")) {
                    iJet = 4;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("e3jge2t")) {
                    iJet = 3;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("ge4jge2t")) {
                    iJet = 4;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if (channelNames[j].Contains("j4_t")) {
                    iJet = 4;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if  (channelNames[j].Contains("j5_t")) {
                    iJet = 5;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else if  (channelNames[j].Contains("jge6_t")) {
                    iJet = 6;
                    syst = Form("%f",1+ 0.1 * iJet); //+10% per jet
                  } else {
                    syst = "-";
                  }
                }
              }

	      TString mysyst = syst;

	      if( !(mysyst.EqualTo("1") || mysyst.EqualTo("-")) ){
		double mysyst_val = mysyst.Atof() - 1.;
		double new_mysyst_val = TMath::Exp( TMath::Sqrt( TMath::Log( 1 + mysyst_val*mysyst_val ) ) );
		//mysyst = TString::Form("%.3f",new_mysyst_val);
		//std::cout << 
		//std::cout << systNames[n] << "\t" << mysyst << "\t" << mysyst_val << "\t" << new_mysyst_val << "\t" << Form("%.3f",new_mysyst_val) << std::endl;
		of << Form("%.3f",new_mysyst_val) << "  ";
	      }
	      else of << mysyst << "  ";

              //of << mysyst << "  ";
              //					of << setprecision(2) << setw(12) <<  valSyst[n*numSamples + i] ;
            }      
          }
        }
            
        of << endl;
      } else {
        cout << systNames[n] << " " << 	sysTypes[n] << " not enabled." << endl;
      }
    }   



  ////////////////////
  ///////
  ///////  MC stats in datacard
  ///////
  ////////////////////
	

  for (int iChan = 0; iChan < numChannels; ++iChan) {
    //if( !channelNames[iChan].Contains("ljets_jge6_tge4") ) continue;
    if( channelEnable[iChan] ){

      //continue;
      //Also look at the data
      TString dataHistName = "data_obs_CFMlpANN_";
      dataHistName += channelNames[iChan];
      TH1 *dataHist = (TH1*)inputFile->Get(dataHistName);

      //And why not signal too!
      TString sigName = hmass;
      sigName += "_CFMlpANN_";
      sigName += channelNames[iChan];
      TH1 *sigHist = (TH1*)inputFile->Get(sigName);

      TH1D* hist_sum;

      bool firstTime = true;
      for (int iSample = 1; iSample < numSamples; ++iSample) {

        TString histNamer = sampleNames[iSample];
        histNamer += "_CFMlpANN_";
        histNamer += channelNames[iChan];

        TString histNamer_sum = sampleNames[iSample];
        histNamer_sum += "_CFMlpANN_";
        histNamer_sum += channelNames[iChan] + "_Clone";


        //std::cout << histNamer << " " << std::endl;

        TH1D* hist_temp = (TH1D*)inputFile->Get(histNamer);
      
        if( firstTime ){
          hist_sum = (TH1D*)hist_temp->Clone(histNamer_sum);
          firstTime = false;
        }
        else            hist_sum->Add(hist_temp);
      }


    for (int iSample = 0; iSample < numSamples; ++iSample) {


      //if( !sampleNames[iSample].EqualTo("ttbarPlusBBbar") ) continue;
      //if( !(channelNames[iChan].Contains("tge4") || channelNames[iChan].Contains("t4") || channelNames[iChan].Contains("t3") || channelNames[iChan].Contains("ge3t")) ) continue;
      //if( !(channelNames[iChan].Contains("tge4") || channelNames[iChan].Contains("t4")) ) continue;

      // if( !sampleNames[iSample].EqualTo("ttbar") &&
      //          !sampleNames[iSample].EqualTo("ttbarPlusBBbar")
      //          ) continue;


      TString histNamer = sampleNames[iSample];
      histNamer += "_CFMlpANN_";
      histNamer += channelNames[iChan];

      //std::cout << histNamer << " " << std::endl;

      TH1D* hist = (TH1D*)inputFile->Get(histNamer);


      int numBins = hist->GetNbinsX();

      for( int bin=0; bin<numBins; bin++ ){

        double data = dataHist->GetBinContent(bin+1);
        double dataErr = dataHist->GetBinError(bin+1);

        double sig = sigHist->GetBinContent(bin+1);

        double bkg = hist->GetBinContent(bin+1);
        double err   = hist->GetBinError(bin+1);
        double total_bkg = hist_sum->GetBinContent(bin+1);
        double total_bkg_err = hist_sum->GetBinError(bin+1);

        //Some skipping logic:

        //Skip empty bins
        if (bkg < 0.01) continue;

        //Skip any bin where the total background error is small enough compared to data to be ignorable
        if (total_bkg_err < dataErr/3.) continue;

        //Skip any process that doesn't contribute significantly to the total error
        double totErrSub = sqrt(total_bkg_err*total_bkg_err - err*err);
        if (totErrSub/total_bkg_err > 0.95) continue;

        //Skip any bins where signal is negligible
        if (sig/total_bkg < 0.02) continue;

        //Skip any bin where the total background error is less than 10%
//         if (total_bkg_err/total_bkg < 0.1) continue;

	TString era = ( is8TeV ) ? "8TeV" : "7TeV";

        TString bin_name;
        bin_name.Form("%d",bin);

        TString systName = sampleNames[iSample] + "_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name;
        if( sampleNames[iSample].Contains("ttH") ){
          systName = "ttH_" + channelNames[iChan] + "_" + era + "_ANNbin" + bin_name;
        }

        of <<   systName << "  " <<     "shape" << "  ";

        //
        // Loop over the samples
        for (int i=0; i<numSamples; i++) {
          //
          // Loop over the channels for each sample
          for (int j=0; j<numChannels; j++) {

            if( channelEnable[j] ) {
                  
              const char *syst = ( j==iChan && i==iSample ) ? Form("1") : Form("-");

              of << syst << "  ";
            
            }
          }
        }
        of << endl;
      }
    }
    } 
  }



    of << "---------------" << endl;
       
    inputFile->Close();
    //        inputFile[1]->Close();
    //
    
        
}
