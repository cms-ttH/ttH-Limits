#!/bin/csh -f

# Nothing fancy here: Assumes you have run a fit and have the
# mlfit.root file lying around.  Runs diffNuisances.py and removed the S+B
# pull info.

set fitFile = mlfit.root

if ( ! -e $fitFile ) then
	echo "Cannot find $fitFile.  Make sure you run the likelihood fit first"
    echo "try combine -M MaxLikelihoodFit -m 125 --rMin -100 --saveWorkspace --minos all datacard.txt"
	exit 2
endif

python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a --stol 9999 --vtol 9999 --stol2 9999 --vtol2 9999 --format latex mlfit.root | awk -F '&' 'BEGIN{print "\\begin{tabular}{|l|r|} \\hline"} NF == 4 {print $1" & "$2" \\\\"}END{print " \\hline \n \\end{tabular}"}' | perl -pe 's/\\Delta/\$\\Delta/' | perl -pe 's/\$  \\\\/\$  \\\\ \\hline/'
