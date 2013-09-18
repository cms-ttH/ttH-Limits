#! /bin/tcsh -ef

set chanName = $1
set datacard = $2
set mass = $3

echo "Running jobs for $chanName ($datacard) at mass = $mass GeV"

# First, let's make a workspace.  No need to do this over and over again
set wsName = ws_ttH_${chanName}_${mass}.root
pushd ../data
echo "==>Convert card to workspace $wsName"
text2workspace.py -m $mass -D data_obs $datacard -b -o $wsName
mv $wsName ../macros/.
popd

# Now, let's get the limit job started.  It can be set to run in the background because nothing else depends on it.
set log = limit_${chanName}_${mass}.log
echo "==> combine -M Asymptotic --minosAlgo stepping -m $mass $wsName >&! $log &  (Running in background...)"
combine -M Asymptotic --minosAlgo stepping -m $mass $wsName >&! $log &

# Alright, while that's working, let's do a set of sequential things (fit -> signal injection -> expected with signal)

# Run the fit.  Don't go on after this because the next step relies on it
set log = fit_${chanName}_${mass}.log
set mlfitName = mlfit_${chanName}_${mass}.root
echo "==> combine -M MaxLikelihoodFit -m $mass --rMin -100 --rMax 100 --minos all --name _${chanName}_${mass} $wsName > & ! $log"
combine -M MaxLikelihoodFit -m $mass --rMin -100 --rMax 100 --minos all --name _${chanName}_${mass} $wsName > & ! $log

# Now that we've got the fit, let's inject some signal into the best fit.
echo "==> Create signal injected Asimov data"
set sigInjName = ttH_${chanName}_${mass}_sigInj.root
root -n -b -q head.C 'makeAsimovWithFitFromData.C+("'$sigInjName'","postFitS","'$mlfitName'","'$wsName'",1.0)'

# Finally, run the limits with signal injected into the Asimov data...
set log = limit_${chanName}_${mass}_sig.log

echo "==> combine -M Asymptotic --minosAlgo stepping -m $mass --toysFile $sigInjName -t -1 $wsName >&! $log"
combine -M Asymptotic --minosAlgo stepping -m $mass --toysFile $sigInjName -t -1 $wsName >&! $log

echo "==> Cleaning up some files"
rm -f higgsCombine*.root

echo '-----DONE!---------'
echo ''






