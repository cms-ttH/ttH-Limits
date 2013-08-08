#!/bin/csh -f

if( $#argv < 3 ) then
    echo "usage ./createCombRootFile.csh [LJ_ROOT_FILE] [DIL_ROOT_FILE] [ERA] ([MASS] [TTBB_UNC])"
    echo "[ERA] = 8TeV or 7TeV"
    exit 2
endif

set lj_file=$1
set dil_file=$2

set isLJroot  = `echo $lj_file  | grep -c .root`
set isDILroot = `echo $dil_file | grep -c .root`

if( ! $isLJroot || ! -e $lj_file ) then
    echo "$lj_file is not a root file!"
    echo "usage ./createCombRootFile.csh [LJ_ROOT_FILE] [DIL_ROOT_FILE] [ERA] ([MASS] [TTBB_UNC])"
    exit 2
endif

if( ! $isDILroot || ! -e $dil_file ) then
    echo "$dil_file is not a root file!"
    echo "usage ./createCombRootFile.csh [LJ_ROOT_FILE] [DIL_ROOT_FILE] [ERA] ([MASS] [TTBB_UNC])"
    exit 2
endif


set type = `echo $3 | tr "[:lower:]" "[:upper:]"`

if ($type != "8TEV" && $type != "7TEV") then
	echo "ERROR: type argument ($3) not recognized."
	echo "Choose either 8TeV or 7TeV"
	exit 2
endif



if( $#argv>3 ) then
    set mass = $4
else
    set mass = 125
endif

set prefix_ttH = "ttH$mass"


if( $#argv>4 ) then
    set ttbb_unc = $5
else
    set ttbb_unc = "1.5"
endif

set is8TeV = `echo $type | grep -c 8`



echo "Hadd together LJ and DIL inputs"
hadd -f combined_lj_dil.root $lj_file $dil_file

echo "Split btag into rate + shape"
root -b -q separate_btag_rate_shape.C'("combined_lj_dil.root","combined_lj_dil_split1.root")'

echo "Split btag shape into shape by category"
root -b -q splitUncertainties_byCat.C'("combined_lj_dil_split1.root","combined_lj_dil_split2.root")'

cp combined_lj_dil_split2.root combined_lj_dil_split3.root

echo "Split Q2 scale systematic by number of partons"
root -b -q copyQ2.C

cp combined_lj_dil_split3.root combined_lj_dil_split3_binMerge.root
if( ! $is8TeV ) then
    echo "Combine edge bins for LJ 8 TeV"
    rm combined_lj_dil_split3_binMerge.root
    root -b -q combine_bins_ANN_8TeV.C'("combined_lj_dil_split3.root","combined_lj_dil_split3_binMerge.root")'
endif

echo "Add bin by bin statistical uncertainty histograms"
root -b -q statUncertainties.C'("combined_lj_dil_split3_binMerge.root","combined_lj_dil_split3_binMerge_withStats.root",'$is8TeV')'


if( $is8TeV ) then
    cp combined_lj_dil_split3_binMerge_withStats.root ttH_8TeV_temp.root
    set rootFile = "ttH_8TeV_temp.root"
    set dataCard_prefix = "ttH_8TeV_temp_LJ_DIL"
    set dataCard        = "ttH_8TeV_temp_LJ_DIL_ljet667_dilep3t1_dilep2t1.dat"
else
    cp combined_lj_dil_split3_binMerge_withStats.root ttH_7TeV_temp.root
    set rootFile = "ttH_7TeV_temp.root"
    set dataCard_prefix = "ttH_7TeV_temp_LJ_DIL"
    set dataCard        = "ttH_7TeV_temp_LJ_DIL_ljet667_dilep3t1_dilep2t1.dat"
endif

#echo "Create datacard with all categories for nominal use"
#root -b -q getCombineInfo_dilep_ljets_awesomer_extraTTBB.C'("ttH125","'$dataCard_prefix'","'$rootFile'",6,6,7,1,7,3,3,'$is8TeV',"'$ttbb_unc'")'

#echo "Remove unwanted systematics"
#perl -pi -e 's/^Q2scale_/\#Q2scale_/g' $dataCard

#echo "Do Maxlikelihood fit"
#combine -M MaxLikelihoodFit -m $mass --rMin -100 --saveWorkspace $dataCard --minos all

#echo "Turn datacard into workspace"
#text2workspace.py -m $mass -D data_obs $dataCard -b -o wsTest.root

#echo "Do limit calculation"
#combine -M Asymptotic --minosAlgo stepping -m 125 $dataCard

#echo "Make yield table with uncertainties (pre-fit)"
#root -b -q printNorms.C'("'${rootFile}'","'${prefix_ttH}'")'


echo "Cleanup stuff"
rm -f higgsCombineTest*
rm -f roostats-*
rm -f combined_lj_dil*
rm -f mlfit.root

echo Done

