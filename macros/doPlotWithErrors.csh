#!/bin/csh -f

if( $#argv < 2 ) then
    echo "usage ./doPlotWithErrors.csh [ROOT_FILE] [ERA] ([MASS] [NTOYS] [DEBUG])"
    echo "[ERA] = 8TeV or 7TeV"
    exit 2
endif

set rootFile=$1

set isRootFile = `echo $rootFile | grep -c .root`

if( ! $isRootFile || ! -e $rootFile ) then
    echo "$rootFile is not a root file!"
    echo "usage ./doPlotWithErrors.csh [ROOT_FILE] [ERA] ([MASS] [NTOYS] [DEBUG])"
    exit 2
endif


set type = `echo $2 | tr "[:lower:]" "[:upper:]"`

if ($type != "8TEV" && $type != "7TEV") then
	echo "ERROR: type argument ($2) not recognized."
	echo "Choose either 8TeV or 7TeV"
	exit 2
endif

set is8TeV = `echo $type | grep -c 8`

if( $is8TeV ) then
    set dataCard_prefix = "ttH_8TeV_temp_LJ_DIL"
    set dataCard        = "ttH_8TeV_temp_LJ_DIL_ljet667_dilep3t1_dilep2t1.dat"
else
    set dataCard_prefix = "ttH_7TeV_temp_LJ_DIL"
    set dataCard        = "ttH_7TeV_temp_LJ_DIL_ljet667_dilep3t1_dilep2t1.dat"
endif

if( $#argv>2 ) then
    set mass = $3
else
    set mass = 125
endif

if( $#argv>3 ) then
    set nToys = $4
else
    set nToys = 1000
endif


if( $#argv>4 ) then
    set useDebug = `echo $5 | grep -c rue`
else
    set useDebug = 0
endif


set prefix_ttH = "ttH$mass"


set catNames = ( ljets_jge6_t2 ljets_j4_t3 ljets_j5_t3 ljets_jge6_t3 ljets_j4_t4 ljets_j5_tge4 ljets_jge6_tge4 e2je2t ge3t )

echo "Using root file $rootFile"

echo "Create datacard with all categories for nominal use"
root -b -q getCombineInfo_dilep_ljets_awesomer_extraTTBB.C'("ttH125","'$dataCard_prefix'","'$rootFile'",6,6,7,1,1,'$is8TeV')'

echo "Turn datacard into workspace"

text2workspace.py -m $mass -D data_obs $dataCard -b -o wsTest.root

echo "Run fit and save workspace"

combine -M MaxLikelihoodFit -m $mass --rMin -100 --saveWorkspace $dataCard --minos all

cp mlfit.root results_mlfit.root

echo "Making directory for plots"
if ( -e Images ) then
	set timestamp = `date +"%Y%m%d_%H%M%S"`
	set backupDir = "Images_old_$timestamp"
	echo "Moving old Images to $backupDir"
	mv Images $backupDir
endif
mkdir Images


echo "Make tables and plots"

set nCats = $#catNames
@ c = 1

while ( $c <= $nCats )
   set catName=$catNames[$c]
   echo $catName
   root -b -q head.C plotWithErrors.C+'("'${rootFile}'","'${catName}'",'$is8TeV','$nToys','$useDebug')'
   @ c += 1
end


echo "Clean up area"

rm roostats-*
rm higgsCombineTest*

echo "Done"

