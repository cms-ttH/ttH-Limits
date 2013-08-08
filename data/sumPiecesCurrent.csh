#!/bin/tcsh -f 

set LJ = RENAMED_ann_LJ_8TeV_2012ABCD.root
set OS = DIL_OS_BASELINE.root
#set SS = histosForLimits_AllLep_2012_53x_SS_AllTag.root
#set SS = histosForLimits_SS_AllLep_AllTag_v2.root
set SS = SS_newANN_newNPSF_AllLep_2012_53x_SS_AllTag.root

set DILCOMBINED = dilCombined_OS_SS.root


# add up OS and SS

hadd -f $DILCOMBINED $OS $SS



# ---- create ----
echo "Putting a job in the background with output to a log file!" 
./createCombRootFile.csh $LJ $DILCOMBINED 8TeV > & ! crashes.log &


