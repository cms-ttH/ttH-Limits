#! /bin/tcsh -ef

set masses = (\
110\
115\
120\
125\
130\
135\
140\
)

if ($#argv != 2) then
	echo "usage ./makeLimits2.csh [DATACARD] [TYPE]"
	echo "[TYPE] = EXP or OBS or INJ"
	exit 2
endif

set dataCard = $1
#Get fancy and ignore case on "type"
set type = `echo $2 | tr "[:lower:]" "[:upper:]"`

if ($type != "OBS" && $type != "EXP" && $type != "INJ") then
	echo "ERROR: type argument ($2) not recognized."
	echo "Choose either EXP or OBS or INJ"
	exit 2
endif

echo "Masses to consider:"
echo $masses
echo "If you want to change the masses, edit the script!"

if ($type == "INJ") then
   set output = limits_${dataCard:r:t}_inj.dat
else
   set output = limits_${dataCard:r:t}.dat
endif

echo "Limits will be written to $output"


echo "# Limits from $dataCard" > $output

if ($type == "INJ") then
   set log = limit_${dataCard:r:t}_inj_prep.log
   set mass = 125
   set ml = mlfit_${dataCard:r:t}_${mass}.root
   set ws = ws_${dataCard:r:t}_${mass}.root
   set asimov = asimov_${dataCard:r:t}_${mass}.root

   echo "  ==> text2workspace.py -m $mass -D data_obs $dataCard -b -o $ws >&! $log"
   text2workspace.py -m $mass -D data_obs $dataCard -b -o $ws >&! $log
   echo "  ==> combine -M MaxLikelihoodFit -m $mass --rMin -100 --rMax 100 --minos all $ws --name _${dataCard:r:t}_${mass} >>&! $log"
   combine -M MaxLikelihoodFit -m $mass --rMin -100 --rMax 100 --minos all $ws --name _${dataCard:r:t}_${mass} >>&! $log
   echo "  ==> root -n -b -q head.C 'makeAsimovWithFitFromData.C+("'$ws'","postFitS","'$ml'","'$ws'",1.0)' >>&! $log"
   root -n -b -q head.C 'makeAsimovWithFitFromData.C+("'$asimov'","postFitS","'$ml'","'$ws'",1.0)' >>&! $log
endif

foreach mass ($masses)
   set log = limit_${dataCard:r:t}_${mass}.log

   # To speed up the following computations, append an "&" to every
   # `combine` line.  USE AT YOUR OWN RISK, AND CONSIDER YOURSELF WARNED.

   if ($type == "EXP") then
      echo "  ==> combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 $dataCard >&! $log"
      combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 $dataCard >&! $log
   else if ($type == "INJ") then
      set log = limit_${dataCard:r:t}_inj_${mass}.log

      echo "  ==> combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 --toysFile $asimov $dataCard >&! $log"
      combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 --toysFile $asimov $dataCard >&! $log
   else
      echo "  ==> combine -M Asymptotic --minosAlgo stepping -m $mass $dataCard >&! $log"
      combine -M Asymptotic --minosAlgo stepping -m $mass $dataCard >&! $log
   endif
end

wait

foreach mass ($masses)
   if ($type == "INJ") then
      set log = limit_${dataCard:r:t}_inj_${mass}.log
      echo $mass >> $output
      grep ' r < ' $log | awk '{print $NF}' | head -1 >> $output
   else
      set log = limit_${dataCard:r:t}_${mass}.log
      echo $mass >> $output
      if ($type == "EXP") then
         echo "NO OBS" >> $output
      else
         grep ' r < ' $log | awk '{print $NF}' | head -1 >> $output
      endif

      grep ' r < ' $log | awk '{print $NF}' | tail -5 >> $output
   endif
end
