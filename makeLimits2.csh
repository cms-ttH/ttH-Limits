#! /bin/csh -f

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
	echo "[TYPE] = EXP or OBS"
	exit 2
endif

set dataCard = $1
#Get fancy and ignore case on "type"
set type = `echo $2 | tr "[:lower:]" "[:upper:]"`

if ($type != "OBS" && $type != "EXP") then
	echo "ERROR: type argument ($2) not recognized."
	echo "Choose either EXP or OBS"
	exit 2
endif

echo "Masses to consider:"
echo $masses
echo "If you want to change the masses, edit the script!"

set output = limits_${dataCard:r:t}.dat

echo "Limits will be written to $output"


echo "# Limits from $dataCard" > $output


foreach mass ($masses)

	set log = limit_${dataCard:r:t}_${mass}.log

	echo Extracting limit...
	if ($type == "EXP") then

		echo "  ==> combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 $dataCard >&! $log"
		combine -M Asymptotic --minosAlgo stepping -m $mass -t -1 $dataCard >&! $log

	else 

		echo "  ==> combine -M Asymptotic --minosAlgo stepping -m $mass $dataCard >&! $log"
		combine -M Asymptotic --minosAlgo stepping -m $mass $dataCard >&! $log

	endif

	echo "#---" >> $output
	echo $mass >> $output
	if ($type == "EXP") then
		echo "NO OBS" >> $output
	else
		grep ' r < ' $log | awk '{print $NF}' | head -1 >> $output
	endif
	
	grep ' r < ' $log | awk '{print $NF}' | tail -5 >> $output


	echo "---"

end
