#! /bin/tcsh -ef

set mass = 125.7

set channels = (\
photons:ttH_Hgg.txt \
hadrons:ttH_hbb_htt_7TeV_8TeV.txt \
leptons:combBCat_QMVA.card.txt \
combination:ttH_all.txt \
)

foreach chan ($channels)

	set chanName = `echo $chan | awk -F ':' '{print $1}'`
	set datacard = `echo $chan | awk -F ':' '{print $2}'`

	set log = job_${chanName}_${mass}.log

    echo "=> Launching $chanName ($datacard) for mass $mass GeV"
	./runAllChan.csh $chanName $datacard $mass >&! $log &


end

