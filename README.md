# Limit Setting Tools

## Setup

In your CMSSW (release 6.1.1 or higher) area, use

    cd $CMSSW_BASE/src
    addpkg HiggsAnalysis/CombinedLimit V03-01-08
    mkdir -p ttH
    git clone https://github.com/cms-ttH/ttH-Limits.git ttH/Limits
    scram b -j32

## Useful Links

* [CMS TWiki on datacard generation](https://twiki.cern.ch/twiki/bin/viewauth/CMS/NovaDatacardMaker)
* [CMS TWiki on summer 2013 limits](https://twiki.cern.ch/twiki/bin/viewauth/CMS/NovaDatacardMakerLimitRef)
