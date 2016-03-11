# Limit Setting Tools

## Setup

In your CMSSW (release 6.1.1 or higher) area, use

    cd $CMSSW_BASE/src
    addpkg HiggsAnalysis/CombinedLimit V03-01-08
    mkdir -p ttH
    git clone https://github.com/cms-ttH/ttH-Limits.git ttH/Limits
    # using https
    git clone https://gitlab.cern.ch/matze/ttH-Limits-data.git ttH/Limits/data
    # using kerberos
    git clone https://:@gitlab.cern.ch:8443/matze/ttH-Limits-data.git ttH/Limits/data
    scram b -j32

The latter two git commands clone our [limit data repository](https://git.cern.ch/web/?p=ttH-Limits-data.git;a=summary),
which has restricted access.
These steps may be changed to use a submodule at a later point.

## Notes

### Reducing H → gg filesizes

After running `scram`, run `reduce_gamma_filesize in.root out.root` to
reduce the filesize of the gamma files to an acceptable level.

## Useful Links

* [CMS TWiki on datacard generation](https://twiki.cern.ch/twiki/bin/viewauth/CMS/NovaDatacardMaker)
* [CMS TWiki on summer 2013 limits](https://twiki.cern.ch/twiki/bin/viewauth/CMS/NovaDatacardMakerLimitRef)
