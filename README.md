# Limit Setting Tools

## Setup

In your CMSSW (release 6.1.1 or higher) area, use

    cd $CMSSW_BASE/src
    mkdir -p ttH
    git clone https://github.com/cms-ttH/ttH-Limits.git -b 94X ttH/Limits
    scram b -j32
