#!/bin/sh

# Light/dumb testing framework for changes to the DatacardMaker
#
# To add testcases, append a line like the following to the end of this
# script:
#
#   test_files <prefix> <root files>
#
# where <prefix> is a shorthand for the channel/combination and <root
# files> are all files needed from the $basedir defined below.
#
# The latter files get hadd'ed in the current directory, a datacard is
# built and differences against reference cards and stdout/stderr from
# programs involved are generated.
#
# Differences in output from the limit chain can be found in the *.diff
# files.
#
# When adding a new channel, the following should be executed to produce
# the reference files (posix sh):
#
#   prefix=MY_FANCY_PREFIX; for i in tmp.*; do mv $i $prefix${i#tmp}; done
#
# after adding the `test_files` line and executing this script once, with
# only the new channel/combination active.
#
# This should also be done when updating the DatacardMaker/channels.

basedir=$CMSSW_BASE/src/BEAN/limitDataCards/Summer2013

test_files() {
   prefix=$1
   shift
   for i in $*; do
      cp $basedir/$i $i
   done
   hadd -f tmp.root $* 2>&1 > tmp.hadd.out
   rm $*

   echo ==== Processing $prefix ====

   mk_datacard tmp.root > tmp.card 2> tmp.card.err
   combine -M Asymptotic --minosAlgo stepping -m 125 -t -1 tmp.card > tmp.limit.out 2> tmp.limit.err

   (diff -u {$prefix,tmp}.card 2>&1 > $prefix.card.diff && echo Card OK) || echo Card changed
   (diff -u {$prefix,tmp}.card.err 2>&1 > $prefix.card.err.diff && echo DatacardMaker stderr OK) || echo DatacardMaker stderr changed
   (diff -u {$prefix,tmp}.limit.out 2>&1 > $prefix.limit.out.diff && echo Combine stdout OK) || echo Combine stdout changed
   (diff -u {$prefix,tmp}.limit.err 2>&1 > $prefix.limit.err.diff && echo Combine stderr OK) || echo Combine stderr changed
}

test_files ditau ditau.root
