#!/bin/sh
# $Id$
##################################################################
# This script defines a set of sample simulation cases that
# is used to verify the integrity of an IFEM-based simulator.
# Copy this file to a sub-folder Test in your App-directory and
# insert the simulations that you want to use as regression tests.
##################################################################

# Define the name of the executable here
mysim=

run () {
# This function runs a simulation with the specified options,
# pipes the terminal output to a log-file, and compares it
# with a previous simulation stored in the Reference folder.
  inp=$1
  log=`basename $1 .inp`.log
  echo Running $inp ...
  shift
  time -f "real time %E, user time %U, system time %S" \
  ../Release/bin/$mysim $inp $* > $log
  if [ ! -e Reference/$log ]; then
    mv $log Reference
  elif cmp -s $log Reference/$log; then
    echo Ok
  else
    echo "Warning: Discrepancies between current run and reference."
    diff $log Reference
  fi
}

if [ ! -d Reference ]; then mkdir Reference; fi

#####################################
# Enter the various cases below:
# Format: run <inputfile> [<options>]
#####################################

