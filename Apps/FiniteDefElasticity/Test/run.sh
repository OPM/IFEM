#!/bin/sh
# $Id$
##################################################################
# This script defines a set of sample simulation cases that
# is used to verify the integrity of the IFEM-based simulator.
# Copy this file to a sub-folder Test in your App-directory and
# insert the simulations that you want to use as regression tests.
##################################################################

# Define the name of the executable here
mysim=NonLinEl

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

# 3D cantilever beam, linear-elastic material
for input in CanTS-p?.inp; do
  ln -s $input TL-$input
  run TL-$input -dense -vtf 1
  ln -s $input UL-$input
  run UL-$input -UL -dense -vtf 1
  rm ?L-$input
done

# 3D Thick cylinder compression, linear-elastic material
run Cyl-p2.inp -UL -vtf 1 -nGauss 3 -nviz 3
run Cyl-p3.inp -UL -vtf 1 -nGauss 4 -nviz 3
run Cyl-p4.inp -UL -vtf 1 -nGauss 4 -nviz 3

# 2D Rubber block, hyperelastic Neo-Hooke material
# - mixed formulation with internal pressure modes
run FBlock-h9x5-Q2P1.inp -2D -MX 1 -nGauss 3 -vtf 1 -lagrange
run FBlock-h8x3-Q3P2.inp -2D -MX 2 -nGauss 4 -vtf 1 -lagrange
run FBlock-h8x2-Q4P3.inp -2D -MX 3 -nGauss 5 -vtf 1 -lagrange
# - mixed formulation with continuous pressure field
run FBlock-h9x5-Q2Q1.inp -2D -mixed -nGauss 3 -vtf 1 -lagrange
run FBlock-h8x3-Q3Q2.inp -2D -mixed -nGauss 4 -vtf 1 -lagrange
run FBlock-h8x2-Q4Q3.inp -2D -mixed -nGauss 5 -vtf 1 -lagrange

# Tension of a 2D elasto-plastic strip
run Necking-Q2P1.inp -2D -MX 1 -vtf 1 -nviz 3 -nGauss 3 -outPrec 6
run Necking-Q2Q1.inp -2D -mixed -vtf 1 -nviz 3 -nGauss 3 -outPrec 6
run Necking-Q2-Q1.inp -2D -Mixed -vtf 1 -nviz 3 -nGauss 3 -outPrec 6

# Cooks membrane, 2D plasticity
run Cook2D-p1-h1.inp -2D -Fbar 1 -nGauss 2 -vtf 1 -nviz 3 -outPrec 6
run Cook2D-p2-h1.inp -2D -Fbar 2 -nGauss 3 -vtf 1 -nviz 3 -outPrec 6
run Cook2D-h1-p2.inp -2D -Fbar 2 -nGauss 3 -vtf 1 -lagrange -outPrec 6
