#!/bin/bash

# Downstream revisions
declare -a downstreams
downstreams=(IFEM-AdvectionDiffusion
             IFEM-Darcy
             IFEM-Poisson
             IFEM-Stokes
             IFEM-NavierStokes
             IFEM-Boussinesq
             IFEM-Elasticity
             IFEM-BeamEx
             IFEM-FiniteDeformation
             IFEM-ThermoElasticity
             IFEM-PoroElasticity
             IFEM-THM
             IFEM-OpenFrac
             IFEM-FSI
             SIMRA-PostProc)

declare -A downstreamRev
downstreamRev[IFEM-Poisson]=master
downstreamRev[IFEM-AdvectionDiffusion]=master
downstreamRev[IFEM-BeamEx]=master
downstreamRev[IFEM-Darcy]=master
downstreamRev[IFEM-Elasticity]=master
downstreamRev[IFEM-FiniteDeformation]=master
downstreamRev[IFEM-NavierStokes]=master
downstreamRev[IFEM-Boussinesq]=master
downstreamRev[IFEM-OpenFrac]=master
downstreamRev[IFEM-PoroElasticity]=master
downstreamRev[IFEM-THM]=master
downstreamRev[IFEM-Stokes]=master
downstreamRev[IFEM-ThermoElasticity]=master
downstreamRev[IFEM-FSI]=master
downstreamRev[SIMRA-PostProc]=master

IFEM_REVISION=$sha1

source `dirname $0`/build-ifem-module.sh

parseRevisions
printHeader IFEM

build_module_and_upstreams IFEM

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  exit 0
fi

# remove cmake rule so apps do not get confused
mv $WORKSPACE/cmake/Modules/FindIFEM.cmake $WORKSPACE

build_downstreams IFEM

# move cmake rule back in place
mv $WORKSPACE/FindIFEM.cmake $WORKSPACE/cmake/Modules

test $? -eq 0 || exit 1
