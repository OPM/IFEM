#!/bin/bash

# Convenience script that compiles and runs all regression tests

ROOT=`readlink -f $0`
ROOT=`dirname $ROOT`/..
PATHS=". Apps/Stokes Apps/FiniteDefElasticity"

rm -f $ROOT/failed.log

globres=0
for p in $PATHS
do
  cd $ROOT/$p
  mkdir Release-Testing
  cd Release-Testing
  cmake $ROOT/$p -DDISABLE_HDF5=1 -DCMAKE_BUILD_TYPE=Release -DIFEM_BUILD_TYPE=Release-Testing
  make -j3
  if ! make test
  then
    globres=1
  fi
done

if test $globres -eq 1
then
  echo "Some tests failed, see failed.log in root of repository for details"
fi

for p in $PATHS
do
  rm $ROOT/$p/Release-Testing -rf
done
