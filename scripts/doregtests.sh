#!/bin/bash

# doregtest.sh
# Arne Morten Kvarving / SINTEF
# Oct 2012

# Convenience script that compiles and runs regression tests
# Params: $1 - build configuration
# Params: $2..$@ - directories to run tests in
#
# Assumes it it is run from top of the source tree
# and that the script sits in the 'scripts' directory
# directly off the top of the tree.

# Configure and build a directory
# Params: $1 - configuration directory
#         $2 - root configuration directory
#         $3 - directory to configure for
Build() {
  if ! test -d $1
  then
    # Copy config off root dir
    config=`cmake -L $2|grep =`
    cmakeopt=""
    IFS=$'\n'
    for line in $config; do
      opt=`echo $line | sed -e 's/=/="/g' -e 's/$/"/g'`
      cmakeopt="-D$opt $cmakeopt"
    done
    mkdir $1
    cd $1
    cmake $3 $cmakeopt
  else
    # Make sure build system is up to date
    cd $1
    cmake $3
  fi
  threads=$CHECK_THREADS
  test -z $thread && threads=3
  if ! make -j$threads
  then
    globres=1
    return
  fi
}

globres=0

ROOT=`readlink -f $0`
ROOT=`dirname $ROOT`/..

rm -f $ROOT/failed.log

# Build root lib
rootconfig=$ROOT/$1
test -d $rootconfig || rootconfig=$ROOT
cd $rootconfig
Build $rootconfig $rootconfig $ROOT

if test $globres -eq 0
then
  for param in ${@:2}; do
    test -d $ROOT/$1 && builddir=$ROOT/$param/$1
    test -d $ROOT/$1 || builddir=$ROOT/$param
    if [ "$param" != "." ]
    then
      Build $builddir $rootconfig $ROOT/$param
    fi
    if test $globres -eq 0
    then
      if ! make test
      then
        globres=1
      fi
    fi
  done
fi

if test $globres -eq 1
then
  echo "Some tests failed, see failed.log in root of repository for details"
fi

exit $globres
