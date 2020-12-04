#!/bin/bash

# Hacks due to not following standards
declare -A MODULE_EXTRA_DIR
MODULE_EXTRA_DIR[IFEM-BeamEx]=IFEM-Elasticity/
MODULE_EXTRA_DIR[IFEM-FiniteDeformation]=IFEM-Elasticity/

declare -A MODULE_APP_DIR
MODULE_APP_DIR[IFEM-BeamEx]=BeamSim
MODULE_APP_DIR[IFEM-Elasticity]=Linear
MODULE_APP_DIR[IFEM-FiniteDeformation]=Nonlinear

declare -A MODULE_EXTRA_APP_DIR
MODULE_EXTRA_APP_DIR[IFEM-Elasticity]=Shell

# Parse revisions from trigger comment and setup arrays
function parseRevisions {
  for upstream in ${upstreams[*]}
  do
    if grep -qi "$upstream=" <<< $ghprbCommentBody
    then
      upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*${upstream,,}=([0-9]+).*/\1/g"`/merge
    fi
  done
  if grep -q "with downstreams" <<< $ghprbCommentBody
  then
    for sidestream in ${sidestreams[*]}
    do
      if grep -qi "$sidestream=" <<< $ghprbCommentBody
      then
        sidestreamRev[$sidestream]=pull/`echo $ghprbCommentBody | sed -r "s/.*${sidestream,,}=([0-9]+).*/\1/g"`/merge
      fi
    done
    for downstream in ${downstreams[*]}
    do
      if grep -qi "$downstream=" <<< $ghprbCommentBody
      then
        downstreamRev[$downstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*${downstream,,}=([0-9]+).*/\1/g"`/merge
      fi
    done
  fi

  # Default to a serial build if no types are given
  if test -z "$BTYPES"
  then
    BTYPES="serial"
  fi

  # Convert to arrays for easy looping
  BTYPES_ARRAY=($BTYPES)
  TOOLCHAINS=($CMAKE_TOOLCHAIN_FILES)
}


# Print revisions and configurations
function printHeader {
  echo "Repository revisions:"
  printf "\t [main library] %25s = %s\n" IFEM $IFEM_REVISION
  for upstream in ${upstreams[*]}
  do
    printf "\t     [upstream] %25s = %s\n" $upstream ${upstreamRev[$upstream]}
  done
  if [ "$1" != "IFEM" ]
  then
    printf "\t  [main module] %25s = %s\n" $1 $sha1
  fi
  if grep -q "with downstreams" <<< $ghprbCommentBody
  then
    for sidestream in ${sidestreams[*]}
    do
      printf "\t   [sidestream] %25s = %s\n" $sidestream ${sidestreamRev[$sidestream]}
    done
    for downstream in ${downstreams[*]}
    do
      printf "\t   [downstream] %25s = %s\n" $downstream ${downstreamRev[$downstream]}
    done
  fi

  echo "Configurations to process:"
  if test -n "$BTYPES"
  then
    for conf in ${!BTYPES_ARRAY[@]}
    do
      if test -n "${TOOLCHAINS[$conf]}"
      then
        echo -e "\t${BTYPES_ARRAY[$conf]} = ${TOOLCHAINS[$conf]}"
      fi
    done
  fi
}


# $1 = Additional cmake parameters
# $2 = 0 to build and install module, 1 to build and test module
# $3 = Source root of module to build
function build_module {
  cmake $3 -DCMAKE_BUILD_TYPE=Release $1
  test $? -eq 0 || exit 1
  # Threaded build
  nproc=`nproc`
  if test $2 -eq 1
  then
    cmake --build . --target testapps -- -j$nproc
    test $? -eq 0 || exit 2
    if test -z "$CTEST_CONFIGURATION"
    then
      ctest -T Test --test-timeout 180 --no-compress-output
    else
      ctest -C $CTEST_CONFIGURATION -T Test --test-timeout 180 --no-compress-output
    fi
    $WORKSPACE/deps/IFEM/jenkins/convert.py -x $WORKSPACE/deps/IFEM/jenkins/conv.xsl -t . > testoutput.xml
  else
    cmake --build . --target install -- -j$nproc
  fi
}

# $1 = Name of module
# $2 = git-rev to use for module
function clone_module {
  # Already cloned by another configuration
  if test -d $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1
  then
    return
  fi

  pushd .

  mkdir -p $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1
  cd $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1
  git init .

  # Hack due to mixed repo locations
  declare -A GH_USER
  GH_USER[IFEM-AdvectionDiffusion]=OPM
  GH_USER[IFEM-BeamEx]=SintefMath
  GH_USER[IFEM-Boussinesq]=SintefMath
  GH_USER[IFEM-Darcy]=OPM
  GH_USER[IFEM-Elasticity]=OPM
  GH_USER[IFEM-FiniteDeformation]=SintefMath
  GH_USER[IFEM-NavierStokes]=SintefMath
  GH_USER[IFEM-OpenFrac]=OPM
  GH_USER[IFEM-PoroElasticity]=OPM
  GH_USER[IFEM-Stokes]=SintefMath
  GH_USER[IFEM-ThermoElasticity]=OPM
  GH_USER[IFEM-FSI]=SintefMath
  GH_USER[IFEM-THM]=SintefMath
  GH_USER[IFEM-Poisson]=OPM
  GH_USER[SIMRA-PostProc]=akva2

  if test -n "$GH_CREDENTIALS"
  then
    git remote add origin https://$GH_CREDENTIALS@github.com/${GH_USER[$1]}/$1
  else
    git remote add origin https://github.com/${GH_USER[$1]}/$1
  fi
  git fetch --depth 1 origin $2:branch_to_build
  git checkout branch_to_build
  test $? -eq 0 || exit 1
  popd
}


# $1 = Module to clone
# $2 = Additional cmake parameters
# $3 = git-rev to use for module
# $4 = Build root
function clone_and_build_module {
  clone_module $1 $3
  pushd .
  mkdir -p $4/build-$1
  cd $4/build-$1
  test_build=0
  if test -n "$5"
  then
    test_build=$5
  fi

  build_module "$2" $test_build $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1/${MODULE_APP_DIR[$1]}
  test $? -eq 0 || exit 1
  if test -n "${MODULE_EXTRA_APP_DIR[$1]}"
  then
    mkdir -p build-${MODULE_EXTRA_APP_DIR[$1]}
    pushd build-${MODULE_EXTRA_APP_DIR[$1]}
    build_module "$2" $test_build $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1/${MODULE_EXTRA_APP_DIR[$1]}
    test $? -eq 0 || exit 1
    popd
  fi
  popd
}


# $1 - build type
# Uses pre-filled arrays upstreams, and associativ array upstreamRev
# which holds the revisions to use for upstreams.
function build_upstreams {
  for upstream in ${upstreams[*]}
  do
    if grep -qi "$upstream=" <<< $ghprbCommentBody
    then
      upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*${upstream,,}=([0-9]+).*/\1/g"`/merge
    fi
    echo "Building upstream $upstream=${upstreamRev[$upstream]} configuration=$1"

    # Build upstream and execute installation
    clone_and_build_module $upstream "-DCMAKE_PREFIX_PATH=$WORKSPACE/$1/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/$1/install -DCMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE" ${upstreamRev[$upstream]} $WORKSPACE/$1
    test $? -eq 0 || exit 1
  done
  test $? -eq 0 || exit 1
}


# $1 - name of the module we are called from
# Uses pre-filled arrays sidestreams, and associative array sidestreamRev
# which holds the default revisions to use for sidestreams
function clone_sidestreams {
  for sidestream in ${sidestreams[*]}
  do
    echo "Cloning sidestream $sidestream=${sidestreamRev[$sidestream]}"
    clone_module $sidestream ${sidestreamRev[$sidestream]}
  done
}


# $1 - name of the module we are called from
# Uses pre-filled arrays downstreams, and associativ array downstreamRev
# which holds the default revisions to use for downstreams
# Uses pre-filled arrays BTYPES_ARRAY and TOOLCHAINS.
function build_downstreams {
  for BTYPE in "${!BTYPES_ARRAY[@]}"
  do
    egrep_cmd="xml_grep --wrap testsuites --cond testsuite $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/build-$1/testoutput.xml"
    for downstream in ${downstreams[*]}
    do
      echo "Building downstream $downstream=${downstreamRev[$downstream]} configuration=${BTYPES_ARRAY[$BTYPE]}"
      # Build downstream and execute installation
      clone_and_build_module $downstream "-DCMAKE_PREFIX_PATH=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAINS[$BTYPE]}" ${downstreamRev[$downstream]} $WORKSPACE/${BTYPES_ARRAY[$BTYPE]} 1
      test $? -eq 0 || exit 1

      # Installation for downstream
      pushd .
      cd $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/build-$downstream
      cmake --build . --target install
      popd
      egrep_cmd="$egrep_cmd $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/build-$downstream/testoutput.xml"
    done

    $egrep_cmd > $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/testoutput.xml
    test $? -eq 0 || exit 1

    # Name testsuite according to configuration
    sed -i -e "s/classname=\"TestSuite\"/classname=\"${BTYPES_ARRAY[$BTYPE]}\"/g" $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/testoutput.xml
    test $? -eq 0 || exit 1
  done
}


# $1 - module name
# Build all upstreams and the context module
# Uses pre-filled arrays BTYPES_ARRAY and TOOLCHAINS.
# Uses variable IFEM_REVISION
function build_module_and_upstreams {
  for BTYPE in "${!BTYPES_ARRAY[@]}"
  do
    if [ "$1" != "IFEM" ]
    then
      pushd .
      mkdir -p ${BTYPES_ARRAY[$BTYPE]}/build-IFEM
      cd ${BTYPES_ARRAY[$BTYPE]}/build-IFEM
      build_module "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAINS[$BTYPE]} -DINSTALL_DOXY=0" 0 $WORKSPACE/deps/IFEM
      test $? -eq 0 || exit 1
      popd

      # remove FindIFEM.cmake - causes problems
      rm -f ${WORKSPACE}/deps/IFEM/cmake/Modules/FindIFEM.cmake
    fi

    CMAKE_TOOLCHAIN_FILE=${TOOLCHAINS[$BTYPE]}
    build_upstreams ${BTYPES_ARRAY[$BTYPE]}

    # Build the final module
    echo "Building main module $1=$sha1 configuration=${BTYPES_ARRAY[$BTYPE]}"

    # Make a clone
    if ! test -d deps/${MODULE_EXTRA_DIR[$1]}$1
    then
      mkdir -p deps/${MODULE_EXTRA_DIR[$1]}
      git clone $WORKSPACE deps/${MODULE_EXTRA_DIR[$1]}$1
    fi

    # Build
    pushd .
    mkdir -p ${BTYPES_ARRAY[$BTYPE]}/build-$1
    cd ${BTYPES_ARRAY[$BTYPE]}/build-$1
    build_module "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_PREFIX_PATH=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}" 1 $WORKSPACE/deps/${MODULE_EXTRA_DIR[$1]}$1/${MODULE_APP_DIR[$1]}
    test $? -eq 0 || exit 1
    cmake --build . --target install
    popd

    # Add testsuite names
    sed -e "s/classname=\"TestSuite\"/classname=\"${BTYPES_ARRAY[$BTYPE]}\"/g" ${WORKSPACE}/${BTYPES_ARRAY[$BTYPE]}/build-$1/testoutput.xml > ${WORKSPACE}/${BTYPES_ARRAY[$BTYPE]}/testoutput.xml
  done
}
