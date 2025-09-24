#!/bin/bash

# This script performs a single IO regression test.
# It is used by the 'make test' target in the buildsystems
# Usually you should use 'make test' rather than calling this script directly
#
# Parameters: -a = Application binary
#             -r = Regression test file
#             -t = File type to generate/test
#             -d = Binary directory
#             -p = Listing app for file type
# Optional parameters:
#             -v = Valgrind binary
#             -m = Number of MPI procs
# A IO regression test file is of the format:
# parameters to application
# blank line
# output from listing application

set -e

function atomic_log_file {
(
  flock -x 200
  cat $1 >> $BINARY_DIR/failed.log

) 200>/var/lock/ifemloglock
}

OPTIND=1
while getopts "a:r:t:m:p:v:d:" OPT
do
  case "${OPT}" in
    a) MYSIM=${OPTARG} ;;
    r) REG_FILE=${OPTARG} ;;
    t) TYPE=${OPTARG} ;;
    m) MPI_PROCS=${OPTARG} ;;
    d) BINARY_DIR=${OPTARG} ;;
    p) LISTING_APP=${OPTARG} ;;
    v) VALGRIND=${OPTARG} ;;
  esac
done


test -f $REG_FILE || exit 1

cd `dirname $REG_FILE`
MAPFILE=`head -n1 $REG_FILE`
tmplog=`mktemp -t ifemlogXXXXXX`
tmpnew=`mktemp -t ifemlogXXXXXX`
tmpold=`mktemp -t ifemlogXXXXXX`
output=`mktemp -t ifemXXXXXX.$TYPE`

if test -n "${VALGRIND}"
then
  vallog=`mktemp -t ifemvallogXXXXXX`
  MYSIM="${VALGRIND} --log-file=${vallog} ${MYSIM}"
fi

if [ "$TYPE" == "hdf5" ]
then
  test -n "$MPI_PROCS" && MYSIM="mpirun -n $MPI_PROCS $MYSIM"
  $MYSIM $MAPFILE -hdf5 $output 2>&1 | tee $tmplog
else
  $MYSIM $MAPFILE -vtf 0 -vtffile $output 2>&1 | tee $tmplog
fi
if test $? -ne 0
then
  atomic_log_file $tmplog
  exit 1
fi
tail -n +3 $REG_FILE > $tmpold
if [ "$TYPE" == "hdf5" ]
then
  $LISTING_APP -r $output > $tmpnew
else
  $LISTING_APP $output > $tmpnew
fi

diff -u $tmpold $tmpnew > $tmplog
res=$?
if test $res -ne 0
then
  faillog=`mktemp -t ifemfailXXXXXX`
  echo "Failure for $REG_FILE" > $faillog
  echo "======================" >> $faillog
  cat $tmplog >> $faillog
  atomic_log_file $faillog
  cat $faillog
  rm -f $faillog
fi
rm $tmpnew $tmpold $tmplog $output

if test -n "$VALGRIND"
then
  if ! grep -q "ERROR SUMMARY: 0 errors" $vallog
  then
    res=1
    atomic_log_file $vallog
  fi
  rm -f $vallog
fi

exit $res
