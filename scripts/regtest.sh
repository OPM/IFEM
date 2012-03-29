#!/bin/bash

# This script performs a single regression.
# It is used by the 'make test' target in the buildsystems
# Usually you should use 'make test' rather than calling this script directly
#
# Parameters: $1 = Application binary
#             $2 = Regression test file
#
# A regression test file is of the format:
# parameters to application
# blank line
# output from program to compare to

mysim=$1

cd `dirname $2`
readarray < $2
$mysim $MAPFILE 2>&1 > templog
globres=1
IFS=$'\n'
for line in `cat $2`
do
  test -z "$line" && continue
  echo "$line" | grep -q ".inp" && continue
  result=0
  if grep -q "$line" templog 
    then result=1 
  fi
  if test $result -eq 0 
  then 
    if test $globres -eq 1
    then
      echo "-------- $2 --------" >> ../../../failed.log
    fi
    globres=0
    echo Failed to find output: $line >> ../../../failed.log
  fi
done

if test $globres -eq 0
then
  cat templog >> ../../../failed.log
  rm templog
  exit 1
fi

rm templog
exit 0
