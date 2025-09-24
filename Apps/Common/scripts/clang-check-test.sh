#!/bin/bash

# This script performs a single analysis using clang-check
# It is used by the 'make test' target in the buildsystems
# Usually you should use 'ctest -C clang-check' rather than calling this script directly
#
# Parameters: $1 = Application binary
#             $2 = Source file to process
#             $3 = Binary directory

clangcheck_cmd=$1
source_file=$2
bin_dir=$3

tmpfil=`mktemp`
$clangcheck_cmd -p $bin_dir -analyze $source_file &> $tmpfil
cat $tmpfil
if test -s $tmpfil
then
  rm $tmpfil
  exit 1
fi

rm $tmpfil
exit 0
