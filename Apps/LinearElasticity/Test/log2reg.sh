#!/bin/sh
# $Id$
# This script creates new reg-files from the specified log-files.
for file in $*; do if [ -f $file ]; then found=yes; echo $file
sed -f ../../../test/log2reg.sed $file | sed \
'1,3 d;4 s/^.*LinEl //;4 s/ *-vtf *1//;4 s/ *-nviz *[1-9]//;4 s/ *-n[uvw] *[1-9]//g;4 s/ *-hdf5//; 4 a\
' > `basename $file .log`.reg
fi; done
if [ -z "$found" ]; then echo "usage $0 <logfile(s)>"; fi
