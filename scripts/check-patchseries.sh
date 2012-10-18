#!/bin/bash

# check-patchseries.sh
# Arne Morten Kvarving / SINTEF
# Oct 2012

# Convenience script that compiles and runs regression tests 
# for a series of commits. Use the CHECK_BASE enviroment variable
# to control the base of the series. It then checks the revisions
# between the base and HEAD of the current branch

# Params: See doregtests.sh
#
# Assumes it it is run from top of the source tree
# and that the script sits in the 'scripts' directory
# directly off the top of the tree.

branch=`git branch|grep \*|awk -F ' ' '{print $2}'`
base=$CHECK_BASE
head=$CHECK_HEAD
if [ -z $base ]
then
  git branch --list --remote|grep trunk && base=trunk
  git branch --list --remote|grep git-svn && base=git-svn
fi
test -z $base && base=trunk
test -z $head && last=HEAD
for rev in `git rev-list $base..$last --reverse`; do
  git checkout $rev
  if ! scripts/doregtests.sh $@
  then
    echo "Something failed for rev $rev"
    git checkout $branch
    exit 1
  fi
done
git checkout $branch
echo "Everything checks out."
