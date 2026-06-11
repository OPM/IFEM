#!/bin/bash
# Checks out previous version of the doc from the gh-pages branch, and
# compares with the newly generated doc while ignoring generated time stamps.

if [ ! -e Release/doc/html/index.html ]; then
  echo Expected file doc/html/index.html is missing.
  echo \"make doc\" probably failed.
  exit 2
fi

git fetch origin gh-pages || exit 3
git checkout -B gh-pages origin/gh-pages || exit 3
git rm -r --ignore-unmatch docs || exit 3
mv Release/doc/html docs
git add docs
if git diff -I "^Generated on .* for IFEM" HEAD --quiet; then
  echo No changes in source code documentation
  git reset --hard HEAD
  exit 0
else
  exit 1
fi
