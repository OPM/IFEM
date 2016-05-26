# opm-common jenkins build scripts:

**build-ifem-module.sh**:
This is a helper script which contains functions for building,
testing and cloning modules.

**build.sh**:
This expects to run on a jenkins instance with IFEM as the 'origin' remote.

It will build and test IFEM. It can be used both for post-merge builds
of the master branch and for a github pull request builder job.
Optionally you it can build all downstreams and execute their tests.

Important environment variables (jenkins conventions):

WORKSPACE - Root of source tree.

ghprbCommentBody - Trigger as used on github.

GH\_CREDENTIALS - Github credentials to use in the form user:password.

CMAKE\_TOOLCHAIN\_FILES - A space separated list of toolchain files specifying
                          cmake parameters for a configuration.

BTYPES - A space separated list of build configuration names.

A typically local run building the current content of IFEM and all
downstreams;

WORKSPACE=`pwd` GH\_CREDENTIALS=user:pwd ghprbCommentBody="jenkins build this with downstreams please" jenkins/build.sh
