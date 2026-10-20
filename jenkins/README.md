# IFEM jenkins build scripts:

**build-ifem-module.sh**:
This is a helper script which contains functions for building,
testing and cloning modules.

**build.sh**:
This expects to run on a jenkins instance with IFEM as the 'origin' remote.

It will build and test IFEM. It can be used both for post-merge builds
of the master branch and for pull request builds.
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

# Jenkins pipeline

The `jenkins/` directory in the IFEM repository is a Jenkins shared library.
`ifemModule()` runs the module's `jenkins/build.sh` once per build
configuration, as parallel stages in a single pipeline, and reports a
commit status (`jenkins/<configuration>`) and test results per configuration.
Each repository contains a `jenkins/Jenkinsfile`:

    @Library('ifem') _

    ifemModule()

The build configurations are machine specific and live on the builder in
`/home/jenkins/build-configurations.sh`, which sets `BTYPES` and
`CMAKE_TOOLCHAIN_FILES` as described above. One stage is created per entry
in `BTYPES`, and each stage runs the build script with only its own
configuration and toolchain file.

Optional arguments, e.g. for a module with a different build script,
private dependencies or a subset of the configurations:

    ifemModule(script: 'jenkins/build.sh',
               configurations: '/home/jenkins/build-configurations.sh',
               configs: ['serial', 'clang'],
               label: 'linux',
               credentialsId: 'github-token',
               statusContexts: true,
               keepWorkspace: false)

Each configuration wipes its workspace before building and, unless
`keepWorkspace` is set, again after the test results are recorded.

Each configuration also records the warnings from gcc or clang, gfortran
and cmake in its build output, and from doxygen in the DoxyWarnings.log of
each build directory, as `Warnings (<configuration>)`. The failed.log of
each build directory is archived.

Private upstreams and downstreams are cloned with GH\_CREDENTIALS set from
a username/password credential, by default the one used by the job's GitHub
branch source.

Builds of a branch, such as post-merge builds of master, have no trigger
comment. `branchTriggers` gives the comment to use for them instead; IFEM
builds master with all downstreams:

    ifemModule(branchTriggers: [master: 'jenkins build this with downstreams please'])

## Jenkins setup

Required plugins: Pipeline, GitHub Branch Source, GitHub PR Comment Build,
Pipeline: GitHub Notify Step, JUnit and Warnings Next Generation.

1. Manage Jenkins → System → Global Trusted Pipeline Libraries: add a library
   named `ifem`, default version `master`, retrieved from the IFEM repository
   (github.com/OPM/IFEM) with Library Path `jenkins/`.
2. Create a multibranch pipeline (or an organization folder) per repository
   using the GitHub branch source, preferably with GitHub App credentials.
   Set the Script Path (Build Configuration for a multibranch pipeline,
   Project Recognizers for an organization folder) to `jenkins/Jenkinsfile`.
3. Under the branch source, add the property strategy
   "Trigger build on pull request comment" with the expression
   `jenkins build this.*please`. To only build pull requests on request,
   use "Named branches get different properties" and add
   "Suppress automatic SCM triggering" to all branches except master,
   so post-merge builds of master still happen.
4. Make sure the repository webhooks (or the GitHub App) deliver
   issue comment events to Jenkins.

The trigger syntax is unchanged, e.g.

    jenkins build this with downstreams ifem-poisson=1 ifem-elasticity=2 please
