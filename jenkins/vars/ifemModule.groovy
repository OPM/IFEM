// Builds an IFEM module in a multibranch pipeline, one parallel stage per
// build configuration. Each stage runs the module's jenkins build script
// with the environment the old pull request builder jobs provided.
//
// The build configurations are machine specific and are read from a script
// on the builder which sets BTYPES and CMAKE_TOOLCHAIN_FILES, see
// jenkins/README.md.
//
// Arguments (all optional):
//   script         - build script, relative to the repository root
//                    (default 'jenkins/build.sh')
//   configurations - script on the builder defining the build configurations
//                    (default '/home/jenkins/build-configurations.sh')
//   configs        - list of configuration names to build
//                    (default: all configurations in BTYPES)
//   label          - agent label to run the configurations on
//   credentialsId  - username/password credentials used to clone private
//                    upstream and downstream repositories
//                    (default: the credentials of the job's GitHub source)
//   statusContexts - report a commit status per configuration
//                    (default true, requires pipeline-githubnotify-step)
//   branchTriggers - map of branch name to the trigger comment to use for
//                    builds of that branch, e.g. post-merge builds of master
//   keepWorkspace  - keep the workspaces after the build, e.g. for debugging
//                    (default false)
//   afterBuild     - map of configuration name to a closure run in its
//                    workspace after a successful build, e.g. for analysis
//   sequential     - build the configurations one after another instead of
//                    in parallel, e.g. to limit the load (default false)
//   configEnv      - map of configuration name to a list of 'NAME=value'
//                    environment variables for its build, which may also
//                    override ghprbCommentBody, e.g. to build downstreams
def call(Map args = [:]) {
  Map opts = [script: args.script ?: 'jenkins/build.sh',
              configurations: args.configurations ?: '/home/jenkins/build-configurations.sh',
              credentialsId: args.credentialsId ?: scmCredentialsId(),
              statusContexts: args.containsKey('statusContexts') ? args.statusContexts : true,
              keepWorkspace: args.keepWorkspace ?: false,
              afterBuild: args.afterBuild ?: [:],
              configEnv: args.configEnv ?: [:]]
  String label = args.label

  String comment = ifemTriggerComment()
  if (comment)
    echo "Triggered by comment: ${comment}"
  else if (!env.CHANGE_ID && args.branchTriggers?.get(env.BRANCH_NAME)) {
    comment = args.branchTriggers[env.BRANCH_NAME]
    echo "Using trigger for branch ${env.BRANCH_NAME}: ${comment}"
  }

  // The configurations may depend on the trigger comment
  List names = args.configs
  if (!names) {
    node(label) {
      withEnv(["ghprbCommentBody=${comment}"]) {
        names = sh(returnStdout: true,
                   script: "#!/bin/bash\nsource ${opts.configurations}\necho \$BTYPES").trim().tokenize()
      }
    }
  }
  if (!names)
    error "No build configurations given in ${opts.configurations}"

  def branches = [:]
  for (int i = 0; i < names.size(); ++i) {
    String name = names[i]
    branches[name] = {
      node(label) {
        stage(name) {
          buildConfig(name, comment, opts)
        }
      }
    }
  }

  if (!args.sequential) {
    parallel branches
    return
  }

  // Build all configurations even if one fails, but stop when aborted
  def failure = null
  for (int i = 0; i < names.size(); ++i) {
    try {
      branches[names[i]]()
    } catch (org.jenkinsci.plugins.workflow.steps.FlowInterruptedException err) {
      throw err
    } catch (err) {
      failure = failure ?: err
    }
  }
  if (failure)
    throw failure
}


def buildConfig(String name, String comment, Map opts) {
  String context = "jenkins/${name}"
  if (opts.statusContexts)
    githubNotify context: context, status: 'PENDING', description: 'Building'

  boolean failed = false
  try {
    // The build scripts reuse anything already cloned into deps/
    deleteDir()
    def scmVars = checkout scm
    String sha1 = env.CHANGE_ID ? "pull/${env.CHANGE_ID}/merge" : scmVars.GIT_COMMIT
    String log = "build-${name}.log"

    // Select the toolchain belonging to this configuration.
    // The output is also kept in a log of its own for the warnings parsers,
    // as the console log mixes the parallel configurations.
    String build = """#!/bin/bash
source ${opts.configurations}
types=(\$BTYPES)
toolchains=(\$CMAKE_TOOLCHAIN_FILES)
for i in \${!types[@]}
do
  if [ "\${types[\$i]}" = "\$BTYPE" ]
  then
    export BTYPES=\$BTYPE
    export CMAKE_TOOLCHAIN_FILES=\${toolchains[\$i]}
    ${opts.credentialsId ? 'export GH_CREDENTIALS=$GH_USER:$GH_TOKEN' : ''}
    set -o pipefail
    ${opts.script} 2>&1 | tee ${log}
    exit
  fi
done
echo "Configuration \$BTYPE not found in ${opts.configurations}"
exit 1
"""

    withEnv(mergeEnv(["sha1=${sha1}",
                      "ghprbCommentBody=${comment}",
                      "BTYPE=${name}"], opts.configEnv[name] ?: [])) {
      if (opts.credentialsId) {
        withCredentials([usernamePassword(credentialsId: opts.credentialsId,
                                          usernameVariable: 'GH_USER',
                                          passwordVariable: 'GH_TOKEN')]) {
          sh build
        }
      } else {
        sh build
      }
    }

    if (opts.afterBuild[name])
      opts.afterBuild[name]()
  } catch (err) {
    failed = true
    if (opts.statusContexts)
      githubNotify context: context, status: 'FAILURE', description: 'Build failed'
    throw err
  } finally {
    // Ensure the test report is collected even if a downstream fails to build
    def results = junit testResults: "${name}/testoutput.xml", allowEmptyResults: true
    if (opts.statusContexts && !failed) {
      if (results.failCount > 0)
        githubNotify context: context, status: 'FAILURE',
                     description: "${results.failCount} of ${results.totalCount} tests failed"
      else
        githubNotify context: context, status: 'SUCCESS',
                     description: "${results.totalCount} tests passed"
    }
    recordWarnings(name)
    archiveArtifacts artifacts: "${name}/**/failed.log", allowEmptyArchive: true

    // The build trees of IFEM and all downstreams are large
    if (!opts.keepWorkspace)
      deleteDir()
  }
}


// Adds the variables of extra to env, replacing those with the same name
@NonCPS
List mergeEnv(List env, List extra) {
  List names = extra.collect { it.toString().tokenize('=')[0] }
  return env.findAll { !names.contains(it.toString().tokenize('=')[0]) } + extra
}


// Records compiler and cmake warnings from the build log of a configuration,
// and doxygen warnings from the DoxyWarnings.log in each build directory.
// Requires the Warnings Next Generation plugin.
def recordWarnings(String name) {
  String log = "build-${name}.log"
  // gcc and clang warnings look the same, use the parser of the compiler used
  boolean isClang = sh(returnStatus: true,
                       script: "grep -q 'CXX compiler identification is Clang' ${log}") == 0
  recordIssues id: "warnings-${name}", name: "Warnings (${name})",
               aggregatingResults: true, enabledForFailure: true,
               skipPublishingChecks: true,
               tools: [isClang ? clang(pattern: log) : gcc(pattern: log),
                       gnuFortran(pattern: log),
                       doxygen(pattern: "${name}/**/DoxyWarnings.log"),
                       cmake(pattern: log)]
}


// The credentials the multibranch job checks out its own repository with
String scmCredentialsId() {
  try {
    return scm.userRemoteConfigs[0].credentialsId
  } catch (err) {
    return null
  }
}
