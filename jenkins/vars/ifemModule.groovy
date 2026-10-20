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
//   credentialsId  - username/password or GitHub App credentials used to
//                    clone private upstream and downstream repositories
//   statusContexts - report a commit status per configuration
//                    (default true, requires pipeline-githubnotify-step)
def call(Map args = [:]) {
  Map opts = [script: args.script ?: 'jenkins/build.sh',
              configurations: args.configurations ?: '/home/jenkins/build-configurations.sh',
              credentialsId: args.credentialsId,
              statusContexts: args.containsKey('statusContexts') ? args.statusContexts : true]
  String label = args.label

  String comment = ifemTriggerComment()
  if (comment)
    echo "Triggered by comment: ${comment}"

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

  parallel branches
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

    // Select the toolchain belonging to this configuration
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
    ${opts.credentialsId ? 'export GH_CREDENTIALS=x-access-token:$GH_TOKEN' : ''}
    exec ${opts.script}
  fi
done
echo "Configuration \$BTYPE not found in ${opts.configurations}"
exit 1
"""

    withEnv(["sha1=${sha1}",
             "ghprbCommentBody=${comment}",
             "BTYPE=${name}"]) {
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
  }
}
