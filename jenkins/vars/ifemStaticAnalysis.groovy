// Builds an IFEM module in the static analysis configurations, one stage per
// configuration built one after another to limit the load, and analyses it:
// cppcheck, clang-tidy and infer on the build of the infer configuration,
// and the test coverage of the coverage_gcov configuration with gcovr.
// The memcheck configuration, running the tests through valgrind, is built
// with all downstreams, unless disabled for modules which do not need it.
// The results are recorded with the Warnings Next Generation and Coverage
// plugins.
//
// Arguments (all optional):
//   module         - name of the module (default: the repository name)
//   configurations - script on the builder defining the static analysis
//                    configurations, see jenkins/README.md
//                    (default '/home/jenkins/build-configurations-sca.sh')
//   suppressions   - cppcheck suppressions list
//                    (default '/var/lib/jenkins/cppcheck-suppressions/<module>')
//   label          - agent label to run the configurations on
//   memcheck       - build the memcheck configuration if the configurations
//                    have one (default true)
def call(Map args = [:]) {
  String module = args.module ?: repositoryName()
  String suppressions = args.suppressions ?: "/var/lib/jenkins/cppcheck-suppressions/${module}"

  ifemModule(configurations: args.configurations ?: '/home/jenkins/build-configurations-sca.sh',
             label: args.label,
             statusContexts: false,
             sequential: true,
             exclude: args.containsKey('memcheck') && !args.memcheck ? ['memcheck'] : [],
             configEnv: [memcheck: ['ghprbCommentBody=with downstreams',
                                    'OPENBLAS_NUM_THREADS=1']],
             afterBuild: [infer: { analyse(module, suppressions) },
                          coverage_gcov: { coverage(module) }])
}


// Runs cppcheck, clang-tidy and infer on the compilation database of the
// infer build. A failing tool marks the build unstable, the others still run.
def analyse(String module, String suppressions) {
  String build = "infer/build-${module}"
  catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
    sh "cppcheck --suppressions-list=${suppressions} --project=\$WORKSPACE/${build}/compile_commands.json " +
       '-q --force --enable=all --xml --xml-version=2 -j16 --output-file=$WORKSPACE/infer/cppcheck-result.xml'
  }
  catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
    sh "#!/bin/bash\nset -o pipefail\n" +
       "run-clang-tidy -checks='clang-analyzer.*' -p \$WORKSPACE/${build}/ -j16 '.*C\$' | tee clang-tidy-report.log"
  }
  catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
    dir(build) {
      sh 'infer run --compilation-database compile_commands.json --keep-going --pmd-xml'
    }
  }
  recordIssues id: 'static-analysis', name: 'Static analysis',
               aggregatingResults: false, enabledForFailure: true,
               skipPublishingChecks: true,
               tools: [cppCheck(pattern: 'infer/cppcheck-result.xml'),
                       clangTidy(pattern: 'clang-tidy-report.log'),
                       pmdParser(id: 'infer', name: 'Infer', pattern: "${build}/infer-out/report.xml")]
}


// Collects the test coverage of the coverage_gcov build. Functions defined
// differently in different translation units, e.g. depending on BLAS,
// are kept as separate functions.
def coverage(String module) {
  catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
    dir("coverage_gcov/build-${module}") {
      sh "gcovr --xml -r \$WORKSPACE -o \$WORKSPACE/coverage.xml -e \$WORKSPACE/coverage_gcov/install " +
         "-e \$WORKSPACE/coverage_gcov/build-${module}/CMakeFiles -e deps/IFEM " +
         '--merge-mode-functions=separate'
    }
  }
  recordCoverage id: 'coverage', name: 'Coverage',
                 tools: [[parser: 'COBERTURA', pattern: 'coverage.xml']],
                 skipPublishingChecks: true
}


// The name of the repository the job builds
String repositoryName() {
  String url = scm.userRemoteConfigs[0].url
  return url.tokenize('/:')[-1] - ~/\.git$/
}
