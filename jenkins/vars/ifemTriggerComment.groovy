// Returns the body of the pull request comment that triggered this build,
// or an empty string if the build was not triggered by a comment.
// Requires the GitHub PR Comment Build plugin (github-pr-comment-build).
String call() {
  def causes = currentBuild.getBuildCauses('com.adobe.jenkins.github_pr_comment_build.GitHubPullRequestCommentCause')
  return causes ? causes[-1].commentBody : ''
}
