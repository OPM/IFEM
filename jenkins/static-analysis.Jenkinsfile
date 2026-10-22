@Library('ifem') _

// Analyse master every night if it has changed
properties([pipelineTriggers([pollSCM('H H(0-5) * * *')])])

ifemStaticAnalysis()
