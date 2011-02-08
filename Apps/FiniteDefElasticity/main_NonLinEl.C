// $Id: main_NonLinEl.C,v 1.6 2011-02-08 09:32:18 kmo Exp $
//==============================================================================
//!
//! \file main_NonLinEl.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric nonlinear elasticity solver.
//!
//==============================================================================

#include "NonLinSIM.h"
#include "SIMFiniteDefEl.h"
#include "LinAlgInit.h"
#include "Profiler.h"
#include "VTF.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Main program for the isogeometric nonlinear elasticity solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -skip2nd : Skip VTF-output of secondary solution fields
  \arg -saveInc \a dtSave : Time increment between each result save to VTF
  \arg -dumpInc \a dtDump [raw] : Time increment between each solution dump
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -2D : Use two-parametric simulation driver
  \arg -Tensor : Use tensorial total Lagrangian formulation (slow...)
  \arg -UL \a mver : Use updated Lagrangian formulation with nonlinear material
  \arg -MX \a pord : Mixed formulation with internal discontinuous pressure
  \arg -mixed : Mixed formulation with continuous pressure and volumetric change
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  SystemMatrix::Type solver = SystemMatrix::SPARSE;
  int form = SIM::TOTAL_LAGRANGE;
  int nGauss = 4;
  int format = 0;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  double dtSave = 0.0;
  double dtDump = 0.0;
  bool dumpWithID = true;
  bool skip2nd = false;
  bool checkRHS = false;
  bool fixDup = false;
  char twoD = false;
  char* infile = 0;

  LinAlgInit linalg(argc,argv);

  std::vector<int> options(2,-1);
  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-dense"))
      solver = SystemMatrix::DENSE;
    else if (!strcmp(argv[i],"-spr"))
      solver = SystemMatrix::SPR;
    else if (!strncmp(argv[i],"-superlu",8))
    {
      solver = SystemMatrix::SPARSE;
      if (isdigit(argv[i][8]))
	SIMbase::num_threads_SLU = atoi(argv[i]+8);
    }
    else if (!strcmp(argv[i],"-samg"))
      solver = SystemMatrix::SAMG;
    else if (!strcmp(argv[i],"-petsc"))
      solver = SystemMatrix::PETSC;
    else if (!strcmp(argv[i],"-nGauss") && i < argc-1)
      nGauss = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-vtf") && i < argc-1)
      format = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nviz") && i < argc-1)
      n[0] = n[1] = n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nu") && i < argc-1)
      n[0] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nv") && i < argc-1)
      n[1] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nw") && i < argc-1)
      n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-skip2nd"))
      skip2nd = true;
    else if (!strcmp(argv[i],"-saveInc") && i < argc-1)
      dtSave = atof(argv[++i]);
    else if (!strcmp(argv[i],"-dumpInc") && i < argc-1)
    {
      dtDump = atof(argv[++i]);
      if (++i < argc && !strcmp(argv[i],"raw"))
	dumpWithID = false;
      else
	--i;
    }
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
      {
	char* endp = 0;
	int endVal = 0;
	ignoredPatches.push_back(strtol(argv[++i],&endp,10));
	if (endp && *endp == ':')
	  endVal = strtol(endp+1,&endp,10);
	while (ignoredPatches.back() < endVal)
	  ignoredPatches.push_back(ignoredPatches.back()+1);
      }
    else if (!strcmp(argv[i],"-vox") && i < argc-1)
      VTF::vecOffset[0] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voy") && i < argc-1)
      VTF::vecOffset[1] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voz") && i < argc-1)
      VTF::vecOffset[2] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-check"))
      form += 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = strcmp(argv[i],"-2Dpstrain") ? 1 : 2;
    else if (!strcmp(argv[i],"-UL"))
    {
      if (form < SIM::UPDATED_LAGRANGE)
	form = SIM::UPDATED_LAGRANGE;
      if (i < argc-1 && isdigit(argv[i+1][0]))
        options[0] = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-MX"))
    {
      form = SIM::MIXED_QnPn1;
      if (i < argc-1 && isdigit(argv[i+1][0]))
	options[1] = atoi(argv[++i]);
      else
	options[1] = 0;
    }
    else if (!strcmp(argv[i],"-mixed"))
      form = SIM::MIXED_QnQn1;
    else if (!strcmp(argv[i],"-Tensor"))
      form = SIM::NONLINEAR;
    else if (!strcmp(argv[i],"-NeoHooke"))
      form = SIM::NEOHOOKE;
    else if (!strcmp(argv[i],"-NeoHookeIV"))
      form = SIM::NEOHOOKE_IV;
    else if (!strncmp(argv[i],"-lag",4))
      SIMbase::discretization = SIMbase::Lagrange;
    else if (!strncmp(argv[i],"-spec",5))
      SIMbase::discretization = SIMbase::Spectral;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-lag] [-spec] [-2D[pstrain]] [-UL [<mVER>]] [-MX [<p>]|-mixe"
	      <<"d] [-nGauss <n>]\n       [-vtf <format>] [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]\n      "
	      <<" [-saveInc <dtSave>] [-dumpInc <dtDump> [raw]]\n      "
	      <<" [-ignore <p1> <p2> ...] [-fixDup] [-checkRHS] [-check]\n";
    return 0;
  }

  if (linalg.myPid == 0)
    std::cout <<"\n >>> Spline FEM Nonlinear Elasticity solver <<<"
	      <<"\n ==============================================\n"
	      <<"\nInput file: "<< infile
	      <<"\nEquation solver: "<< solver
	      <<"\nNumber of Gauss points: "<< nGauss
	      <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
	      <<"\nNumber of visualization points: "
	      << n[0] <<" "<< n[1];
  if (twoD)
    n[2] = 1;
  else if (linalg.myPid == 0)
    std::cout <<" "<< n[2];

  if (linalg.myPid == 0)
  {
    if (dtSave > 0.0)
      std::cout <<"\nTime between each result save: "<< dtSave;
    if (dtDump > 0.0)
      std::cout <<"\nTime between each primary solution dump: "<< dtDump;
    if (SIMbase::discretization == SIMbase::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (SIMbase::discretization == SIMbase::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    if (fixDup)
      std::cout <<"\nCo-located nodes will be merged";
    if (checkRHS)
      std::cout <<"\nCheck that each patch has a right-hand coordinate system";
    if (!ignoredPatches.empty())
    {
      std::cout <<"\nIgnored patches:";
      for (size_t i = 0; i < ignoredPatches.size(); i++)
	std::cout <<" "<< ignoredPatches[i];
    }
    std::cout << std::endl;
  }
  utl::profiler->start("Model input");

  // Create the finite deformation elasticity model
  SIMbase* model;
  if (twoD)
    model = new SIMFiniteDefEl2D(form%100,twoD==1,options);
  else
    model = new SIMFiniteDefEl3D(checkRHS,form%100,options);

  // Read in solver and model definitions
  NonLinSIM simulator(model);
  if (!simulator.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  model->printProblem(std::cout);
  if (!model->preprocess(ignoredPatches,fixDup))
    return 2;

  // Save FE model to VTF file for visualization
  if (!simulator.saveModel(infile,format,n))
    return 3;

  std::ostream* oss = 0;
  if (dtDump > 0.0 && !dumpWithID)
  {
    // Write (refined?) model to g2-file
    strcat(strtok(infile,"."),".g2");
    std::cout <<"\nWriting updated g2-file "<< infile << std::endl;
    std::ofstream osg(infile);
    model->dumpGeometry(osg);
    // Open ASCII file for solution dump
    strcat(strtok(infile,"."),".sol");
    oss = new std::ofstream(infile);
    *oss <<"#NPoints="<< model->getNoNodes() <<"\n";
  }

  // Define the initial configuration
  NonLinSIM::SolvePrm params;
  simulator.init(params);

  const double epsT = 1.0e-6;
  if (dtDump <= 0.0) dtDump = params.stopTime + 1.0;
  double nextDump = params.time.t + dtDump;
  double nextSave = params.time.t + dtSave;

  int iStep = 0; // Save initial state to VTF
  if (params.multiSteps())
    if (!simulator.saveStep(-(++iStep),params.time.t,n,skip2nd))
      return 4;

  if (form >= 100)
    return 0; // model check

  // Initialize the linear solver
  model->initSystem(solver,1,1);
  model->setAssociatedRHS(0,0);
  model->setQuadratureRule(nGauss);

  // Invoke the time/load-step loop
  while (simulator.advanceStep(params))
  {
    // Solve the nonlinear FE problem at this load step
    if (!simulator.solveStep(params,SIM::STATIC))
      return 5;

    // Dump primary solution for inspection or external processing
    if (params.time.t + epsT*params.time.dt > nextDump)
    {
      if (dumpWithID)
	simulator.dumpStep(params.step,params.time.t,std::cout);
      else
	simulator.dumpStep(params.step,params.time.t,*oss,false);
      nextDump = params.time.t + dtDump;
    }

    // Save solution variables to VTF for visualization
    if (params.time.t + epsT*params.time.dt > nextSave)
      if (simulator.saveStep(++iStep,params.time.t,n,skip2nd))
	nextSave = params.time.t + dtSave;
      else
	return 6;
  }

  if (oss) delete oss;
  return 0;
}
