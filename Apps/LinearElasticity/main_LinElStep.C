// $Id$
//==============================================================================
//!
//! \file main_LinElStep.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric linear elasticity solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMLinEl3D.h"
#include "LinAlgInit.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Main program for the NURBS-based isogeometric linear elasticity solver.

  This program is basicly the same as main_LinEl3D.C, except that this one runs
  in a loop over a specified number of load steps to compute the linear solution
  for each step, when the load function is assumed to vary with time.
  This version is for 3D patches only, and no eigenvalue solution is included.

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
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -nstep \a n : Number of time steps
  \arg -dt \a dt : Time increment size
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SystemMatrix::Type solver = SystemMatrix::SPARSE;
  int nGauss = 4;
  int format = 0;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  char* infile = 0;
  TimeDomain time;
  int nStep = 1;

  LinAlgInit linalg(argc,argv);

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-dense"))
      solver = SystemMatrix::DENSE;
    else if (!strcmp(argv[i],"-spr"))
      solver = SystemMatrix::SPR;
    else if (!strcmp(argv[i],"-superlu"))
      solver = SystemMatrix::SPARSE;
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
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-nstep") && i < argc-1)
      nStep = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-dt") && i < argc-1)
      time.dt = atof(argv[++i]);
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu|-samg|-petsc]\n"
	      <<"       [-nGauss <n>] [-nstep <n>] [-dt <dt>]\n"
	      <<"       [-vtf <format>] [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup] [-checkRHS]\n";
    return 0;
  }

  // Load vector visualization is not available when using additional viz-points
  if (n[0] > 2 || n[1] > 2 || n[2] > 2) vizRHS = false;

  std::cout <<"\n >>> IFEM Linear Elasticity solver <<<"
	    <<"\n ===========================================\n"
    InitIFEM(argc, argv);
  std::cout <<"\nInput file: "<< infile
	    <<"\nEquation solver: "<< solver
	    <<"\nNumber of Gauss points: "<< nGauss
	    <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
	    <<"\nNumber of visualization points: "
	    << n[0] <<" "<< n[1] <<" "<< n[2];
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

  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Read in model definitions and establish the FE structures
  SIMbase* model = new SIMLinEl3D(checkRHS);
  if (!model->read(infile) || !model->preprocess(ignoredPatches,fixDup))
    return 1;

  model->setQuadratureRule(nGauss);
  utl::profiler->stop("Model input");

  // Write VTF-file with model geometry
  int nBlock = 0;
  if (!model->writeGlv(infile,n,format))
    return 2;

  // Write Dirichlet boundary conditions
  if (!model->writeGlvBC(n,nBlock))
    return 3;

  Matrix eNorm;
  Vector gNorm, displ, load;

  model->initSystem(solver,1,1);
  model->setAssociatedRHS(0,0);

  for (int iStep = 1; iStep <= nStep; iStep++)
  {
    time.t += time.dt;
    std::cout <<"\n\nStep "<< iStep <<", time = "<< time.t << std::endl;
    model->setMode(SIM::STATIC);

    if (!model->assembleSystem(time))
      return 4;
    else if (vizRHS)
      model->extractLoadVec(load);

    // Solve the linear system of equations
    if (!model->solveSystem(displ,1))
      return 5;

    // Evaluate solution norms
    model->setMode(SIM::RECOVERY);
    if (!model->solutionNorms(Vectors(1,displ),eNorm,gNorm))
      return 6;
    std::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< gNorm(1);
    if (gNorm.size() > 1)
      std::cout <<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< gNorm(2);
    if (gNorm.size() > 2)
      std::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< gNorm(3);
    if (gNorm.size() > 3)
      std::cout <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< gNorm(4)
		<<"\nExact relative error (%) : "<< gNorm(4)/gNorm(3)*100.0;
    std::cout << std::endl;

    utl::profiler->start("Postprocessing");

    // Write boundary tractions, if any
    if (!model->writeGlvT(iStep,nBlock))
      return 7;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",n,iStep,nBlock))
      return 8;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(displ,n,iStep,nBlock))
      return 9;

    // Write element norms, if feasable
    if (n[0] == 2 && n[1] == 2 && n[2] == 2)
      if (!model->writeGlvN(eNorm,iStep,nBlock))
	return 10;

    // Write step information to VTF-file
    if (!model->writeGlvStep(iStep,time.t))
      return 11;

    utl::profiler->stop("Postprocessing");
  }

  model->closeGlv();
  delete model;

  return 0;
}
