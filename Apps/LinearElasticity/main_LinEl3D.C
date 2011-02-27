// $Id$
//==============================================================================
//!
//! \file main_LinEl3D.C
//!
//! \date Jan 28 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric linear elasticity solver.
//!
//==============================================================================

#include "SIMLinEl3D.h"
#include "SIMLinEl2D.h"
#include "LinAlgInit.h"
#include "Profiler.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Main program for the NURBS-based isogeometric linear elasticity solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -dumpASC : Dump model and solution to ASCII files for external processing
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -2D : Use two-parametric simulation driver (plane stress)
  \arg -2Dpstrain : Use two-parametric simulation driver (plane strain)
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SystemMatrix::Type solver = SystemMatrix::SPARSE;
  int nGauss = 4;
  int format = -1;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  int iop = 0;
  int nev = 10;
  int ncv = 20;
  double shf = 0.0;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpASCII = false;
  char twoD = false;
  char* infile = 0;

  LinAlgInit linalg(argc,argv);

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
    else if (!strcmp(argv[i],"-dumpASC"))
      dumpASCII = true;
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
    else if (!strcmp(argv[i],"-eig") && i < argc-1)
      iop = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nev") && i < argc-1)
      nev = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ncv") && i < argc-1)
      ncv = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-shift") && i < argc-1)
      shf = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = strcmp(argv[i],"-2Dpstrain") ? 1 : 2;
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
	      <<" <inputfile> [-dense|-spr|-superlu<nt>|-samg|-petsc]\n"
	      <<"       [-free] [-lag] [-spec] [-2D[pstrain]] [-nGauss <n>]\n"
	      <<"       [-vtf <format>] [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]\n"
	      <<"       [-eig <iop>] [-nev <nev>] [-ncv <ncv] [-shift <shf>]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup]"
	      <<" [-checkRHS] [-check] [-dumpASC]\n";
    return 0;
  }

  if (twoD) n[2] = 1;
  // Load vector visualization is not available when using additional viz-points
  if (n[0] > 2 || n[1] > 2 || n[2] > 2) vizRHS = false;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (iop != 4 && iop != 6) SIMbase::ignoreDirichlet = false;

  if (linalg.myPid == 0)
  {
    std::cout <<"\n >>> Spline FEM Linear Elasticity solver <<<"
	      <<"\n ===========================================\n"
	      <<"\nInput file: "<< infile
	      <<"\nEquation solver: "<< solver
	      <<"\nNumber of Gauss points: "<< nGauss;
    if (format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "
		<< n[0] <<" "<< n[1];
      if (!twoD) std::cout <<" "<< n[2];
    }

    if (iop > 0 && iop < 100)
      std::cout <<"\nEigenproblem solver: "<< iop
		<<"\nNumber of eigenvalues: "<< nev
		<<"\nNumber of Arnoldi vectors: "<< ncv
		<<"\nShift value: "<< shf;
    if (SIMbase::discretization == SIMbase::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (SIMbase::discretization == SIMbase::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    if (SIMbase::ignoreDirichlet)
      std::cout <<"\nSpecified boundary conditions are ignored";
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
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Read in model definitions and establish the FE data structures
  SIMbase* model;
  if (twoD)
    model = new SIMLinEl2D(SIM::LINEAR,twoD==1);
  else
    model = new SIMLinEl3D(checkRHS);

  if (!model->read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (linalg.myPid == 0)
    model->printProblem(std::cout);

  if (!model->preprocess(ignoredPatches,fixDup))
    return 1;

  model->setQuadratureRule(nGauss);

  Matrix eNorm;
  Vector gNorm, load;
  Vectors displ(1);
  std::vector<Mode> modes;
  std::vector<Mode>::const_iterator it;

  switch (iop) {
  case 0:
  case 5:
    // Static solution: Assemble [Km] and {R}
    model->setMode(SIM::STATIC);
    model->initSystem(solver,1,1);
    model->setAssociatedRHS(0,0);
    if (!model->assembleSystem())
      return 2;
    else if (vizRHS)
      model->extractLoadVec(load);

    // Solve the linear system of equations
    if (!model->solveSystem(displ.front(),1))
      return 3;

    // Evaluate solution norms
    model->setMode(SIM::RECOVERY);
    if (!model->solutionNorms(displ,eNorm,gNorm))
      return 4;

    if (linalg.myPid == 0)
    {
      std::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< gNorm(1);
      std::cout <<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< gNorm(2);
      if (gNorm.size() > 2)
	std::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< gNorm(3);
      if (gNorm.size() > 3)
	std::cout <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< gNorm(4)
		  <<"\nExact relative error (%) : "<< gNorm(4)/gNorm(3)*100.0;
      std::cout << std::endl;
    }

    if (iop == 0) break;

    // Linearized buckling: Assemble [Km] and [Kg]
    model->setMode(SIM::BUCKLING);
    model->initSystem(solver,2,0);
    if (!model->assembleSystem(displ))
      return 5;

    // Solve the generalized eigenvalue problem
    if (!model->systemModes(modes,nev,ncv,iop,shf))
      return 6;
    break;

  case 1:
  case 2:
    // Assemble and solve the regular eigenvalue problem
    model->setMode(SIM::STIFF_ONLY);
    model->initSystem(solver,1,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes,nev,ncv,iop,shf))
      return 6;
    break;

  case 100:
    break; // Model check

  default:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->initSystem(solver,2,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes,nev,ncv,iop,shf))
      return 6;
  }

  utl::profiler->start("Postprocessing");

  if (format >= 0)
  {
    // Write VTF-file with model geometry
    int iStep = 1, nBlock = 0;
    if (!model->writeGlv(infile,n,format))
      return 7;

    // Write boundary tractions, if any
    if (!model->writeGlvT(iStep,nBlock))
      return 8;

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(n,nBlock))
      return 8;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",n,iStep,nBlock))
      return 9;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(displ.front(),n,iStep,nBlock))
      return 10;

    // Write eigenmodes
    for (it = modes.begin(); it != modes.end(); it++)
      if (!model->writeGlvM(*it, iop==3 || iop==4 || iop==6, n, nBlock))
	return 11;

    // Write element norms (when no additional visualization points are used)
    if (n[0] == 2 && n[1] == 2 && n[2] <= 2)
      if (!model->writeGlvN(eNorm,iStep,nBlock))
	return 12;

    model->closeGlv();
  }

  if (dumpASCII)
  {
    // Write (refined) model to g2-file
    std::ofstream osg(strcat(strtok(infile,"."),".g2"));
    std::cout <<"\nWriting updated g2-file "<< infile << std::endl;
    model->dumpGeometry(osg);
    if (!displ.front().empty())
    {
      // Write solution (control point values) to ASCII files
      std::ofstream osd(strcat(strtok(infile,"."),".dis"));
      std::cout <<"\nWriting deformation to file "<< infile << std::endl;
      model->dumpPrimSol(displ.front(),osd,false);
      std::ofstream oss(strcat(strtok(infile,"."),".sol"));
      std::cout <<"\nWriting solution to file "<< infile << std::endl;
      model->dumpSolution(displ.front(),oss);
    }
    if (!modes.empty())
    {
      // Write eigenvectors to ASCII files
      std::ofstream ose(strcat(strtok(infile,"."),".eig"));
      std::cout <<"\nWriting eigenvectors to file "<< infile << std::endl;
      for (it = modes.begin(); it != modes.end(); it++)
      {
	ose <<"# Eigenvector_"<< it->eigNo <<" Eigenvalue="<< it->eigVal <<"\n";
	model->dumpPrimSol(it->eigVec,ose,false);
      }
    }
  }

  utl::profiler->stop("Postprocessing");
  delete model;
  return 0;
}
