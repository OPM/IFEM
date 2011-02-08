// $Id: main_LinEl.C,v 1.20 2010-05-06 17:31:18 kmo Exp $
//==============================================================================
//!
//! \file main_LinEl.C
//!
//! \date Jan 28 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the NURBS-based linear elasticity solver.
//!
//==============================================================================

#include "LinearEl.h"
#include "VolumePatch.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Main program for the NURBS-based linear elasticity solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense : Use the dense LAPACK matrix equation solver
  \arg -superlu : Use the sparse superLU equation solver
  \arg -samg : Use the sparse algebraic multi-grid equation solver
  \arg -swapJ : Swap the sign of the Jacobian to account for inverse ordering
  \arg -splineMethod \a iop : Spline evaluation method (1 or 2)
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -dumpGeo : Dumps the (possibly refined) geometry to a g2-file
  \arg -dumpMat : Dumps the system matrices to file
  \arg -dumpDis : Dumps the solution vector to files, one for each component
  \arg -dumpGrid : Dumps a grid file with 8-noded hexahedron elements
  \arg -dumpGrid2 : Dumps a grid file with 27-noded hexahedron elements
*/

int main (const int argc, char** argv)
{
  SystemMatrix::Type solver = SystemMatrix::SPR;
  int nGauss = 4;
  int format = 0;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  int iop = 0;
  int nev = 10;
  int ncv = 20;
  double shf = 0.0;
  bool free = false;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpGeo = false;
  bool dumpMat = false;
  bool dumpDis = false;
  int dumpGrid = 0;
  char* infile = 0;

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-dense"))
      solver = SystemMatrix::DENSE;
    else if (!strcmp(argv[i],"-superlu"))
      solver = SystemMatrix::SPARSE;
    else if (!strcmp(argv[i],"-samg"))
      solver = SystemMatrix::SAMG;
    else if (!strcmp(argv[i],"-swapJ"))
      VolumePatch::swapJac = true;
    else if (!strcmp(argv[i],"-splineMethod") && i < argc-1)
      VolumePatch::splineEvalMethod = atoi(argv[++i]);
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
    else if (!strcmp(argv[i],"-shift"))
      shf = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      free = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-dumpGeo"))
      dumpGeo = true;
    else if (!strcmp(argv[i],"-dumpMat"))
      dumpMat = true;
    else if (!strcmp(argv[i],"-dumpDis"))
      dumpDis = true;
    else if (!strcmp(argv[i],"-dumpGrid"))
      dumpGrid = 1;
    else if (!strcmp(argv[i],"-dumpGrid2"))
      dumpGrid = 2;
    else if (!infile)
      infile = argv[i];

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0] <<" <inputfile> [-dense|-superlu|-samg]\n"
	      <<"       [-splineMethod <iop>] [-nGauss <n>] [-vtf <format>]\n"
	      <<"       [-nviz <nviz>] [-nu <nu>] [-nv <nv>] [-nw <nw>]\n"
	      <<"       [-eig <iop>] [-nev <nev>] [-ncv <ncv] [-shift <shf>]\n"
	      <<"       [-free] [-check] [-checkRHS] [-vizRHS] [-fixDup]\n"
	      <<"       [-ignore <p1> <p2> ...] [-swapJ] [-dumpGeo]\n"
	      <<"       [-dumpMat] [-dumpDis] [-dumpGrid|-dumpGrid2]\n";
    return 0;
  }


  // Load vector visualization is not available when using additional viz-points
  if (n[0] > 2 || n[1] > 2 || n[2] > 2) vizRHS = false;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (abs(iop) != 4 && abs(iop) != 6) free = false;

  std::cout <<"\n >>> ICADA Spline FEM prototype <<<"
	    <<"\n ==================================\n"
	    <<"\nInput file: "<< infile
	    <<"\nEquation solver: "<< solver
	    <<"\nspline evaluation method: "<< VolumePatch::splineEvalMethod
	    <<"\nNumber of Gauss points: "<< nGauss
	    <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
	    <<"\nNumber of visualization points: "
	    << n[0] <<" "<< n[1] <<" "<< n[2];
  if (iop != 0 && iop < 100)
    std::cout <<"\nEigenproblem solver: "<< iop
	      <<"\nNumber of eigenvalues: "<< nev
	      <<"\nNumber of Arnoldi vectors: "<< ncv
	      <<"\nShift value: "<< shf;
  if (free)
    std::cout <<"\nSpecified boundary conditions are ignored";
  if (VolumePatch::swapJac)
    std::cout <<"\nSwapped Jacobian determinant";
  if (fixDup)
    std::cout <<"\nCo-located nodes will be merged";
  if (checkRHS)
    std::cout <<"\nChecking that each patch has a right-hand coordinate system";
  if (!ignoredPatches.empty())
  {
    std::cout <<"\nIgnored patches:";
    for (size_t i = 0; i < ignoredPatches.size(); i++)
      std::cout <<" "<< ignoredPatches[i];
  }
  std::cout << std::endl;

  // Read in model definitions and establish the FE structures
  LinearEl model(infile,checkRHS,free);
  if (!model.preprocess(ignoredPatches,fixDup))
    return 1;

  if (dumpGeo) // Write the refined geometry to g2-file
    model.dumpGeometry(strcat(strtok(infile,"."),".g2"));

  Result solution;
  Matrix eNorm;

  switch (iop) {
  case 0:
  case 5:
  case 54:
  case 56:
    // Assemble and solve the linear system of equations
    if (!model.assembleKandR(solver,nGauss,vizRHS ? &solution.load : 0))
      return 2;
    else if (!model.solve(solution.displ))
      return 3;
    else if (iop == 0)
      if (!model.solutionNorms(solution.norms,nGauss,solution.displ))
	return 4;
      else
	break;

    // Linearized buckling: Assemble stiffness matrices [K] and [Kg]
    if (!model.assembleKandKg(solver,nGauss,solution.displ))
      return 2;
    else if (dumpMat)
      model.dumpMat("Km.mat","Kg.mat");
    else if (!model.modes(iop > 50 ? iop%50 : iop, nev,ncv,shf,solution.modes))
      return 3;
    break;

  case -1:
  case -2:
  case  1:
  case  2:
    // Assemble and solve the regular eigenvalue problem
    if (!model.assembleKonly(solver,nGauss))
      return 2;
    else if (!model.modes(iop,nev,ncv,shf,solution.modes))
      return 3;
    break;

  case 100:
  case 101:
  case 102:
    break; // Model check

  default:
    // Assemble and solve the generalized eigenvalue problem
    solution.freq = true;
    if (!model.assembleKandM(solver,nGauss))
      return 2;
    else if (!model.modes(iop,nev,ncv,shf,solution.modes))
      return 3;
  }

  // Write VTF-file with model and results
  if (dumpGrid == 0 || iop < 100)
    if (!model.writeGlv(infile,solution,n,format,iop==100)) return 4;

  // Write a global FE grid file for external processing
  if (dumpGrid > 0)
    if (!model.writeGlobalGrid(infile, n, dumpGrid == 2 ? 27 : 8)) return 5;

  // Write solution vectors to files for external viewing
  if (dumpDis)
  {
    strcat(strtok(infile,"."),"_dis");
    if (!solution.displ.empty())
      model.dumpSolution(infile,solution.displ);

    infile[strlen(infile)-3] = 'm';
    infile[strlen(infile)-2] = '0';
    infile[strlen(infile)-1] = '\0';
    for (size_t j = 1; j < solution.modes.size() && j < 10; j++)
    {
      infile[strlen(infile)-1]++;
      model.dumpSolution(infile,solution.modes[j-1].eigVec);
    }
  }

  return 0;
}
