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
#include "SIMLinElKL.h"
#include "SIMLinElBeamC1.h"
#include "AdaptiveSIM.h"
#include "LinAlgInit.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
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
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
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
  \arg -2Daxisymm : Use two-parametric simulation driver (axi-symmetric solid)
  \arg -KL : Use two-parametric simulation driver for Kirchhoff-Love plate
  \arg -Beam : Use one-parametric simulation driver for C1-continous beam
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
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
  bool dumpHDF5 = false;
  bool dumpASCII = false;
  bool twoD = false;
  bool KLp = false;
  bool Beam = false;
  char* infile = 0;

  const LinAlgInit& linalg = LinAlgInit::Init(argc,argv);

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
    else if (!strncmp(argv[i],"-lag",4))
      SIMbase::discretization = ASM::Lagrange;
    else if (!strncmp(argv[i],"-spec",5))
      SIMbase::discretization = ASM::Spectral;
    else if (!strncmp(argv[i],"-LR",3))
      SIMbase::discretization = ASM::LRSpline;
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
    else if (!strcmp(argv[i],"-hdf5"))
      dumpHDF5 = true;
    else if (!strcmp(argv[i],"-dumpASC"))
      dumpASCII = true;
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
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
    else if (!strcmp(argv[i],"-Beam"))
      Beam = true;
    else if (!strcmp(argv[i],"-KL"))
    {
      twoD = KLp = true;
      SIMbase::discretization = ASM::SplineC1;
    }
    else if (!strncmp(argv[i],"-2Dpstra",8))
      twoD = SIMLinEl2D::planeStrain = true;
    else if (!strncmp(argv[i],"-2Daxi",6))
      twoD = SIMLinEl2D::axiSymmetry = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = true;
    else if (!strncmp(argv[i],"-adap",5))
    {
      SIMbase::discretization = ASM::LRSpline;
      iop = 10;
    }
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-free] [-lag] [-spec] [-LR] [-2D[pstrain|axisymm]] [-KL|-Bea"
	      <<"m] [-adap] [-nGauss <n>]\n       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
	      <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup]"
	      <<" [-checkRHS] [-check] [-dumpASC]\n";
    return 0;
  }

  if (Beam) n[1] = n[2] = 1;
  if (twoD) n[2] = 1;
  // Load vector visualization is not available when using additional viz-points
  if (n[0] > 2 || n[1] > 2 || n[2] > 2) vizRHS = false;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (iop != 4 && iop != 6) SIMbase::ignoreDirichlet = false;

  if (linalg.myPid == 0)
  {
    std::cout <<"\n >>> Spline FEM Linear Elasticity solver <<<"
	      <<"\n ===========================================\n"
	      <<"\n Executing command:\n";
    for (int i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< solver
	      <<"\nNumber of Gauss points: "<< nGauss;
    if (format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< n[0];
      if (!Beam)
      {
	std::cout <<" "<< n[1];
	if (!twoD) std::cout <<" "<< n[2];
      }
    }

    if (iop > 0 && iop < 10)
      std::cout <<"\nEigenproblem solver: "<< iop
		<<"\nNumber of eigenvalues: "<< nev
		<<"\nNumber of Arnoldi vectors: "<< ncv
		<<"\nShift value: "<< shf;
    if (SIMbase::discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (SIMbase::discretization == ASM::Spectral)
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
  if (Beam)
    model = new SIMLinElBeamC1();
  else if (KLp)
    model = new SIMLinElKL();
  else if (twoD)
    model = new SIMLinEl2D();
  else
    model = new SIMLinEl3D(checkRHS);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = 0;
  if (iop == 10)
    theSim = aSim = new AdaptiveSIM(model);

  if (!theSim->read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (linalg.myPid == 0)
    model->printProblem(std::cout);

  if (!model->preprocess(ignoredPatches,fixDup))
    return 1;

  model->setQuadratureRule(nGauss);

  Matrix eNorm, ssol;
  Vector gNorm, load;
  Vectors displ(1), projs;
  std::vector<Mode> modes;
  std::vector<Mode>::const_iterator it;
  int iStep = 1, nBlock = 0;

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

    model->setMode(SIM::RECOVERY);
    if (SIMbase::discretization == ASM::Spline ||
	SIMbase::discretization == ASM::SplineC1)
    {
      // Project the FE stresses onto the splines basis
      if (!model->project(ssol,displ.front(),SIMbase::GLOBAL))
	return 4;
      else
	projs.push_back(ssol);
    }

    // Integrate solution norms and errors
    if (!model->solutionNorms(displ,projs,eNorm,gNorm))
      return 4;

    if (linalg.myPid == 0)
    {
      std::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< gNorm(1);
      std::cout	<<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< gNorm(2);
      if (model->haveAnaSol() && gNorm.size() >= 4)
	std::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< gNorm(3)
		  <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< gNorm(4)
		  <<"\nExact relative error (%) : "<< gNorm(4)/gNorm(3)*100.0;
      size_t j = model->haveAnaSol() ? 5 : 3;
      for (size_t i = 0; i < projs.size() && j < gNorm.size(); i++)
      {
        std::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm(j++);
	std::cout <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm(j++);
	std::cout <<"\n- relative error (% of |u^r|) : "
		  << gNorm(j-1)/gNorm(j-2)*100.0;
	if (model->haveAnaSol() && j++ <= gNorm.size())
	  std::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm(j-1)
		    <<"\n- relative error (% of |u|)   : "
		    << gNorm(j-1)/gNorm(3)*100.0;
      }
      std::cout << std::endl;

      model->dumpResults(displ.front(),0.0,std::cout,true,6);
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

  case 10:
    // Adaptive simulation
    while (true)
      if (!aSim->solveStep(infile,solver,iStep))
        return 5;
      else if (!aSim->writeGlv(infile,format,n,iStep,nBlock))
	return 6;
      else if (!aSim->adaptMesh(++iStep))
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

  if (dumpHDF5 && !displ.front().empty())
  {
    strtok(infile,".");
    if (linalg.myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< infile <<".hdf5"<< std::endl;
    DataExporter exporter(true);
    exporter.registerField("u","solution",DataExporter::SIM,
			   DataExporter::SECONDARY);
    exporter.setFieldValue("u",model,&displ.front());
    exporter.registerWriter(new HDF5Writer(infile));
    exporter.registerWriter(new XMLWriter(infile));
    exporter.dumpTimeLevel();
  }

  if (iop != 10 && format >= 0)
  {
    // Write VTF-file with model geometry
    if (!model->writeGlvG(n,nBlock,infile,format))
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

    // Write projected solution fields to VTF-file
    if (projs.size() > 0)
      if (!model->writeGlvP(projs.front(),n,iStep,nBlock))
	return 11;

    // Write eigenmodes
    for (it = modes.begin(); it != modes.end(); it++)
      if (!model->writeGlvM(*it, iop==3 || iop==4 || iop==6, n, nBlock))
	return 12;

    // Write element norms
    if (!model->writeGlvN(eNorm,iStep,nBlock))
      return 13;

    model->writeGlvStep(1);
  }
  model->closeGlv();

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
  delete theSim;
  return 0;
}
