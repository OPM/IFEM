// $Id$
//==============================================================================
//!
//! \file main_NonLinEl.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric finite deformation solver.
//!
//==============================================================================

#include "SIMFiniteDefEl.h"
#include "NonLinSIM.h"
#include "ASMmxBase.h"
#include "LinAlgInit.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include "VTF.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

extern std::vector<int> mixedDbgEl; //!< List of elements for additional output


/*!
  \brief Main program for the isogeometric finite deformation solver.

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
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -skip2nd : Skip VTF/HDF5-output of secondary solution fields
  \arg -noEnergy : Skip integration of energy norms
  \arg -saveInc \a dtSave : Time increment between each result save to VTF/HDF5
  \arg -dumpInc \a dtDump [raw] : Time increment between each solution dump
  \arg -outPrec \a nDigit : Number of digits in solution component printout
  \arg -ztol \a eps : Zero tolerance for printing of solution norms
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -2D : Use two-parametric simulation driver (plane stress)
  \arg -2Dpstrain : Use two-parametric simulation driver (plane strain)
  \arg -UL : Use updated Lagrangian formulation with nonlinear material
  \arg -MX \a pord : Mixed formulation with internal discontinuous pressure
  \arg -mixed : Mixed formulation with continuous pressure and volumetric change
  \arg -Mixed : Same as -mixed, but use C^(p-1) continuous displacement basis
  \arg -Fbar \a nvp : Use the F-bar formulation
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  SystemMatrix::Type solver = SystemMatrix::SPARSE;
  int form = SIM::TOTAL_LAGRANGE;
  int nGauss = 4;
  int format = -1;
  int outPrec = 3;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  bool doHDF5 = false;
  double dtSave = 0.0;
  double dtDump = 0.0;
  double zero_tol = 1.0e-8;
  bool dumpWithID = true;
  bool skip2nd = false;
  bool energy = true;
  bool checkRHS = false;
  bool fixDup = false;
  bool twoD = false;
  char* infile = 0;

  const LinAlgInit& linalg = LinAlgInit::Init(argc,argv);

  std::vector<int> options(2,0);
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
      SIMbase::discretization = SIMbase::Lagrange;
    else if (!strncmp(argv[i],"-spec",5))
      SIMbase::discretization = SIMbase::Spectral;
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
    else if (!strcmp(argv[i],"-noEnergy"))
      energy = false;
    else if (!strcmp(argv[i],"-hdf5"))
      doHDF5 = true;
    else if (!strcmp(argv[i],"-saveInc") && i < argc-1)
      dtSave = atof(argv[++i]);
    else if (!strcmp(argv[i],"-outPrec") && i < argc-1)
      outPrec = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ztol") && i < argc-1)
      zero_tol = atof(argv[++i]);
    else if (!strcmp(argv[i],"-dumpInc") && i < argc-1)
    {
      dtDump = atof(argv[++i]);
      if (++i < argc && !strcmp(argv[i],"raw"))
	dumpWithID = false;
      else
	--i;
    }
    else if (!strcmp(argv[i],"-dbgElm"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(mixedDbgEl,argv[++i]);
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
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
    else if (!strcmp(argv[i],"-2Dpstrain"))
      twoD = SIMLinEl2D::planeStrain = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = true;
    else if (!strcmp(argv[i],"-UL"))
    {
      if (form < SIM::UPDATED_LAGRANGE)
	form = SIM::UPDATED_LAGRANGE;
    }
    else if (!strcmp(argv[i],"-MX"))
    {
      form = SIM::MIXED_QnPn1;
      if (i < argc-1 && isdigit(argv[i+1][0]))
	options[1] = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-Mixed"))
    {
      form = SIM::MIXED_QnQn1;
      ASMmxBase::useCpminus1 = true;
    }
    else if (!strcmp(argv[i],"-mixed"))
      form = SIM::MIXED_QnQn1;
    else if (!strcmp(argv[i],"-Fbar"))
    {
      form = SIM::FBAR;
      if (i < argc-1 && isdigit(argv[i+1][0]))
        options[1] = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-Tensor"))
      form = SIM::NONLINEAR;
    else if (!strcmp(argv[i],"-NeoHooke"))
      form = SIM::NEOHOOKE;
    else if (!strcmp(argv[i],"-NeoHookeIV"))
      form = SIM::NEOHOOKE_IV;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-lag] [-spec] [-nGauss <n>] [-2D[pstrain]] [-UL|-MX [<p>]"
	      <<"|-[M|m]ixed|-Fbar <nvp>]\n       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n      "
	      <<" [-saveInc <dtSave>] [-skip2nd] [-dumpInc <dtDump> [raw]]"
	      <<" [-outPrec <nd>]\n       [-ztol <eps>] [-ignore <p1> <p2> ...]"
	      <<" [-fixDup] [-checkRHS] [-check]\n";
    return 0;
  }

  if (linalg.myPid == 0)
  {
    std::cout <<"\n >>> Spline FEM Finite Deformation Nonlinear solver <<<"
	      <<"\n ======================================================\n"
	      <<"\n Executing command:\n";
    for (int i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< solver
	      <<"\nNumber of Gauss points: "<< nGauss;
    if (format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "
		<< n[0] <<" "<< n[1];
      if (!twoD) std::cout <<" "<< n[2];
    }
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

  // Create the finite deformation continuum model
  options[0] = form%100;
  SIMbase* model;
  if (twoD)
    model = new SIMFiniteDefEl2D(options);
  else
    model = new SIMFiniteDefEl3D(checkRHS,options);

  // Read in solver and model definitions
  NonLinSIM simulator(model);
  if (!simulator.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  model->printProblem(std::cout);
  if (SIMbase::discretization == SIMbase::Spline)
  {
    if (ASMmxBase::useCpminus1)
      std::cout <<"Using C^(p-1) continuous displacement basis\n";
    else if (form == SIM::MIXED_QnQn1)
      std::cout <<"Using C^(p-2) continuous displacement basis\n";
  }

  // Preprocess the model and establish data structures for the algebraic system
  if (!model->preprocess(ignoredPatches,fixDup))
    return 2;

  if (format >= 0)
  {
    // Save FE model to VTF file for visualization
    if (twoD) n[2] = 1;
    if (!simulator.saveModel(infile,format,n))
      return 3;
  }

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
  if (format >= 0 && params.multiSteps())
    if (!simulator.saveStep(-(++iStep),params.time.t,n,skip2nd))
      return 4;

  if (form >= 100)
    return 0; // model check

  DataExporter* writer = 0;
  if (doHDF5)
  {
    strtok(infile,".");
    if (linalg.myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< infile <<".hdf5"<< std::endl;

    writer = new DataExporter(true);
    writer->registerField("u","solution",DataExporter::SIM, skip2nd ? DataExporter::PRIMARY : DataExporter::SECONDARY);
    writer->setFieldValue("u",model,(void*)&simulator.getSolution());
    writer->registerWriter(new HDF5Writer(infile));
    writer->registerWriter(new XMLWriter(infile));
  }

  // Initialize the linear solver
  model->initSystem(solver,1,1);
  model->setAssociatedRHS(0,0);
  model->setQuadratureRule(nGauss);

  // Invoke the time/load-step loop
  while (simulator.advanceStep(params))
  {
    // Solve the nonlinear FE problem at this load step
    if (!simulator.solveStep(params,SIM::STATIC,"displacement",energy,
			     zero_tol, outPrec > 3 ? outPrec : 0))
      return 5;

    // Print solution components at the user-defined points
    simulator.dumpResults(params.time.t,std::cout,outPrec);

    if (params.time.t + epsT*params.time.dt > nextDump)
    {
      // Dump primary solution for inspection or external processing
      if (dumpWithID)
	simulator.dumpStep(params.step,params.time.t,std::cout);
      else
	simulator.dumpStep(params.step,params.time.t,*oss,false);

      nextDump = params.time.t + dtDump;
    }

    if (params.time.t + epsT*params.time.dt > nextSave)
    {
      // Save solution variables to VTF for visualization
      if (format >= 0)
        if (!simulator.saveStep(++iStep,params.time.t,n,skip2nd))
	  return 6;

      // Save solution variables to HDF5
      if (writer)
	writer->dumpTimeLevel();

      nextSave = params.time.t + dtSave;
    }
  }

  if (writer) delete writer;
  if (oss) delete oss;
  return 0;
}
