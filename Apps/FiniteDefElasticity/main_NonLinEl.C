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
#include "NonlinearDriver.h"
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

#ifndef USE_OPENMP
extern std::vector<int> mixedDbgEl; //!< List of elements for additional output
#endif


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
  \arg -2Daxisymm : Use two-parametric simulation driver (axi-symmetric solid)
  \arg -UL : Use updated Lagrangian formulation with nonlinear material
  \arg -MX \a pord : Mixed formulation with internal discontinuous pressure
  \arg -mixed : Mixed formulation with continuous pressure and volumetric change
  \arg -Mixed : Same as -mixed, but use C^(p-1) continuous displacement basis
  \arg -Fbar \a nvp : Use the F-bar formulation
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  SIMoptions dummy;
  std::vector<int> options(2,0), ignoredPatches;
  int i, form = SIM::TOTAL_LAGRANGE;
  int outPrec = 3;
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

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-skip2nd"))
      skip2nd = true;
    else if (!strcmp(argv[i],"-noEnergy"))
      energy = false;
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
#ifndef USE_OPENMP
    else if (!strcmp(argv[i],"-dbgElm"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(mixedDbgEl,argv[++i]);
#endif
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
    else if (!strncmp(argv[i],"-2Dpstra",8))
      twoD = SIMLinEl2D::planeStrain = true;
    else if (!strncmp(argv[i],"-2Daxi",6))
      twoD = SIMLinEl2D::axiSymmetry = true;
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
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-lag|-spec] [-nGauss <n>] [-2D[pstrain|axis]] [-UL|-MX [<p>]"
	      <<"|-[M|m]ixed|-Fbar <nvp>]\n       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n      "
	      <<" [-saveInc <dtSave>] [-skip2nd] [-dumpInc <dtDump> [raw]]"
	      <<" [-outPrec <nd>]\n       [-ztol <eps>] [-ignore <p1> <p2> ...]"
	      <<" [-fixDup] [-checkRHS] [-check]\n";
    return 0;
  }

  if (linalg.myPid == 0)
  {
    std::cout <<"\n >>> IFEM Finite Deformation Nonlinear solver <<<"
	      <<"\n ================================================\n"
	      <<"\n Executing command:\n";
    for (i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< dummy.solver
	      <<"\nNumber of Gauss points: "<< dummy.nGauss[0];
    if (dummy.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (dummy.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "
		<< dummy.nViz[0] <<" "<< dummy.nViz[1];
      if (!twoD) std::cout <<" "<< dummy.nViz[2];
    }
    if (dummy.dtSave > 0.0)
      std::cout <<"\nTime between each result save: "<< dummy.dtSave;
    if (dtDump > 0.0)
      std::cout <<"\nTime between each primary solution dump: "<< dtDump;
    if (dummy.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (dummy.discretization == ASM::Spectral)
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
    if (outPrec != 3)
      std::cout <<"\nNorm- and component output precision: "<< outPrec;
    if (zero_tol != 1.0e-8)
      std::cout <<"\nNorm output zero tolerance: "<< zero_tol;
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
  model->opt.discretization = dummy.discretization;
  NonlinearDriver simulator(model);
  if (!simulator.read(infile))
    return 1;

  // Parse the obsolete options again to let them override input file tags
  dummy.discretization = model->opt.discretization; // but not this option
  for (i = 1; i < argc; i++)
    if (!model->opt.parseOldOptions(argc,argv,i))
      if (!strcmp(argv[i],"-ignore"))
	while (i < argc-1 && isdigit(argv[i+1][0])) ++i;
  model->opt.discretization = dummy.discretization; // XML-tag is used, if set

  if (linalg.myPid == 0)
  {
    std::cout <<"\n\nEquation solver: "<< model->opt.solver
	      <<"\nNumber of Gauss points: "<< model->opt.nGauss[0];
    if (model->opt.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (model->opt.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "
		<< model->opt.nViz[0] <<" "<< model->opt.nViz[1];
      if (!twoD) std::cout <<" "<< model->opt.nViz[2];
    }
    if (model->opt.dtSave > 0.0)
      std::cout <<"\nTime between each result save: "<< model->opt.dtSave;
    if (model->opt.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (model->opt.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    std::cout << std::endl;

    model->printProblem(std::cout);
    if (model->opt.discretization >= ASM::Spline)
    {
      if (ASMmxBase::useCpminus1)
	std::cout <<"Using C^(p-1) continuous displacement basis\n";
      else if (form == SIM::MIXED_QnQn1)
	std::cout <<"Using C^(p-2) continuous displacement basis\n";
    }
  }

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  if (!model->preprocess(ignoredPatches,fixDup))
    return 2;

  if (model->opt.format >= 0)
  {
    // Save FE model to VTF file for visualization
    if (twoD) model->opt.nViz[2] = 1;
    if (!simulator.saveModel(infile))
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
    if (form < 100)
    {
      strcat(strtok(infile,"."),".sol");
      oss = new std::ofstream(infile);
      *oss <<"#NPoints="<< model->getNoNodes() <<"\n";
    }
  }

  if (form >= 100)
    return 0; // model check

  // Define the initial configuration
  simulator.init();

  DataExporter* writer = NULL;
  if (model->opt.dumpHDF5(infile))
  {
    // Opend HDF5 result database
    if (linalg.myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< model->opt.hdf5
		<<".hdf5"<< std::endl;
    writer = new DataExporter(true);
    writer->registerField("u","solution", DataExporter::SIM, skip2nd ?
			  DataExporter::PRIMARY : DataExporter::SECONDARY);
    writer->setFieldValue("u",model,&simulator.getSolution());
    writer->registerWriter(new HDF5Writer(model->opt.hdf5));
    writer->registerWriter(new XMLWriter(model->opt.hdf5));
  }

  // Now invoke the main solution driver
  int status = simulator.solveProblem(skip2nd,energy,writer,oss,
				      dtDump,zero_tol,outPrec);

  if (writer) delete writer;
  if (oss) delete oss;
  return status;
}
