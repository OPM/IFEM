// $Id$
//==============================================================================
//!
//! \file SIMoptions.C
//!
//! \date Feb 13 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for encapsulation of general simulation options.
//!
//==============================================================================

#include "SIMoptions.h"
#include "SystemMatrix.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "IFEM.h"
#include "LogStream.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <cstring>
#include <cstdlib>
#include <fstream>


SIMoptions::SIMoptions ()
{
  discretization = ASM::Spline;
  solver = SystemMatrix::SPARSE;
#ifdef USE_OPENMP
  num_threads_SLU = omp_get_max_threads();
#else
  num_threads_SLU = 1;
#endif

  eig = 0;
  nev = 10;
  ncv = 20;
  shift = 0.0;

  format  = -1;
  saveInc =  1;
  dtSave  =  0.0;
  pSolOnly = false;
  enableController = false;

  nGauss[0] = nGauss[1] = 4;
  nViz[0] = nViz[1] = nViz[2] = 2;

  printPid = 0;
}


void SIMoptions::setLinearSolver (const std::string& eqsolver)
{
  if (eqsolver == "dense")
    solver = SystemMatrix::DENSE;
  else if (eqsolver == "spr")
    solver = SystemMatrix::SPR;
  else if (eqsolver == "superlu")
    solver = SystemMatrix::SPARSE;
  else if (eqsolver == "samg")
    solver = SystemMatrix::SAMG;
  else if (eqsolver == "petsc")
    solver = SystemMatrix::PETSC;
}


bool SIMoptions::parseEigSolTag (const TiXmlElement* elem)
{
  const char* value;
  if ((value = utl::getValue(elem,"mode")))
    eig = atoi(value);
  else if ((value = utl::getValue(elem,"nev")))
    nev = atoi(value);
  else if ((value = utl::getValue(elem,"ncv")))
    ncv = atoi(value);
  else if ((value = utl::getValue(elem,"shift")))
    shift = atof(value);

  return true;
}


bool SIMoptions::parseDiscretizationTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"discretization")) {
    std::string discr;
    if (utl::getAttribute(elem,"type",discr,true)) {
      if (discr == "lagrange")
	discretization = ASM::Lagrange;
      else if (discr == "spectral")
	discretization = ASM::Spectral;
      else if (discr == "splines")
	discretization = ASM::Spline;
      else if (discr == "lrsplines")
	discretization = ASM::LRSpline;
    }
  }

  else if (!strcasecmp(elem->Value(),"nGauss") && elem->FirstChild()) {
    std::string value(elem->FirstChild()->Value());
    char* cval = strtok(const_cast<char*>(value.c_str())," ");
    for (int i = 0; i < 2 && cval; i++, cval = strtok(nullptr," "))
      for (int j = i; j < 2; j++)
        nGauss[j] = atoi(cval);
  }

  return true;
}


bool SIMoptions::parseOutputTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"vtfformat")) {
    if (elem->FirstChild()) {
      if (!strcasecmp(elem->FirstChild()->Value(),"ascii"))
	format = 0;
      else if (!strcasecmp(elem->FirstChild()->Value(),"binary"))
	format = 1;
    }
    if (utl::getAttribute(elem,"nviz",nViz[0]))
      nViz[2] = nViz[1] = nViz[0];
    utl::getAttribute(elem,"nu",nViz[0]);
    utl::getAttribute(elem,"nv",nViz[1]);
    utl::getAttribute(elem,"nw",nViz[2]);
  }

  else if (!strcasecmp(elem->Value(),"logging")) {
    int pid = 0;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
#endif
    utl::getAttribute(elem,"output_pid",printPid);
    if (printPid != -1 && printPid != IFEM::getOptions().printPid) {
      IFEM::getOptions().printPid = printPid;
      IFEM::cout.setPIDs(printPid,pid);
      if (printPid != pid)
        IFEM::cout.setNull();
      IFEM::cout <<"IFEM: Printing output from PID "<< printPid
                 <<" to console."<< std::endl;
    }
    utl::getAttribute(elem,"output_prefix",log_prefix);
    if (!log_prefix.empty() && log_prefix != IFEM::getOptions().log_prefix) {
      if ((pid == 0 && printPid == -1) || pid == IFEM::getOptions().printPid)
        IFEM::cout <<"IFEM: Logging output to files with prefix "
                   << log_prefix << std::endl;
      IFEM::getOptions().log_prefix = log_prefix;
      char cPid[12];
      sprintf(cPid,"_p%04d.log",pid);
      IFEM::cout.addExtraLog(new std::ofstream(log_prefix+cPid));
    }
  }

  else if (!strcasecmp(elem->Value(),"stride")) {
    const char* value = utl::getValue(elem,"stride");
    if (value) saveInc = atoi(value);
    utl::getAttribute(elem,"dt",dtSave);
  }

  else if (!strcasecmp(elem->Value(),"hdf5")) {
    if (elem->FirstChild()) {
      hdf5 = elem->FirstChild()->Value();
      size_t pos = hdf5.find_last_of('.');
      if (pos < hdf5.size())
	hdf5.erase(pos);
    }
    else // use the default output file name
      hdf5 = "(default)";
  }

  else if (!strcasecmp(elem->Value(),"primarySolOnly"))
    pSolOnly = true;

  else if (!strcasecmp(elem->Value(),"projection")) {
    std::string type;
    if (utl::getAttribute(elem,"type",type))
      this->parseProjectionMethod(type.c_str());
    for (const TiXmlNode* ch = elem->FirstChild(); ch; ch = ch->NextSibling())
      this->parseProjectionMethod(ch->Value());
  }

  return true;
}


bool SIMoptions::dumpHDF5 (const char* defaultName)
{
  if (hdf5.empty()) return false;

  if (hdf5 == "(default)") {
    hdf5 = defaultName;
    hdf5.erase(hdf5.find_last_of("."));
  }

  return true;
}


/*!
  These options may also be specified as tags on the XML input file. However,
  specified command-line options will override similar option on the input file.
*/

bool SIMoptions::parseOldOptions (int argc, char** argv, int& i)
{
  if (!strcmp(argv[i],"-dense"))
    solver = SystemMatrix::DENSE;
  else if (!strcmp(argv[i],"-spr"))
    solver = SystemMatrix::SPR;
  else if (!strncmp(argv[i],"-superlu",8))
  {
    solver = SystemMatrix::SPARSE;
    if (isdigit(argv[i][8]))
      num_threads_SLU = atoi(argv[i]+8);
  }
  else if (!strcmp(argv[i],"-samg"))
    solver = SystemMatrix::SAMG;
  else if (!strcmp(argv[i],"-petsc"))
    solver = SystemMatrix::PETSC;
  else if (!strncmp(argv[i],"-lag",4))
    discretization = ASM::Lagrange;
  else if (!strncmp(argv[i],"-spec",5))
    discretization = ASM::Spectral;
  else if (!strncmp(argv[i],"-LR",3))
    discretization = ASM::LRSpline;
  else if (!strcmp(argv[i],"-nGauss") && i < argc-1)
    nGauss[0] = nGauss[1] = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-vtf") && i < argc-1)
    format = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-nviz") && i < argc-1)
    nViz[0] = nViz[1] = nViz[2] = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-nu") && i < argc-1)
    nViz[0] = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-nv") && i < argc-1)
    nViz[1] = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-nw") && i < argc-1)
    nViz[2] = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-hdf5"))
  {
    if (i < argc-1 && argv[i+1][0] != '-')
      hdf5 = strtok(argv[++i],".");
    else // use the default output file name
      hdf5 = "(default)";
  }
  else if (!strcmp(argv[i],"-saveInc") && i < argc-1)
    dtSave = atof(argv[++i]);
  else if (!strcmp(argv[i],"-eig") && i < argc-1)
    eig = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-nev") && i < argc-1)
    nev = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-ncv") && i < argc-1)
    ncv = atoi(argv[++i]);
  else if (!strcmp(argv[i],"-shift") && i < argc-1)
    shift = atof(argv[++i]);
  else if (!strcasecmp(argv[i],"-controller"))
    enableController = true;
  else if (argv[i][0] == '-')
    return this->parseProjectionMethod(argv[i]+1);
  else
    return false;

  return true;
}


bool SIMoptions::parseProjectionMethod (const char* ptype)
{
  if (!strcasecmp(ptype,"global") || !strcasecmp(ptype,"grvl"))
    project[GLOBAL] = "Greville point projection";
  else if (!strcasecmp(ptype,"dgl2"))
    project[DGL2] = "Discrete global L2-projection";
  else if (!strcasecmp(ptype,"cgl2"))
    project[CGL2] = "Continuous global L2-projection";
  else if (!strcasecmp(ptype,"scr"))
    project[SCR] = "Superconvergent recovery";
  else if (!strcasecmp(ptype,"vdsa"))
    project[VDSA] = "VDSA projected";
  else if (!strcasecmp(ptype,"quasi"))
    project[QUASI] = "Quasi-interpolated";
  else if (!strcasecmp(ptype,"lsq"))
    project[LEASTSQ] = "Least-square projected";
  else
    return false;

  return true;
}


bool SIMoptions::ignoreOldOptions (int argc, char** argv, int& i)
{
  static SIMoptions dummy;
  return dummy.parseOldOptions(argc,argv,i);
}


utl::LogStream& SIMoptions::print (utl::LogStream& os, bool addBlankLine) const
{
  if (addBlankLine) os <<"\n";

  os <<"\nEquation solver: "<< solver;

  if (eig > 0)
    os <<"\nEigenproblem solver: "<< eig
       <<"\nNumber of eigenvalues: "<< nev
       <<"\nNumber of Arnoldi vectors: "<< ncv
       <<"\nShift value: "<< shift;

  os <<"\nNumber of Gauss points: "<< nGauss[0];
  if (nGauss[1] != nGauss[0]) os <<" "<< nGauss[1];

  switch (discretization) {
  case ASM::Lagrange:
    os <<"\nLagrangian basis functions are used"; break;
  case ASM::Spectral:
    os <<"\nSpectral basis functions are used"; break;
  case ASM::LRSpline:
    os <<"\nLR-spline basis functions are used"; break;
  case ASM::SplineC1:
    os <<"\nSpline basis with C1-continuous patch interfaces is used"; break;
  default: break;
  }

  if (!project.empty()) {
    ProjectionMap::const_iterator it = project.begin();
    os <<"\nEnabled projection(s): "<< it->second;
    for (++it; it != project.end(); ++it)
      os <<"\n                       "<< it->second;
  }

  if (format >= 0) {
    os <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
       <<"\nNumber of visualization points: "<< nViz[0];
    for (int j = 1; j < 3 && nViz[j] > 1; j++) os <<" "<< nViz[j];
  }

  if (!hdf5.empty())
    os <<"\nHDF5 result database: "<< hdf5 <<".hdf5";
  else if (format < 0)
    return os;

  if (dtSave > 0.0)
    os <<"\nTime between each result save: "<< dtSave;
  if (saveInc > 1)
    os <<"\nIncrements between each result save: "<< saveInc;

  if (pSolOnly)
    os <<"\nSecondary solution variables are not saved.";

  return os;
}
