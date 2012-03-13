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
#include <cstring>
#include <cstdlib>


SIMoptions::SIMoptions () : nGauss({ 4, 4 }), nViz({ 2, 2, 2})
{
  discretization = ASM::Spline;
  solver = SystemMatrix::SPARSE;
  num_threads_SLU = 1;

  eig = 0;
  nev = 10;
  ncv = 20;
  shift = 0.0;

  format  = -1;
  saveInc =  1;
  dtSave  =  0.0;
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
  if (!strcasecmp(elem->Value(),"mode")) {
    if (elem->FirstChild())
      eig = atoi(elem->FirstChild()->Value());
  }

  else if (!strcasecmp(elem->Value(),"nev")) {
    if (elem->FirstChild())
      nev = atoi(elem->FirstChild()->Value());
  }
  else if (!strcasecmp(elem->Value(),"ncv")) {
    if (elem->FirstChild())
      ncv = atoi(elem->FirstChild()->Value());
  }

  else if (!strcasecmp(elem->Value(),"shift")) {
    if (elem->FirstChild())
      shift = atof(elem->FirstChild()->Value());
  }

  return true;
}


bool SIMoptions::parseDiscretizationTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"discretization")) {
    std::string discr;
    if (discretization == ASM::Spline &&
	utl::getAttribute(elem,"type",discr,true)) {
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

  else if (!strcasecmp(elem->Value(),"nGauss")) {
    if (elem->FirstChild()) {
      std::string value(elem->FirstChild()->Value());
      char* cval = strtok(const_cast<char*>(value.c_str())," ");
      for (int i = 0; i < 2 && cval; i++, cval = strtok(NULL," "))
	for (int j = i; j < 2; j++)
	  nGauss[j] = atoi(cval);
    }
    /* The above was not very elegant. Is there a better way of reading values
       from a node like <nGauss>2 4</nGauss> ? At least the below did not work
    if (elem->FirstChild())
      nGauss[0] = atoi(elem->FirstChild()->Value());
    if (elem->FirstChild()->NextSibling())
      nGauss[1] = atoi(elem->FirstChild()->NextSibling()->Value());
    else
      nGauss[1] = nGauss[0];
    */
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

  else if (!strcasecmp(elem->Value(),"stride")) {
    if (elem->FirstChild())
      saveInc = atoi(elem->FirstChild()->Value());
    utl::getAttribute(elem,"dt",dtSave);
  }

  else if (!strcasecmp(elem->Value(),"hdf5")) {
    if (elem->FirstChild())
    {
      hdf5 = elem->FirstChild()->Value();
      size_t pos = hdf5.find_last_of('.');
      if (pos < hdf5.size())
	hdf5.erase(pos);
    }
    else // use the default output file name
      hdf5 = "(default)";
  }

  else if (!strcasecmp(elem->Value(),"projection"))
    for (const TiXmlNode* ch = elem->FirstChild(); ch; ch = ch->NextSibling())
      if (!strcasecmp(ch->Value(),"global"))
	project[GLOBAL] = "Greville point projection";
      else if (!strcasecmp(ch->Value(),"dgl2"))
	project[DGL2] = "Discrete global L2-projection";
      else if (!strcasecmp(ch->Value(),"cgl2"))
	project[CGL2] = "Continuous global L2-projection";
      else if (!strcasecmp(ch->Value(),"scr"))
	project[SCR] = "Continuous global L2-projection";
      else if (!strcasecmp(ch->Value(),"vdsa"))
	project[VDSA] = "VDSA projected";
      else if (!strcasecmp(ch->Value(),"quasi"))
	project[QUASI] = "Quasi-interpolated";
      else if (!strcasecmp(ch->Value(),"lsq"))
	project[LEASTSQ] = "Least-square projected";

  return true;
}


bool SIMoptions::dumpHDF5 (char* defaultName)
{
  if (hdf5.empty()) return false;

  if (hdf5 == "(default)")
    hdf5 = strtok(defaultName,".");

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
  else if (!strcasecmp(argv[i],"-dgl2"))
    project[DGL2] = "Discrete global L2-projection";
  else if (!strcasecmp(argv[i],"-cgl2"))
    project[CGL2] = "Continuous global L2-projection";
  else if (!strcasecmp(argv[i],"-scr"))
    project[SCR]  = "Superconvergent recovery";
  else if (!strcasecmp(argv[i],"-vdsa"))
    project[VDSA] = "VDSA projected";
  else if (!strcasecmp(argv[i],"-quasi"))
    project[QUASI] = "Quasi-interpolated";
  else if (!strcasecmp(argv[i],"-lsq"))
    project[LEASTSQ] = "Least-square projected";
  else
    return false;

  return true;
}
