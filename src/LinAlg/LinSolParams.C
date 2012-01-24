// $Id$
//==============================================================================
//!
//! \file LinSolParams.C
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//!
//==============================================================================

#include "LinSolParams.h"
#include "Utilities.h"
#include <fstream>
#ifdef HAS_PETSC
#include "petscversion.h"

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
#include "petscpcmg.h"
#define DEFAULTSOLVER MATSOLVERPETSC
#else
#include "petscmg.h"
#define DEFAULTSOLVER MAT_SOLVER_PETSC
#endif

#endif
#include "tinyxml.h"


void LinSolParams::setDefault ()
{
#ifdef HAS_PETSC
  // Use GMRES with ILU preconditioner as default
  method    = KSPGMRES;
  hypretype = "boomeramg";
  prec      = PCILU;
  package   = DEFAULTSOLVER;
  levels    = 0;
  overlap   = 0;
  nullspc   = NONE;

  atol      = 1.0e-6;
  rtol      = 1.0e-6;
  dtol      = 1.0e3;
  maxIts    = 1000;
#endif
}


void LinSolParams::copy (const LinSolParams& spar)
{
#ifdef HAS_PETSC
  // Copy linear solver parameters
  method    = spar.method;
  hypretype = spar.hypretype;
  prec      = spar.prec;
  package   = spar.package;
  levels    = spar.levels;
  overlap   = spar.overlap;
  nullspc   = spar.nullspc;

  atol   = spar.atol;
  rtol   = spar.rtol;
  dtol   = spar.dtol;
  maxIts = spar.maxIts;
#endif
}


bool LinSolParams::read (std::istream& is, int nparam)
{
  char* cline = 0;
  for (int i = 0; i < nparam && (cline = utl::readLine(is)); i++)
#ifdef HAS_PETSC
    if (!strncasecmp(cline,"type",4)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      method.assign(c,j);
    }

    else if (!strncasecmp(cline,"pc",2)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      prec.assign(c,j);
    }

    else if (!strncasecmp(cline,"hypretype",9)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      hypretype.assign(c,j);
    }

    else if (!strncasecmp(cline,"package",7)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      package.assign(c,j);
    }

    else if (!strncasecmp(cline,"levels",6)) {
      char* c = strchr(cline,'=');
      levels = atoi(++c);
    }

    else if (!strncasecmp(cline,"overlap",7)) {
      char* c = strchr(cline,'=');
      overlap = atoi(++c);
    }

    else if (!strncasecmp(cline,"atol",4)) {
      char* c = strchr(cline,'=');
      atol = atof(++c);
    }

    else if (!strncasecmp(cline,"rtol",4)) {
      char* c = strchr(cline,'=');
      rtol = atof(++c);
    }

    else if (!strncasecmp(cline,"dtol",4)) {
      char* c = strchr(cline,'=');
      dtol = atof(++c);
    }
    else if (!strncasecmp(cline,"maxits",6)) {
      char* c = strchr(cline,'=');
      maxIts = atoi(++c);
    }

    else if (!strncasecmp(cline,"nullspace",6)) {
      char* c = strchr(cline,'=');
      if (!strncasecmp(c,"CONSTANT",8))
	nullspc = CONSTANT;
      else if (!strncasecmp(c,"RIGID_BODY",10))
	nullspc = RIGID_BODY;
      else
	nullspc = NONE;
    }

    else {
      std::cerr <<" *** LinSolParams::read: Unknown keyword: "
		<< cline << std::endl;
      return false;
    }
#else
    ;
#endif
  return true;
}


#ifdef HAS_PETSC
static const char* IsOK (const TiXmlElement* child, const char* value)
{
  if (strcasecmp(child->Value(),value) || child->FirstChild())
    return 0;

  return child->FirstChild()->Value();
}


bool LinSolParams::read (const TiXmlElement* elem)
{
  const char* value = 0;
  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if ((value = IsOK(child,"type")))
      method = value;
    else if ((value = IsOK(child,"pc")))
      prec = value;
    else if ((value = IsOK(child,"hypretype")))
      hypretype = value;
    else if ((value = IsOK(child,"package")))
      package = value;
    else if ((value = IsOK(child,"levels")))
      levels = atoi(value);
    else if ((value = IsOK(child,"overlap")))
      overlap = atoi(value);
    else if ((value = IsOK(child,"atol")))
      atol = atof(value);
    else if ((value = IsOK(child,"rtol")))
      rtol = atof(value);
    else if ((value = IsOK(child,"dtol")))
      dtol = atof(value);
    else if ((value = IsOK(child,"maxits")))
      maxIts = atoi(value);
    else if ((value = IsOK(child,"nullspace"))) {
      if (!strcasecmp(value,"constant"))
	nullspc = CONSTANT;
      else if (!strcasecmp(value,"rigid_body"))
	nullspc = RIGID_BODY;
    }
    else
      std::cerr <<" *** LinSolParams::read: Unknown keyword: "
		<< child->Value() << std::endl;
    child = child->NextSiblingElement();
  }
  return true;
}
#else
// No need to loop trough the XML entries when PETSc is not included.
bool LinSolParams::read (const TiXmlElement*) { return true; }
#endif


bool LinSolParams::read (const char* filename)
{
  std::ifstream fs(filename);
  return this->read(fs,1000);
}


#ifdef HAS_PETSC
void LinSolParams::setParams(KSP& ksp) const
{
  // Set linear solver method
  KSPSetType(ksp,method.c_str());
  KSPSetTolerances(ksp,rtol,atol,dtol,maxIts);
  //KSPSetFromOptions(ksp);

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  //PCSetType(pc,prec.c_str());
  //PCFactorSetLevels(pc,levels);
  //PCMGSetLevels(pc,levels,PETSC_NULL);
  //PCMGSetType(pc,PC_MG_MULTIPLICATIVE);
  PCSetType(pc,prec.c_str());
#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5))
    PCHYPRESetType(pc,hypretype.c_str());
#endif
  if (overlap > 0) {
    PCASMSetType(pc,PC_ASM_BASIC);
    PCASMSetOverlap(pc,overlap);
  }
  PCFactorSetMatSolverPackage(pc,package.c_str());
  PCSetFromOptions(pc);
  PCSetUp(pc);

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
}
#endif
