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
  overlap   = 1;
  npart[0] = npart[1] = npart[2] = 0;
  nullspc   = NONE;
  asmlu     = false;

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
  npart[0]  = spar.npart[0];
  npart[1]  = spar.npart[1];
  npart[2]  = spar.npart[2];               
  nullspc   = spar.nullspc;
  asmlu     = spar.asmlu;

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
      if ((asmlu = (prec.substr(0,5) == "asmlu")))
	prec = "asm";
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

    else if (!strncasecmp(cline,"nx",2)) {
      char* c = strchr(cline,'=');
      npart[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"ny",2)) {
      char* c = strchr(cline,'=');
      npart[1] = atoi(++c);
    }

    else if (!strncasecmp(cline,"nz",2)) {
      char* c = strchr(cline,'=');
      npart[2] = atoi(++c);
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


bool LinSolParams::read (const TiXmlElement* child)
{
#ifdef HAS_PETSC
  const char* value = 0;
  if ((value = utl::getValue(child,"type")))
    method = value;
  else if ((value = utl::getValue(child,"pc")))
  {
    asmlu = !strncasecmp(value,"asmlu",5);
    prec  = asmlu ? "asm" : value;
  }
  else if ((value = utl::getValue(child,"hypretype")))
    hypretype = value;
  else if ((value = utl::getValue(child,"package")))
    package = value;
  else if ((value = utl::getValue(child,"levels")))
    levels = atoi(value);
  else if ((value = utl::getValue(child,"overlap")))
    overlap = atoi(value);
  else if ((value = utl::getValue(child,"atol")))
    atol = atof(value);
  else if ((value = utl::getValue(child,"rtol")))
    rtol = atof(value);
  else if ((value = utl::getValue(child,"dtol")))
    dtol = atof(value);
  else if ((value = utl::getValue(child,"maxits")))
    maxIts = atoi(value);
  else if ((value = utl::getValue(child,"nx")))
    npart[0] = atoi(value);
  else if ((value = utl::getValue(child,"ny")))
    npart[1] = atoi(value);
  else if ((value = utl::getValue(child,"nz")))
    npart[2] = atoi(value);
  else if ((value = utl::getValue(child,"nullspace"))) {
    if (!strcasecmp(value,"constant"))
      nullspc = CONSTANT;
    else if (!strcasecmp(value,"rigid_body"))
      nullspc = RIGID_BODY;
  }
  else
  {
    std::cerr <<" *** LinSolParams::read: Unknown keyword: "
	      << child->Value() << std::endl;
    return false;
  }
#endif
  return true;
}


bool LinSolParams::read (const char* filename)
{
  std::ifstream fs(filename);
  return this->read(fs,1000);
}


#ifdef HAS_PETSC
void LinSolParams::setParams (KSP& ksp, std::vector<std::vector<PetscInt> >& locSubdDofs,
			      std::vector<std::vector<PetscInt> >& subdDofs) const
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
  if (!strncasecmp(prec.c_str(),"asm",3) ||!strncasecmp(prec.c_str(),"gasm",4)) {
    PCASMSetType(pc,PC_ASM_BASIC);
    PCASMSetOverlap(pc,overlap);
  }

  if (!locSubdDofs.empty() && !subdDofs.empty()) {
    const size_t nsubds = subdDofs.size();

    IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
    for (size_t i = 0;i < nsubds;i++) {
      ISCreateGeneral(PETSC_COMM_WORLD,locSubdDofs[i].size(),&(locSubdDofs[i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
      ISCreateGeneral(PETSC_COMM_WORLD,subdDofs[i].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
    }
    PCASMSetLocalSubdomains(pc,nsubds,isSubdDofs,isLocSubdDofs);
  }

  PCFactorSetMatSolverPackage(pc,package.c_str());
  PCSetFromOptions(pc);
  PCSetUp(pc);

  // If LU factorization is used on each subdomain
  if (asmlu) {
    KSP* subksp;
    PC   subpc;
    PetscInt first, nlocal;
    PCASMGetSubKSP(pc,&nlocal,&first,&subksp);

    for (int i = 0; i < nlocal; i++) {
      KSPGetPC(subksp[i],&subpc);
      PCSetType(subpc,PCLU);
      KSPSetType(subksp[i],KSPPREONLY);
    }
  }

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
}
#endif
