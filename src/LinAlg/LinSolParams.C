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
#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <iterator>
#ifdef HAS_PETSC
#include "petscversion.h"
#include "petscpcmg.h"
#endif
#include "tinyxml.h"


void LinSolParams::setDefault ()
{
#ifdef HAS_PETSC
  // Use GMRES with ILU preconditioner as default
  method      = KSPGMRES;
  prec        = PCILU;
  subprec.resize(1,PCILU);
  hypretype.resize(1,"boomeramg");
  package.resize(1,MATSOLVERPETSC);
  levels.resize(1,0);
  mglevels.resize(1,1);
  presmoother.resize(1,PCILU);
  noPreSmooth.resize(1,1);
  postsmoother.resize(1,PCILU);
  noPostSmooth.resize(1,1);
  maxCoarseSize.resize(1,1000);
  overlap.resize(1,1);
  nx.resize(1,0);
  ny.resize(1,0);
  nz.resize(1,0);
  nullspc.resize(1,NONE);
  asmlu.resize(1,false);
  nblock   = 1;
  schur    = false;
  schurPrec = SIMPLE;
  ncomps.resize(2); ncomps[0]   = 2; ncomps[1]   = 1;
    
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
  method        = spar.method;
  prec          = spar.prec;
  subprec       = spar.subprec;
  hypretype     = spar.hypretype;
  package       = spar.package;
  overlap       = spar.overlap;
  levels        = spar.levels;
  mglevels      = spar.mglevels;
  presmoother   = spar.presmoother;
  postsmoother  = spar.postsmoother;
  noPreSmooth   = spar.noPreSmooth;
  noPostSmooth  = spar.noPostSmooth;
  maxCoarseSize = spar.maxCoarseSize;
  overlap       = spar.overlap;
  nx            = spar.nx;
  ny            = spar.ny;
  nz            = spar.nz;
  nullspc       = spar.nullspc;
  asmlu         = spar.asmlu;
  nblock        = spar.nblock;
  schur         = spar.schur;
  schurPrec     = spar.schurPrec;
  ncomps        = spar.ncomps;

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
      if ((asmlu[0] = (prec.substr(0,5) == "asmlu")))
	prec = "asm";
    }

    else if (!strncasecmp(cline,"hypretype",9)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      hypretype[0].assign(c,j);
    }

    else if (!strncasecmp(cline,"package",7)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      package[0].assign(c,j);
    }

    else if (!strncasecmp(cline,"levels",6)) {
      char* c = strchr(cline,'=');
      levels[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"overlap",7)) {
      char* c = strchr(cline,'=');
      overlap[0] = atoi(++c);
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
      nx[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"ny",2)) {
      char* c = strchr(cline,'=');
      ny[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"nz",2)) {
      char* c = strchr(cline,'=');
      nz[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"ncomponents",11)) {
      char* c = strchr(cline,'=');
      std::istringstream this_line(c);
      std::istream_iterator<int> begin(this_line), end;
      ncomps.assign(begin, end);
      nblock = ncomps.size();
    }

    else if (!strncasecmp(cline,"nullspace",6)) {
      char* c = strchr(cline,'=');
      if (!strncasecmp(c,"CONSTANT",8))
	nullspc[0] = CONSTANT;
      else if (!strncasecmp(c,"RIGID_BODY",10))
	nullspc[0] = RIGID_BODY;
      else
	nullspc[0] = NONE;
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
  else if ((value = utl::getValue(child,"pc"))) {
    prec = value;
    if (!strncasecmp(prec.c_str(),"asmlu",5)) {
      prec = "asm";
      asmlu[0] = true;
    }
  }
  else if ((value = utl::getValue(child,"schurpc"))) {
    if (!strncasecmp(value,"SIMPLE",6))
      schurPrec = SIMPLE;
    else if (!strncasecmp(value,"MSIMPLER",8))
      schurPrec = MSIMPLER;
    else if (!strncasecmp(value,"PCD",3))
      schurPrec = PCD;
  }
  else if ((value = utl::getValue(child,"subpc"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    subprec.assign(begin,end);
    asmlu.resize(subprec.size());
    for (size_t i = 0;i < subprec.size();i++)
      if (!strncasecmp(subprec[i].c_str(),"asmlu",5)) {
	asmlu[i] = true;
	subprec[i] = "asm";
      }
  }
  else if ((value = utl::getValue(child,"hypretype"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    hypretype.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"package"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    package.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"levels"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    levels.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"mglevels"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    mglevels.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"noPreSmooth"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    noPreSmooth.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"noPostSmooth"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    noPostSmooth.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"presmoother"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    presmoother.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"postsmoother"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    postsmoother.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"overlap"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    overlap.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"nx"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    nx.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"ny"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    ny.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"nz"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    nz.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"atol")))
    atol = atof(value);
  else if ((value = utl::getValue(child,"rtol")))
    rtol = atof(value);
  else if ((value = utl::getValue(child,"dtol")))
    dtol = atof(value);
  else if ((value = utl::getValue(child,"maxits")))
    maxIts = atoi(value);
  else if ((value = utl::getValue(child,"maxCoarseSize"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    maxCoarseSize.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"ncomponents"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    ncomps.assign(begin, end);
    nblock = ncomps.size();
  }
  else if ((value = utl::getValue(child,"nullspace"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    std::istream_iterator<std::string> it;
    for (it = begin;it != end;it++) {
      if (!strcasecmp(it->c_str(),"constant"))
	nullspc.push_back(CONSTANT);
      else if (!strcasecmp(value,"rigid_body"))
	nullspc.push_back(RIGID_BODY);
    }
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
  //PCFactorSetLevels(pc,levels[0]);
  //PCMGSetLevels(pc,levels[0],PETSC_NULL);
  //PCMGSetType(pc,PC_MG_MULTIPLICATIVE);
  PCSetType(pc,prec.c_str());
#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5))
    PCHYPRESetType(pc,hypretype[0].c_str());
#endif
  if (!strncasecmp(prec.c_str(),"asm",3) ||!strncasecmp(prec.c_str(),"gasm",4)) {
    PCASMSetType(pc,PC_ASM_BASIC);
    PCASMSetOverlap(pc,overlap[0]);
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

  PCFactorSetMatSolverPackage(pc,package[0].c_str());
  // RUNAR
  //PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);
  PCSetFromOptions(pc);
  PCSetUp(pc);

  // If LU factorization is used on each subdomain
  if (asmlu[0]) {
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
