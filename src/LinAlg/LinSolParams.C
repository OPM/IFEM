// $Id: LinSolParams.C,v 1.5 2010-12-06 09:08:28 rho Exp $
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
#include "petscmg.h"
#endif


void LinSolParams::setDefault ()
{
#ifdef HAS_PETSC
  // Use GMRES with ILU preconditioner as default
  method  = KSPGMRES;
  prec    = PCILU;
  package = MAT_SOLVER_PETSC; 
  levels  = 0;
  overlap = 0;

  atol   = 1.0e-6;
  rtol   = 1.0e-6;
  dtol   = 1.0e-6;
  maxIts = 1000;
#endif
}


void LinSolParams::copy (const LinSolParams& spar)
{
#ifdef HAS_PETSC
  // Copy linear solver parameters
  method  = spar.method;
  prec    = spar.prec;
  package = spar.package;
  levels  = spar.levels;
  overlap = spar.overlap;

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

    else
    {
      std::cerr <<" *** LinSolParams::read: Unknown keyword: "
		<< cline << std::endl;
      return false;
    }
#else
    ;
#endif
  return true;
}


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
  if (overlap > 0)
    PCASMSetOverlap(pc,overlap);
  PCFactorSetMatSolverPackage(pc,package.c_str());
  PCSetFromOptions(pc);
  PCSetUp(pc);

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
}
#endif
