// $Id$
//==============================================================================
//!
//! \file SIMLinElKL.C
//!
//! \date Sep 16 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love plates.
//!
//==============================================================================

#include "SIMLinElKL.h"
#include "LinIsotropic.h"
#include "KirchhoffLovePlate.h"
#include "AnalyticSolutions.h"
#include "Utilities.h"
#include "Property.h"
#include "AnaSol.h"
#include <string.h>


SIMLinElKL::SIMLinElKL () : SIMLinEl2D(SIM::NONE)
{
  nf[0] = 1;
  myProblem = new KirchhoffLovePlate();
}


bool SIMLinElKL::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    double g = atof(strtok(keyWord+7," "));
    if (klp)
      klp->setGravity(g);
    if (myPid == 0)
      std::cout <<"\nGravitation constant: "<< g << std::endl;
  }

  else if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+10);
    if (myPid == 0)
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
	this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      double thk = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      mVec.push_back(new LinIsotropic(E,nu,rho,true));
      tVec.push_back(thk);
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho <<" "<< thk << std::endl;
    }

    if (!mVec.empty())
      klp->setMaterial(mVec.front());
    if (!tVec.empty() && tVec.front() != 0.0)
      klp->setThickness(tVec.front());
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"NAVIERPLATE",7))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double t  = atof(strtok(NULL," "));
      double pz = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: NavierPlate a="<< a <<" b="<< b
		<<" t="<< t <<" pz="<< pz <<" E="<< E <<" nu="<< nu
		<< std::endl;
      mySol = new AnaSol(new NavierPlate(a,b,t,pz,E,nu));
    }
    else
    {
      std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
		<< cline <<" (ignored)"<< std::endl;
      return true;
    }
  }

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMLinElKL::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  klp->setMaterial(mVec[propInd]);
  if (tVec[propInd] != 0.0)
    klp->setThickness(tVec[propInd]);

  return true;
}
