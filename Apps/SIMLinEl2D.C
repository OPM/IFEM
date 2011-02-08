// $Id: SIMLinEl2D.C,v 1.18 2011-02-08 09:06:02 kmo Exp $
//==============================================================================
//!
//! \file SIMLinEl2D.C
//!
//! \date Dec 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 2D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMLinEl2D.h"
#include "LinearElasticity.h"
#include "NonlinearElasticity.h"
#include "AnalyticSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Property.h"
#include "Tensor.h"
#include <string.h>


SIMLinEl2D::SIMLinEl2D (int form, bool planeStress)
{
  switch (form)
    {
    case SIM::LINEAR:
      myProblem = new LinearElasticity(2,planeStress);
      break;
    case SIM::NONLINEAR:
      myProblem = new NonlinearElasticity(2,planeStress);
      break;
    }

  asol = 0;
}


SIMLinEl2D::~SIMLinEl2D ()
{
  if (asol) delete asol;
}


bool SIMLinEl2D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    double gx = atof(strtok(keyWord+7," "));
    double gy = atof(strtok(NULL," "));
    elp->setGravity(gx,gy,0.0);
    if (myPid == 0)
      std::cout <<"\nGravitation vector: "
		<< gx <<" "<< gy << std::endl;
  }

  else if (!strncasecmp(keyWord,"ISOTROPHIC",10))
  {
    int nmat = atoi(keyWord+10);
    if (myPid == 0)
      std::cout <<"\nNumber of isotrophic materials: "<< nmat << std::endl;

    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int   code = atoi(strtok(cline," "));
      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = atof(strtok(NULL," "));
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho << std::endl;
      if (code == 0 || i == 0)
	elp->setMaterial(E,nu,rho);
      for (unsigned int j = 0; j < myProps.size() && code > 0; j++)
	if (myProps[j].pindx == (size_t)code &&
	    myProps[j].pcode == Property::UNDEFINED)
	{
	  myProps[j].pindx = mVec.size();
	  myProps[j].pcode = Property::MATERIAL;
	}
      if (code > 0)
	mVec.push_back(IsoMat(E,nu,rho));
    }
  }

  else if (!strncasecmp(keyWord,"CONSTANT_PRESSURE",17))
  {
    int npres = atoi(keyWord+17);
    if (myPid == 0)
      std::cout <<"\nNumber of pressures: "<< npres << std::endl;

    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      int pdir = atoi(strtok(NULL," "));
      double p = atof(strtok(NULL," "));
      if (myPid == 0)
	std::cout <<"\tPressure code "<< code <<" direction "<< pdir
		  <<": "<< p << std::endl;

      this->setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new PressureField(p,pdir);
    }
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"HOLE",4))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      asol = new Hole(a,F0,nu);
      std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      asol = new Lshape(a,F0,nu);
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      asol = new CanTS(L,H,F0);
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		<<" F0="<< F0 << std::endl;
    }
    else if (!strncasecmp(cline,"CANTM",5))
    {
      double H  = atof(strtok(NULL," "));
      double M0 = atof(strtok(NULL," "));
      asol = new CanTM(H,M0);
      std::cout <<"\nAnalytical solution: CanTM H="<< H
		<<" M0="<< M0 << std::endl;
    }
    else if (!strncasecmp(cline,"CURVED",6))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double u0 = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      asol = new CurvedBeam(u0,a,b,E);
      std::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
                <<" u0="<< u0 <<" E="<< E << std::endl;
    }
    else
      std::cerr <<"  ** SIMLinEl2D::parse: Unknown analytical solution "
		<< cline << std::endl;

    // Define the analytical boundary traction field
    int code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && asol)
    {
      std::cout <<"Pressure code "<< code <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*asol);
    }
  }

  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of properties onto
  // the topological entities (blocks and faces) of the model.

  else if (!strncasecmp(keyWord,"PRESSURE",8))
  {
    Property press;
    press.pcode = Property::NEUMANN;
    press.ldim = 1;

    int npres = atoi(keyWord+8);
    std::cout <<"\nNumber of pressures: "<< npres << std::endl;
    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      press.pindx = 1+i;
      press.patch = atoi(strtok(cline," "));

      int pid = this->getLocalPatchIndex(press.patch);
      if (pid < 0) return false;
      if (pid < 1) continue;

      press.lindx = atoi(strtok(NULL," "));
      if (press.lindx < 1 || press.lindx > 4)
      {
	std::cerr <<" *** SIMLinEl2D::parse: Invalid edge index "
		  << (int)press.lindx << std::endl;
	return false;
      }

      if (asol)
      {
	std::cout <<"\tTraction on P"<< press.patch
		  <<" E"<< (int)press.lindx << std::endl;
	myTracs[1+i] = new TractionField(*asol);
      }
      else
      {
	int pdir = atoi(strtok(NULL," "));
	double p = atof(strtok(NULL," "));
	std::cout <<"\tPressure on P"<< press.patch
		  <<" E"<< (int)press.lindx <<" direction "<< pdir <<": ";
	if ((cline = strtok(NULL," ")))
	  myTracs[1+i] = new PressureField(utl::parseRealFunc(cline,p),pdir);
	else
	{
	  std::cout << p;
	  myTracs[1+i] = new PressureField(p,pdir);
	}
	std::cout << std::endl;
      }

      press.patch = pid;
      myProps.push_back(press);
    }
  }

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    int nmat = atoi(keyWord+8);
    std::cout <<"\nNumber of materials: "<< nmat << std::endl;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      double E   = atof(strtok(cline," "));
      double nu  = atof(strtok(NULL," "));
      double rho = atof(strtok(NULL," "));
      while ((cline = strtok(NULL," ")))
	if (!strncasecmp(cline,"ALL",3))
        {
	  std::cout <<"\tMaterial for all patches: "
		    << E <<" "<< nu <<" "<< rho << std::endl;
	  elp->setMaterial(E,nu,rho);
	}
	else
        {
	  int patch = atoi(cline);
	  int pid = this->getLocalPatchIndex(patch);
	  if (pid < 0) return false;
	  if (pid < 1) continue;

	  std::cout <<"\tMaterial for P"<< patch
		    <<": "<< E <<" "<< nu <<" "<< rho << std::endl;
	  myProps.push_back(Property(Property::MATERIAL,mVec.size(),pid,3));
	  mVec.push_back(IsoMat(E,nu,rho));
	  if (i == 0)
	    elp->setMaterial(E,nu,rho);
	}
    }
  }

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMLinEl2D::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd].E,mVec[propInd].nu,mVec[propInd].rho);
  return true;
}


bool SIMLinEl2D::initNeumann (size_t propInd)
{
  TracFuncMap::const_iterator tit = myTracs.find(propInd);
  if (tit == myTracs.end()) return false;

  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  elp->setTraction(tit->second);
  return true;
}
