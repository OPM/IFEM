// $Id$
//==============================================================================
//!
//! \file SIMFiniteDefEl.C
//!
//! \date Dec 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#include "SIMFiniteDefEl.h"
#include "LinIsotropic.h"
#include "LinearMaterial.h"
#include "NeoHookeMaterial.h"
#include "PlasticMaterial.h"
#include "PlasticityUL.h"
#include "NonlinearElasticityULMixed.h"
#include "NonlinearElasticityULMX.h"
#include "NeoHookeElasticity.h"
#include "Utilities.h"
#include "Property.h"


SIMFiniteDefEl2D::SIMFiniteDefEl2D (const std::vector<int>& options)
  : SIMLinEl2D(SIM::NONE)
{
  int form = options.size() > 0 ? options[0] : 0; // problem formulation
  int pOrd = options.size() > 1 ? options[1] : 0; // pressure field order

  switch (form)
    {
    case SIM::PLASTICITY:
      myProblem = new PlasticityUL(2);
      break;
    case SIM::MIXED_QnQn1:
      nf[1] = 2; // continuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMixed(2);
      break;
    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMX(2,pOrd);
      break;
    case SIM::UPDATED_LAGRANGE:
      myProblem = new NonlinearElasticityUL(2);
      break;
    case SIM::TOTAL_LAGRANGE:
      myProblem = new NonlinearElasticityTL(2);
      break;
    case SIM::NONLINEAR: // Old tensor-based TL-formulation
      myProblem = new NonlinearElasticity(2);
      break;
    default:
      std::cerr <<" *** SIMFiniteDefEl2D: Unknown problem formulation "
		<< form << std::endl;
    }
}


bool SIMFiniteDefEl2D::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+9);
    if (myPid == 0)
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    char* cline = 0;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = atof(strtok(NULL," "));
      int matVer = (cline = strtok(NULL," ")) ? atoi(cline) : -1;
      if (matVer >= 0)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinearMaterial(new LinIsotropic(E,nu,rho)));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" ("<< matVer <<")"<< std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"PLASTIC",7))
  {
    int nmat = atoi(keyWord+7);
    if (myPid == 0)
      std::cout <<"\nNumber of plastic materials: "<< nmat << std::endl;

    char* cline = 0;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      RealArray pMAT;
      double rho = atof(strtok(NULL," "));
      while ((cline = strtok(NULL," ")))
      {
	pMAT.push_back(rho);
	rho = atof(cline);
      }
      mVec.push_back(new PlasticPrm(pMAT,rho));
      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout <<" rho = "<< rho << std::endl;
      }
    }
  }

  else
    return this->SIMLinEl2D::parse(keyWord,is);

  if (!mVec.empty())
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
    if (elp) elp->setMaterial(mVec.front());
  }

  return true;
}


SIMFiniteDefEl3D::SIMFiniteDefEl3D (bool checkRHS,
				    const std::vector<int>& options)
  : SIMLinEl3D(checkRHS,SIM::NONE)
{
  int form = options.size() > 0 ? options[0] : 0; // problem formulation
  int pOrd = options.size() > 1 ? options[1] : 0; // pressure field order

  switch (form)
    {
    case SIM::PLASTICITY:
      myProblem = new PlasticityUL(2);
      break;
    case SIM::MIXED_QnQn1:
      nf[1] = 2; // continuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMixed();
      break;
    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMX(3,pOrd);
      break;
    case SIM::UPDATED_LAGRANGE:
      myProblem = new NonlinearElasticityUL();
      break;
    case SIM::TOTAL_LAGRANGE:
      myProblem = new NonlinearElasticityTL();
      break;
    case SIM::NONLINEAR: // Old tensor-based TL-formulation
      myProblem = new NonlinearElasticity();
      break;
    case SIM::NEOHOOKE: // NeoHookean TL-formulation (not working)
      myProblem = new NeoHookeElasticity();
      break;
    case SIM::NEOHOOKE_IV: // NeoHookean isochoric/volumetric TL-formulation
      myProblem = new NeoHookeElasticityIV();
      break;
    default:
      std::cerr <<" *** SIMFiniteDefEl3D: Unknown problem formulation "
		<< form << std::endl;
    }
}


bool SIMFiniteDefEl3D::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+9);
    if (myPid == 0)
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    char* cline = 0;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = atof(strtok(NULL," "));
      int matVer = (cline = strtok(NULL," ")) ? atoi(cline) : -1;
      if (matVer >= 0)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinearMaterial(new LinIsotropic(E,nu,rho)));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" ("<< matVer <<")"<< std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"PLASTIC",7))
  {
    int nmat = atoi(keyWord+7);
    if (myPid == 0)
      std::cout <<"\nNumber of plastic materials: "<< nmat << std::endl;

    char* cline = 0;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      RealArray pMAT;
      double rho = atof(strtok(NULL," "));
      while ((cline = strtok(NULL," ")))
      {
	pMAT.push_back(rho);
	rho = atof(cline);
      }
      mVec.push_back(new PlasticPrm(pMAT,rho));
      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout <<" rho = "<< rho << std::endl;
      }
    }
  }

  else
    return this->SIMLinEl3D::parse(keyWord,is);

  if (!mVec.empty())
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
    if (elp) elp->setMaterial(mVec.front());
  }

  return true;
}
