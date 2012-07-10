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
#include "Elasticity.h"
#include "AlgEqSystem.h"
#include "Utilities.h"
#include "Property.h"
#include "tinyxml.h"


SIMFiniteDefEl2D::SIMFiniteDefEl2D (const std::vector<int>& options) : nlo(2)
{
  nlo.form = options.size() > 0 ? options[0] : 0; // problem formulation
  nlo.pOrd = options.size() > 1 ? options[1] : 0; // pressure field order
  if (nlo.form == SIM::MIXED_QnQn1) nf[1] = 2;
}


SIMFiniteDefEl2D::~SIMFiniteDefEl2D ()
{
  for (IntegrandMap::const_iterator i = myInts.begin(); i != myInts.end(); i++)
    if (i->second != myProblem)
      delete i->second;
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
      if (matVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (matVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.back() = new LinearMaterial(mVec.back());

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
      while ((cline = strtok(NULL," ")))
	pMAT.push_back(atof(cline));
      mVec.push_back(new PlasticMaterial(pMAT));

      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout << std::endl;
      }
    }
  }

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    std::cerr <<" *** SIMFiniteDefEl2D::parse: The keyword MATERIAL"
	      <<" is not supported\n     for finite deformation analysis."
	      <<" You must use ISOTROPIC instead."<< std::endl;
    return false;
  }
  else
    return this->SIMLinEl2D::parse(keyWord,is);

  if (!myProblem)
    myProblem = nlo.getIntegrand();

  if (!mVec.empty())
    static_cast<Elasticity*>(myProblem)->setMaterial(mVec.front());

  return true;
}


bool SIMFiniteDefEl2D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"finitedeformation"))
    return this->SIMLinEl2D::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"formulation"))
    {
      nlo.parse(child);
      nf[1] = nlo.form == SIM::MIXED_QnQn1 ? 2 : 0;
    }

    else if (!strcasecmp(child->Value(),"isotropic"))
    {
      int code = this->parseMaterialSet(child,mVec.size());

      int matVer = -1;
      double E = 1000.0, nu = 0.3, rho = 1.0;
      utl::getAttribute(child,"K",E);
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"G",nu);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      utl::getAttribute(child,"version",matVer);
      if (matVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (matVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.back() = new LinearMaterial(mVec.back());

      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" ("<< matVer <<")"<< std::endl;
    }

    else if (!strcasecmp(child->Value(),"plastic"))
    {
      int code = this->parseMaterialSet(child,mVec.size());

      RealArray pMAT;
      if (child->FirstChild())
      {
	std::string value(child->FirstChild()->Value());
	char* cval = strtok(const_cast<char*>(value.c_str())," ");
	for (; cval; cval = strtok(NULL," "))
	  pMAT.push_back(atof(cval));
      }
      if (pMAT.size() < 11) pMAT.resize(11,0.0);
      utl::getAttribute(child,"Bmod" ,pMAT[0]);
      utl::getAttribute(child,"Emod" ,pMAT[0]);
      utl::getAttribute(child,"Smod" ,pMAT[1]);
      utl::getAttribute(child,"nu"   ,pMAT[1]);
      utl::getAttribute(child,"rho"  ,pMAT[3]);
      utl::getAttribute(child,"Hiso" ,pMAT[4]);
      utl::getAttribute(child,"Hkin" ,pMAT[5]);
      utl::getAttribute(child,"yield",pMAT[7]);
      utl::getAttribute(child,"Y0"   ,pMAT[8]);
      utl::getAttribute(child,"Yinf" ,pMAT[8]);
      utl::getAttribute(child,"beta" ,pMAT[9]);
      utl::getAttribute(child,"istrt",pMAT[10]);
      mVec.push_back(new PlasticMaterial(pMAT));

      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout << std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"contact"))
    {
      if (!this->createFEMmodel())
	return false;

      if (!this->parseContactTag(child,myModel,myEntitys))
	return false;

      preserveNOrder = true; // because extra nodes have been added
    }

  if (!myProblem)
    myProblem = nlo.getIntegrand();

  if (!mVec.empty())
    static_cast<Elasticity*>(myProblem)->setMaterial(mVec.front());

  return true;
}


bool SIMFiniteDefEl2D::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!this->SIMLinEl2D::preprocess(ignored,fixDup))
    return false;

  myInts.insert(std::make_pair(0,myProblem));
  if (this->withContact())
  {
    opt.num_threads_SLU *= -1; // do not lock the sparsity pattern
    return this->preprocessContact(myInts,*mySam,this->getNoSpaceDim());
  }

  return true;
}


bool SIMFiniteDefEl2D::createContactSet (const std::string& slaveSet, int& code)
{
  if (!(code = this->getUniquePropertyCode(slaveSet)))
    return false;

  this->setPropertyType(code,Property::NEUMANN_GENERIC);
  return true;
}


void SIMFiniteDefEl2D::preprocessBeforeAsmInit (int& ngnod)
{
  this->renumberContactBodies(*g2l);
  for (PatchVec::const_iterator it = myModel.begin(); it != myModel.end(); it++)
    this->addLagrangeMultipliers(*it,ngnod);
}


bool SIMFiniteDefEl2D::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    this->SIMContact::updateDirichlet(time);

  return this->SIMbase::updateDirichlet(time,prevSol);
}


bool SIMFiniteDefEl2D::updateConfiguration (const Vector& solution)
{
  return this->updateContactBodies(solution);
}


bool SIMFiniteDefEl2D::assembleDiscreteTerms (const IntegrandBase* problem)
{
  return this->assembleMortarTangent(problem,
				     myEqSys->getMatrix(),
				     myEqSys->getVector());
}


SIMFiniteDefEl3D::SIMFiniteDefEl3D (bool checkRHS,
				    const std::vector<int>& options)
  : SIMLinEl3D(checkRHS), nlo(3)
{
  nlo.form = options.size() > 0 ? options[0] : 0; // problem formulation
  nlo.pOrd = options.size() > 1 ? options[1] : 0; // pressure field order
  if (nlo.form == SIM::MIXED_QnQn1) nf[1] = 2;
}


SIMFiniteDefEl3D::~SIMFiniteDefEl3D ()
{
  for (IntegrandMap::const_iterator i = myInts.begin(); i != myInts.end(); i++)
    if (i->second != myProblem)
      delete i->second;
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
      if (matVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinIsotropic(E,nu,rho));
      if (matVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.back() = new LinearMaterial(mVec.back());

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
      while ((cline = strtok(NULL," ")))
	pMAT.push_back(atof(cline));
      mVec.push_back(new PlasticMaterial(pMAT));

      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout << std::endl;
      }
    }
  }

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    std::cerr <<" *** SIMFiniteDefEl3D::parse: The MATERIAL keyword"
	      <<" is not supported\n     for finite deformation analysis."
	      <<" You must use ISOTROPIC instead."<< std::endl;
    return false;
  }
  else
    return this->SIMLinEl3D::parse(keyWord,is);

  if (!myProblem)
    myProblem = nlo.getIntegrand();

  if (!mVec.empty())
    static_cast<Elasticity*>(myProblem)->setMaterial(mVec.front());

  return true;
}


bool SIMFiniteDefEl3D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"finitedeformation"))
    return this->SIMLinEl3D::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"formulation"))
    {
      nlo.parse(child);
      nf[1] = nlo.form == SIM::MIXED_QnQn1 ? 2 : 0;
    }

    else if (!strcasecmp(child->Value(),"isotropic"))
    {
      int code = this->parseMaterialSet(child,mVec.size());

      int matVer = -1;
      double E = 1000.0, nu = 0.3, rho = 1.0;
      utl::getAttribute(child,"K",E);
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"G",nu);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      utl::getAttribute(child,"version",matVer);
      if (matVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else
	mVec.push_back(new LinIsotropic(E,nu,rho));
      if (matVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
	mVec.back() = new LinearMaterial(mVec.back());

      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" ("<< matVer <<")"<< std::endl;
    }

    else if (!strcasecmp(child->Value(),"plastic"))
    {
      int code = this->parseMaterialSet(child,mVec.size());

      RealArray pMAT;
      if (child->FirstChild())
      {
	std::string value(child->FirstChild()->Value());
	char* cval = strtok(const_cast<char*>(value.c_str())," ");
	for (; cval; cval = strtok(NULL," "))
	  pMAT.push_back(atof(cval));
      }
      if (pMAT.size() < 11) pMAT.resize(11,0.0);
      utl::getAttribute(child,"Bmod" ,pMAT[0]);
      utl::getAttribute(child,"Emod" ,pMAT[0]);
      utl::getAttribute(child,"Smod" ,pMAT[1]);
      utl::getAttribute(child,"nu"   ,pMAT[1]);
      utl::getAttribute(child,"rho"  ,pMAT[3]);
      utl::getAttribute(child,"Hiso" ,pMAT[4]);
      utl::getAttribute(child,"Hkin" ,pMAT[5]);
      utl::getAttribute(child,"yield",pMAT[7]);
      utl::getAttribute(child,"Y0"   ,pMAT[8]);
      utl::getAttribute(child,"Yinf" ,pMAT[8]);
      utl::getAttribute(child,"beta" ,pMAT[9]);
      utl::getAttribute(child,"istrt",pMAT[10]);
      mVec.push_back(new PlasticMaterial(pMAT));

      if (myPid == 0)
      {
        std::cout <<"\tMaterial code "<< code <<":";
	for (size_t i = 0; i < pMAT.size(); i++)
	  std::cout <<" "<< pMAT[i];
	std::cout << std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"contact"))
    {
      if (!this->createFEMmodel())
	return false;

      if (!this->parseContactTag(child,myModel,myEntitys))
	return false;

      preserveNOrder = true; // because extra nodes have been added
    }

  if (!myProblem)
    myProblem = nlo.getIntegrand();

  if (!mVec.empty())
    static_cast<Elasticity*>(myProblem)->setMaterial(mVec.front());

  return true;
}


bool SIMFiniteDefEl3D::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!this->SIMLinEl3D::preprocess(ignored,fixDup))
    return false;

  myInts.insert(std::make_pair(0,myProblem));
  if (this->withContact())
  {
    opt.num_threads_SLU *= -1; // do not lock the sparsity pattern
    return this->preprocessContact(myInts,*mySam,this->getNoSpaceDim());
  }

  return true;
}


bool SIMFiniteDefEl3D::createContactSet (const std::string& slaveSet, int& code)
{
  if (!(code = this->getUniquePropertyCode(slaveSet)))
    return false;

  this->setPropertyType(code,Property::NEUMANN_GENERIC);
  return true;
}


void SIMFiniteDefEl3D::preprocessBeforeAsmInit (int& ngnod)
{
  this->renumberContactBodies(*g2l);
  for (PatchVec::const_iterator it = myModel.begin(); it != myModel.end(); it++)
    this->addLagrangeMultipliers(*it,ngnod);
}


bool SIMFiniteDefEl3D::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    this->SIMContact::updateDirichlet(time);

  return this->SIMbase::updateDirichlet(time,prevSol);
}


bool SIMFiniteDefEl3D::updateConfiguration (const Vector& solution)
{
  return this->updateContactBodies(solution);
}


bool SIMFiniteDefEl3D::assembleDiscreteTerms (const IntegrandBase* problem)
{
  return this->assembleMortarTangent(problem,
				     myEqSys->getMatrix(),
				     myEqSys->getVector());
}
