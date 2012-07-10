// $Id$
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
#include "LinIsotropic.h"
#include "LinearElasticity.h"
#include "AnalyticSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "Property.h"
#include "AnaSol.h"
#include "tinyxml.h"


bool SIMLinEl2D::planeStrain = false;
bool SIMLinEl2D::axiSymmetry = false;


SIMLinEl2D::~SIMLinEl2D ()
{
  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode > 0) myVectors.erase(aCode);

  for (size_t i = 0; i < mVec.size(); i++)
    delete mVec[i];
}


void SIMLinEl2D::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode > 0) myVectors.erase(aCode);
  aCode = 0;

  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (elp)
  {
    elp->setMaterial(NULL);
    elp->setBodyForce(NULL);
    elp->setTraction((VecFunc*)NULL);
    elp->setTraction((TractionFunc*)NULL);
  }

  for (size_t i = 0; i < mVec.size(); i++)
    delete mVec[i];
  mVec.clear();

  this->SIMbase::clearProperties();
}


bool SIMLinEl2D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  int nmat = 0;
  double gx = 0.0, gy = 0.0;

  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    gx = atof(strtok(keyWord+7," "));
    gy = atof(strtok(NULL," "));
    if (myPid == 0)
      std::cout <<"\nGravitation vector: "
		<< gx <<" "<< gy << std::endl;
  }

  else if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    nmat = atoi(keyWord+10);
    if (myPid == 0)
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
	this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = atof(strtok(NULL," "));
      mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho << std::endl;
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
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"HOLE",4))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Hole(a,F0,nu));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu));
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		<<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0));
    }
    else if (!strncasecmp(cline,"CANTM",5))
    {
      double H  = atof(strtok(NULL," "));
      double M0 = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: CanTM H="<< H
		<<" M0="<< M0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTM(H,M0));
    }
    else if (!strncasecmp(cline,"CURVED",6))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double u0 = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
		<<" u0="<< u0 <<" E="<< E << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinEl2D::parse: Invalid analytical solution "
		<< cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getStressSol())
    {
      std::cout <<"Pressure code "<< code <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*mySol->getStressSol());
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

      if (mySol && mySol->getStressSol())
      {
	std::cout <<"\tTraction on P"<< press.patch
		  <<" E"<< (int)press.lindx << std::endl;
	myTracs[1+i] = new TractionField(*mySol->getStressSol());
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
    nmat = atoi(keyWord+8);
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
	  mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
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
	  mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
	}
    }
  }

  else
    return this->SIM2D::parse(keyWord,is);

  if (gx != 0.0 || gy != 0.0)
    this->getIntegrand()->setGravity(gx,gy);
  if (nmat > 0 && !mVec.empty())
    this->getIntegrand()->setMaterial(mVec.front());

  return true;
}


bool SIMLinEl2D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIM2D::parse(elem);

  double gx = 0.0, gy = 0.0;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity")) {
      utl::getAttribute(child,"x",gx);
      utl::getAttribute(child,"y",gy);
      if (myPid == 0)
        std::cout <<"\tGravitation vector: "<< gx <<" "<< gy << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      int code = this->parseMaterialSet(child,mVec.size());

      double E = 1000.0, nu = 0.3, rho = 1.0;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      const TiXmlElement* opt = child->FirstChildElement();
      for (; opt; opt = opt->NextSiblingElement())
        if (!strcasecmp(opt->Value(),"planestrain"))
          planeStrain = true;
        else if (!strcasecmp(opt->Value(),"axisymmetric"))
          axiSymmetry = true;

      mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
    }

    else if (!strcasecmp(child->Value(),"bodyforce")) {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,12);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (child->FirstChild() && code > 0) {
        utl::getAttribute(child,"type",type,true);
        std::cout <<"\tBodyforce code "<< code;
        if (!type.empty()) std::cout <<" ("<< type <<")";
        VecFunc* f = utl::parseVecFunc(child->FirstChild()->Value(),type);
        if (f) this->setVecProperty(code,Property::BODYLOAD,f);
        std::cout << std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "hole") {
        double a = 0.0, F0 = 0.0, nu = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"F0",F0);
        utl::getAttribute(child,"nu",nu);
        std::cout <<"\tAnalytical solution: Hole a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
        if (!mySol)
          mySol = new AnaSol(new Hole(a,F0,nu));
      }
      else if (type == "lshape") {
        double a = 0.0, F0 = 0.0, nu = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"F0",F0);
        utl::getAttribute(child,"nu",nu);
        std::cout <<"\tAnalytical solution: Lshape a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
        if (!mySol)
          mySol = new AnaSol(new Lshape(a,F0,nu));
      }
      else if (type == "cants") {
        double L = 0.0, H = 0.0, F0 = 0.0;
        utl::getAttribute(child,"L",L);
        utl::getAttribute(child,"H",H);
        utl::getAttribute(child,"F0",F0);
        std::cout <<"\tAnalytical solution: CanTS L="<< L <<" H="<< H
                  <<" F0="<< F0 << std::endl;
        if (!mySol)
          mySol = new AnaSol(new CanTS(L,H,F0));
      }
      else if (type == "cantm") {
        double L = 0.0, H = 0.0, M0 = 0.0;
        utl::getAttribute(child,"L",L);
        utl::getAttribute(child,"H",H);
        utl::getAttribute(child,"M0",M0);
        std::cout <<"\tAnalytical solution: CanTM H="<< H
                  <<" M0="<< M0 << std::endl;
        if (!mySol)
          mySol = new AnaSol(new CanTM(H,M0));
      }
      else if (type == "curved") {
        double a = 0.0, b = 0.0, u0 = 0.0, E = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"b",b);
        utl::getAttribute(child,"u0",u0);
        utl::getAttribute(child,"E",E);
        std::cout <<"\tAnalytical solution: Curved Beam a="<< a <<" b="<< b
                  <<" u0="<< u0 <<" E="<< E << std::endl;
        if (!mySol)
          mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
      }
      else if (type == "expression") {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(child,false);
      }
      else
        std::cerr <<"  ** SIMLinEl2D::parse: Invalid analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      int code = 0;
      utl::getAttribute(child,"code",code);
      if (code > 0 && mySol && mySol->getStressSol()) {
        std::cout <<"\tNeumann code "<< code
                  <<": Analytical traction"<< std::endl;
        this->setPropertyType(code,Property::NEUMANN);
        myTracs[code] = new TractionField(*mySol->getStressSol());
      }
    }

  if (gx != 0.0 || gy != 0.0)
    this->getIntegrand()->setGravity(gx,gy);
  if (!mVec.empty())
    this->getIntegrand()->setMaterial(mVec.front());

  return true;
}


Elasticity* SIMLinEl2D::getIntegrand ()
{
  if (myProblem)
    return dynamic_cast<Elasticity*>(myProblem);

  Elasticity* elp = new LinearElasticity(2,axiSymmetry);
  myProblem = elp;

  return elp;
}


bool SIMLinEl2D::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!myProblem)
  {
    this->getIntegrand();
    if (myPid == 0) this->printProblem(std::cout);
  }

  if (mySol) // Define analytical boundary condition fields
    for (PropertyVec::iterator p = myProps.begin(); p != myProps.end(); p++)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        if (!mySol->getVectorSol())
          p->pcode = Property::UNDEFINED;
        else if (aCode == abs(p->pindx))
          p->pcode = Property::DIRICHLET_INHOM;
        else if (aCode == 0)
        {
          aCode = abs(p->pindx);
          myVectors[aCode] = mySol->getVectorSol();
          p->pcode = Property::DIRICHLET_INHOM;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        if (mySol->getStressSol())
        {
          p->pcode = Property::NEUMANN;
          myTracs[p->pindx] = new TractionField(*mySol->getStressSol());
        }
        else
          p->pcode = Property::UNDEFINED;
      }

  return this->SIM2D::preprocess(ignored,fixDup);
}


bool SIMLinEl2D::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd]);
  return true;
}


bool SIMLinEl2D::initBodyLoad (size_t patchInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  VecFuncMap::const_iterator it = myVectors.end();
  for (size_t i = 0; i < myProps.size(); i++)
    if (myProps[i].pcode == Property::BODYLOAD && myProps[i].patch == patchInd)
      if ((it = myVectors.find(myProps[i].pindx)) != myVectors.end()) break;

  elp->setBodyForce(it == myVectors.end() ? NULL : it->second);
  return true;
}


bool SIMLinEl2D::initNeumann (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  VecFuncMap::const_iterator  vit = myVectors.find(propInd);
  TracFuncMap::const_iterator tit = myTracs.find(propInd);

  if (vit != myVectors.end())
    elp->setTraction(vit->second);
  else if (tit != myTracs.end())
    elp->setTraction(tit->second);
  else
    return false;

  return true;
}


std::ostream& SIMLinEl2D::printNorms (const Vectors& norms, std::ostream& os)
{
  if (norms.empty()) return os;

  NormBase* norm = this->getNormIntegrand();
  const Vector& gnorm = norms.front();

  os <<"Energy norm "<< norm->getName(1,1) <<": "<< gnorm(1)
     <<"\nExternal energy "<< norm->getName(1,2) <<": "<< gnorm(2);

  if (mySol)
    os <<"\nExact norm "<< norm->getName(1,3) <<": "<< gnorm(3)
       <<"\nExact error "<< norm->getName(1,4) <<": "<< gnorm(4)
       <<"\nExact relative error (%) : "<< 100.0*gnorm(4)/gnorm(3);

  delete norm;

  return os << std::endl;
}
