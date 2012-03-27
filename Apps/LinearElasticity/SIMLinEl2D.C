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


SIMLinEl2D::SIMLinEl2D (int form) : SIM2D(2)
{
  if (form == SIM::NONE)
    return;
  else
    myProblem = new LinearElasticity(2,axiSymmetry);
}


SIMLinEl2D::~SIMLinEl2D ()
{
  for (size_t i = 0; i < mVec.size(); i++)
    delete mVec[i];
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
    elp->setGravity(gx,gy);
    if (myPid == 0)
      std::cout <<"\nGravitation vector: "
		<< gx <<" "<< gy << std::endl;
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
      double rho = atof(strtok(NULL," "));
      mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho << std::endl;
    }
    if (!mVec.empty())
      elp->setMaterial(mVec.front());
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
      mySol = new AnaSol(new Hole(a,F0,nu));
      std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      mySol = new AnaSol(new Lshape(a,F0,nu));
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      mySol = new AnaSol(new CanTS(L,H,F0));
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		<<" F0="<< F0 << std::endl;
    }
    else if (!strncasecmp(cline,"CANTM",5))
    {
      double H  = atof(strtok(NULL," "));
      double M0 = atof(strtok(NULL," "));
      mySol = new AnaSol(new CanTM(H,M0));
      std::cout <<"\nAnalytical solution: CanTM H="<< H
		<<" M0="<< M0 << std::endl;
    }
    else if (!strncasecmp(cline,"CURVED",6))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double u0 = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
      std::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
                <<" u0="<< u0 <<" E="<< E << std::endl;
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinEl2D::parse: Unknown analytical solution "
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

    if (!mVec.empty())
      elp->setMaterial(mVec.front());
  }

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMLinEl2D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIM2D::parse(elem);

  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity")) {
      double gx = 0.0, gy = 0.0;
      utl::getAttribute(child,"x",gx);
      utl::getAttribute(child,"y",gy);
      elp->setGravity(gx,gy);
      if (myPid == 0)
        std::cout <<"\nGravitation vector: " << gx <<" "<< gy << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      int code = 0;
      if (utl::getAttribute(child,"code",code) && code > 0)
        setPropertyType(code,Property::MATERIAL,mVec.size());
      double E = 1000.0, nu = 0.3, rho = 1.0;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
    }

    else if (!strcasecmp(child->Value(),"constantpressure") ||
             !strcasecmp(child->Value(),"linearpressure")) {
      if (child->FirstChild() && child->FirstChild()->Value()) {
        double p = atof(child->FirstChild()->Value());
        int code = 0, pdir = 1;
        utl::getAttribute(child,"code",code);
        utl::getAttribute(child,"dir",pdir);
        setPropertyType(code,Property::NEUMANN);
        if (!strcasecmp(child->Value(),"linearpressure")) {
          RealFunc* pfl = new ConstTimeFunc(new LinearFunc(p));
          myTracs[code] = new PressureField(pfl,pdir);
        }
        else
          myTracs[code] = new PressureField(p,pdir);
        if (myPid == 0)
          std::cout <<"\tPressure code "<< code <<" direction "<< pdir
                    <<": "<< p << std::endl;
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
        mySol = new AnaSol(new Hole(a,F0,nu));
        std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
      }
      else if (type == "lshape") {
        double a = 0.0, F0 = 0.0, nu = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"F0",F0);
        utl::getAttribute(child,"nu",nu);
        mySol = new AnaSol(new Lshape(a,F0,nu));
        std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
      }
      else if (type == "cants") {
        double L = 0.0, H = 0.0, F0 = 0.0;
        utl::getAttribute(child,"L",L);
        utl::getAttribute(child,"H",H);
        utl::getAttribute(child,"F0",F0);
        mySol = new AnaSol(new CanTS(L,H,F0));
        std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
                  <<" F0="<< F0 << std::endl;
      }
      else if (type == "cantm") {
        double L = 0.0, H = 0.0, M0 = 0.0;
        utl::getAttribute(child,"L",L);
        utl::getAttribute(child,"H",H);
        utl::getAttribute(child,"M0",M0);
        mySol = new AnaSol(new CanTM(H,M0));
        std::cout <<"\nAnalytical solution: CanTM H="<< H
                  <<" M0="<< M0 << std::endl;
      }
      else if (type == "curved") {
        double a = 0.0, b = 0.0, u0 = 0.0, E = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"b",b);
        utl::getAttribute(child,"u0",u0);
        utl::getAttribute(child,"E",E);
        mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
        std::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
                  <<" u0="<< u0 <<" E="<< E << std::endl;
      }
      else if (type == "expression") {
        std::cout <<"\nAnalytical solution: Expression"<< std::endl;
        mySol = new AnaSol(child,false);
      }
      else
        std::cerr <<"  ** SIMLinEl2D::parse: Unknown analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      int code = 0;
      utl::getAttribute(child,"code",code);
      if (code > 0 && mySol && mySol->getStressSol()) {
        std::cout <<"Pressure code "<< code
                  <<": Analytical traction"<< std::endl;
        setPropertyType(code,Property::NEUMANN);
        myTracs[code] = new TractionField(*mySol->getStressSol());
      }
    }

  if (!mVec.empty())
    elp->setMaterial(mVec.front());

  return true;
}


bool SIMLinEl2D::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd]);
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
