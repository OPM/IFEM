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
#include <string.h>

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
    cline = strtok(keyWord+6," ");
    int code = -1;
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
      std::string variables, stress;
      int lines = atoi(strtok(NULL, " "));
      char* c = strtok(NULL, " ");
      if (c)
        code = atoi(c);
      else
        code = 0;
      for (int i = 0; i < lines; i++) {
        std::string function = utl::readLine(is);
        size_t pos;
        if ((pos = function.find("Variables=")) != std::string::npos) {
          variables += function.substr(pos+10);
          if (variables[variables.size()-1] != ';')
            variables += ";";
        }
        if ((pos = function.find("Stress=")) != std::string::npos) {
          stress = function.substr(pos+7);
	  mySol = new AnaSol(new STensorFuncExpr(stress,variables));
	  std::cout <<"\nAnalytical solution:";
	  if (!variables.empty())
	    std::cout <<"\n\tVariables = "<< variables;
	  std::cout <<"\n\tStress = "<< stress << std::endl;
	  break;
        }
      }
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


bool SIMLinEl2D::parse(const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"linearelasticity"))
    return SIM2D::parse(elem);

  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp)
    return false;

  std::vector<const TiXmlElement*> parsed = handlePriorityTags(elem);
  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if (find(parsed.begin(),parsed.end(),child) != parsed.end()) {
      child = child->NextSiblingElement();
      continue;
    }
    if (!strcasecmp(child->Value(),"gravity")) {
      double gx=0, gy=0;
      if (child->Attribute("x"))
        gx = atof(child->Attribute("x"));
      if (child->Attribute("y"))
        gy = atof(child->Attribute("y"));
      elp->setGravity(gx,gy);
      if (myPid == 0)
        std::cout <<"\nGravitation vector: " << gx <<" "<< gy << std::endl;
    } else if (!strcasecmp(child->Value(),"isotropic")) {
      int code=0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      if (code > 0)
        setPropertyType(code,Property::MATERIAL,mVec.size());
      double E=1000.f, nu=0.3f, rho = 1.f;
      if (child->Attribute("E"))
        E = atof(child->Attribute("E"));
      if (child->Attribute("rho"))
        rho = atof(child->Attribute("rho"));
      if (child->Attribute("nu"))
        nu = atof(child->Attribute("nu"));
      mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
      if (!mVec.empty())
        elp->setMaterial(mVec.front());
    } else if (!strcasecmp(child->Value(),"constantpressure")) {
      int code=0, pdir=1;
      double p=0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      if (child->Attribute("dir"))
        pdir = atoi(child->Attribute("dir"));
      if (child->FirstChild() && child->FirstChild()->Value())
        p = atof(child->FirstChild()->Value());
      if (myPid == 0)
	std::cout <<"\tPressure code "<< code <<" direction "<< pdir
		  <<": "<< p << std::endl;
      setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new PressureField(p,pdir);
    } else if (!strcasecmp(child->Value(),"anasol")) {
      if (child->Attribute("type") &&
          !strcasecmp(child->Attribute("type"),"hole")) {
        double a=0, F0 = 0, nu = 0;
        if (child->Attribute("a"))
          a = atof(child->Attribute("a"));
        if (child->Attribute("F0"))
          F0 = atof(child->Attribute("F0"));
        if (child->Attribute("nu"))
          nu = atof(child->Attribute("nu"));
        mySol = new AnaSol(new Hole(a,F0,nu));
        std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"Lshape")) {
        double a=0, F0 = 0, nu = 0;
        if (child->Attribute("a"))
          a = atof(child->Attribute("a"));
        if (child->Attribute("F0"))
          F0 = atof(child->Attribute("F0"));
        if (child->Attribute("nu"))
          nu = atof(child->Attribute("nu"));
        mySol = new AnaSol(new Lshape(a,F0,nu));
        std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"cants")) {
        double L=0, F0 = 0, H = 0;
        if (child->Attribute("L"))
          L = atof(child->Attribute("L"));
        if (child->Attribute("F0"))
          F0 = atof(child->Attribute("F0"));
        if (child->Attribute("H"))
          H = atof(child->Attribute("H"));
        mySol = new AnaSol(new CanTS(L,H,F0));
        std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
                  <<" F0="<< F0 << std::endl;
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"cantm")) {
        double M0 = 0, H = 0;
        if (child->Attribute("M0"))
          M0 = atof(child->Attribute("M0"));
        if (child->Attribute("H"))
          H = atof(child->Attribute("H"));
        mySol = new AnaSol(new CanTM(H,M0));
        std::cout <<"\nAnalytical solution: CanTM H="<< H
                  <<" M0="<< M0 << std::endl;
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"curved")) {
        double a=0, b=0, u0 = 0, E = 0;
        if (child->Attribute("a"))
          a = atof(child->Attribute("a"));
        if (child->Attribute("b"))
          b = atof(child->Attribute("b"));
        if (child->Attribute("u0"))
          u0 = atof(child->Attribute("u0"));
        if (child->Attribute("E"))
          E = atof(child->Attribute("E"));
        mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
        std::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
                  <<" u0="<< u0 <<" E="<< E << std::endl;
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"expression")) {
        std::string variables, stress;
        const TiXmlElement* var = child->FirstChildElement("variables");
        if (var && var->FirstChild() && var->FirstChild()->Value()) {
          variables = var->FirstChild()->Value();
          if (variables[variables.size()-1] != ';')
            variables += ";";
        }
        const TiXmlElement* str = child->FirstChildElement("stress");
        if (str && str->FirstChild() && str->FirstChild()->Value())
          stress = str->FirstChild()->Value();
        mySol = new AnaSol(new STensorFuncExpr(stress,variables));
	std::cout <<"\nAnalytical solution:";
	if (!variables.empty())
	  std::cout <<"\n\tVariables = "<< variables;
	std::cout <<"\n\tStress = "<< stress << std::endl;
      } else {
        std::cerr <<"  ** SIMLinEl2D::parse: Unknown analytical solution "
          << (child->Attribute("type")?child->Attribute("type"):"") << std::endl;
      }

      // Define the analytical boundary traction field
      int code=0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      if (code > 0 && mySol && mySol->getStressSol()) {
        std::cout <<"Pressure code "<< code <<": Analytical traction"<< std::endl;
        setPropertyType(code,Property::NEUMANN);
        myTracs[code] = new TractionField(*mySol->getStressSol());
      }
    }
    else
      SIM2D::parse(child);

    child = child->NextSiblingElement();
  }

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
