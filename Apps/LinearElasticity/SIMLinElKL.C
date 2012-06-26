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
#include "AlgEqSystem.h"
#include "ASMbase.h"
#include "SAMpatch.h"
#include "Functions.h"
#include "Utilities.h"
#include "Property.h"
#include "AnaSol.h"
#include "Vec3Oper.h"
#include "tinyxml.h"


SIMLinElKL::SIMLinElKL ()
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

  else if (!strncasecmp(keyWord,"POINTLOAD",9))
  {
    int nload = atoi(keyWord+9);
    std::cout <<"\nNumber of point loads: "<< nload;

    myLoads.resize(nload);
    for (int i = 0; i < nload && (cline = utl::readLine(is)); i++)
    {
      myLoads[i].patch = atoi(strtok(cline," "));
      myLoads[i].xi[0] = atof(strtok(NULL," "));
      myLoads[i].xi[1] = atof(strtok(NULL," "));
      myLoads[i].pload = atof(strtok(NULL," "));
      if (myPid == 0)
	std::cout <<"\n\tPoint "<< i+1 <<": P"<< myLoads[i].patch
		  <<" xi = "<< myLoads[i].xi[0] <<" "<< myLoads[i].xi[1]
		  <<" load = "<< myLoads[i].pload;
    }
  }

  else if (!strncasecmp(keyWord,"PRESSURE",8))
  {
    int npres = atoi(keyWord+8);
    std::cout <<"\nNumber of pressures: "<< npres << std::endl;

    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      double p = atof(strtok(NULL," "));
      std::cout <<"\tPressure code "<< code <<": ";
      cline = strtok(NULL," ");
      myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,p));
      std::cout << std::endl;
      if (code > 0)
	this->setPropertyType(code,Property::BODYLOAD);
    }
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"NAVIERPLATE",7))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double t  = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      double pz = atof(strtok(NULL," "));
      std::cout <<"\nAnalytic solution: NavierPlate a="<< a <<" b="<< b
		<<" t="<< t <<" E="<< E <<" nu="<< nu <<" pz="<< pz;
      if ((cline = strtok(NULL," ")))
      {
	double xi  = atof(cline);
	double eta = atof(strtok(NULL," "));
	std::cout <<" xi="<< xi <<" eta="<< eta;
	if ((cline = strtok(NULL," ")))
	{
	  double c = atof(cline);
	  double d = atof(strtok(NULL," "));
	  std::cout <<" c="<< c <<" d="<< d;
	  if (!mySol)
	    mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d));
	}
	else if (!mySol)
	  mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta));
      }
      else if (!mySol)
	mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (!mySol)
	mySol = new AnaSol(is,lines,false);
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


bool SIMLinElKL::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"kirchhofflove"))
    return this->SIM2D::parse(elem);

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity")) {
      double g = 0.0;
      utl::getAttribute(child,"g",g);
      klp->setGravity(g);
      if (myPid == 0)
        std::cout <<"\nGravitation constant: "<< g << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      std::string set;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,0);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E = 1000.0, nu = 0.3, rho = 1.0, thk = 0.1;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      utl::getAttribute(child,"thickness",thk);

      mVec.push_back(new LinIsotropic(E,nu,rho,true));
      tVec.push_back(thk);
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" " << thk << std::endl;
      klp->setMaterial(mVec.front());
      if (tVec.front() != 0.0)
        klp->setThickness(tVec.front());
    }

    else if (!strcasecmp(child->Value(),"pointload")) {
      PointLoad load; int patch;
      utl::getAttribute(child,"patch",patch);
      if (child->FirstChild()) {
        load.patch = patch;
        load.pload = atof(child->FirstChild()->Value());
        utl::getAttribute(child,"xi",load.xi[0]);
        utl::getAttribute(child,"eta",load.xi[1]);
        if (myPid == 0)
          std::cout <<"\n\tPoint: P"<< load.patch
                    <<" xi = "<< load.xi[0] <<" "<< load.xi[1]
                    <<" load = "<< load.pload;
        myLoads.push_back(load);
      }
    }

    else if (!strcasecmp(child->Value(),"pressure")) {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);

      if (child->FirstChild() && code > 0) {
        utl::getAttribute(child,"type",type,true);
        std::cout <<"\tPressure code "<< code;
        if (!type.empty()) std::cout <<" ("<< type <<")";
        myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),type);
        this->setPropertyType(code,Property::BODYLOAD);
      }
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "navierplate") {
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0, t = 0.0;
        double E = 10000.0, nu = 0.3, pz = 1.0, xi = 0.0, eta = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"b",b);
        utl::getAttribute(child,"c",c);
        utl::getAttribute(child,"d",d);
        utl::getAttribute(child,"t",t);
        utl::getAttribute(child,"E",E);
        utl::getAttribute(child,"nu",nu);
        utl::getAttribute(child,"pz",pz);
        utl::getAttribute(child,"xi",xi);
        utl::getAttribute(child,"eta",eta);
        std::cout <<"\nAnalytic solution: NavierPlate a="<< a <<" b="<< b
                  <<" t="<< t <<" E="<< E <<" nu="<< nu <<" pz="<< pz;
        if (xi != 0.0 && eta != 0.0) {
          std::cout <<" xi="<< xi <<" eta="<< eta;
          if (c != 0.0 && d != 0.0) {
            std::cout <<" c="<< c <<" d="<< d;
            if (!mySol)
              mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d));
          }
          else if (!mySol)
            mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta));
        }
        else if (!mySol)
          mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz));
      }
      else if (type == "expression") {
        std::cout <<"\nAnalytical solution: Expression"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
                  << type <<" (ignored)"<< std::endl;
    }

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


bool SIMLinElKL::initBodyLoad (size_t patchInd)
{
  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  SclFuncMap::const_iterator it = myScalars.find(0);
  for (size_t i = 0; i < myProps.size(); i++)
    if (myProps[i].pcode == Property::BODYLOAD && myProps[i].patch == patchInd)
      if ((it = myScalars.find(myProps[i].pindx)) != myScalars.end()) break;

  klp->setPressure(it == myScalars.end() ? 0 : it->second);
  return true;
}


bool SIMLinElKL::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!this->SIMLinEl2D::preprocess(ignored,fixDup))
    return false;

  // Preprocess the nodal point loads
  for (PloadVec::iterator p = myLoads.begin(); p != myLoads.end();)
  {
    int pid = this->getLocalPatchIndex(p->patch);
    if (pid < 1 || myModel[pid-1]->empty())
      p = myLoads.erase(p);
    else if ((p->inod = myModel[pid-1]->evalPoint(p->xi,p->xi,p->X)) < 1)
    {
      p = myLoads.erase(p);
      std::cerr <<"  ** SIMLinElKL::preprocess: Load point ("
		<< p->xi[0] <<','<< p->xi[1]
		<<") on patch #"<< p->patch <<" is not a nodal point"
		<<" (ignored)." << std::endl;
    }
    else
    {
      int ipt = 1 + (int)(p-myLoads.begin());
      if (ipt == 1) std::cout <<'\n';
      std::cout <<"Load point #"<< ipt <<": patch #"<< p->patch
		<<" (u,v)=("<< p->xi[0] <<','<< p->xi[1]
		<<"), node #"<< p->inod <<", X = "<< p->X << std::endl;
      p++;
    }
  }

  return true;
}


bool SIMLinElKL::assembleDiscreteTerms (const IntegrandBase* problem)
{
  if (problem != myProblem)
    return true; // Do this only for the main integrand

  SystemVector* b = myEqSys->getVector();
  for (size_t i = 0; i < myLoads.size() && b; i++)
    if (!mySam->assembleSystem(*b,&myLoads[i].pload,myLoads[i].inod))
      return false;

  return true;
}


double SIMLinElKL::externalEnergy (const Vectors& psol) const
{
  double energy = this->SIMbase::externalEnergy(psol);

  // External energy from the nodal point loads
  for (size_t i = 0; i < myLoads.size(); i++)
    energy += myLoads[i].pload * psol.front()(myLoads[i].inod);

  return energy;
}
