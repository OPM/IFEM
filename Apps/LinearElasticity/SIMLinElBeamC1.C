// $Id$
//==============================================================================
//!
//! \file SIMLinElBeamC1.C
//!
//! \date Sep 16 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of C1-continous beams.
//!
//==============================================================================

#include "SIMLinElBeamC1.h"
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


SIMLinElBeamC1::SIMLinElBeamC1 () : SIM1D(SIM::NONE)
{
  nf = 1;
  myProblem = new KirchhoffLovePlate(1);
}


bool SIMLinElBeamC1::parse (char* keyWord, std::istream& is)
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
      double rho = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      double thk = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      mVec.push_back(new LinIsotropic(E,0.0,rho,true));
      tVec.push_back(thk);
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< rho <<" "<< thk << std::endl;
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
      myLoads[i].xi    = atof(strtok(NULL," "));
      myLoads[i].pload = atof(strtok(NULL," "));
      if (myPid == 0)
	std::cout <<"\n\tPoint "<< i+1 <<": P"<< myLoads[i].patch
		  <<" xi = "<< myLoads[i].xi <<" load = "<< myLoads[i].pload;
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
    if (mySol) return true;

    cline = strtok(keyWord+6," ");
    /*
    if (!strncasecmp(cline,"NAVIERBEAM",7))
    {
      double a  = atof(strtok(NULL," "));
      double t  = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      double pz = atof(strtok(NULL," "));
      std::cout <<"\nAnalytic solution: NavierPlate a="<< a
		<<" t="<< t <<" E="<< E <<" pz="<< pz;
      if ((cline = strtok(NULL," ")))
      {
	double xi = atof(cline);
	std::cout <<" xi="<< xi
	if ((cline = strtok(NULL," ")))
	{
	  double c = atof(cline);
	  std::cout <<" c="<< c <<" d="<< d;
	  mySol = new AnaSol(new NavierBeam(a,t,E,nu,pz,xi,c));
	}
	else
	  mySol = new AnaSol(new NavierBeam(a,t,E,nu,pz,xi));
      }
      else
	mySol = new AnaSol(new NavierBeam(a,t,E,pz));
    }
    */
    if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinElBeamC1::parse: Unknown analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }
  }

  else
    return this->SIM1D::parse(keyWord,is);

  return true;
}


bool SIMLinElBeamC1::parse(const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"eulerbernoulli"))
    return this->SIM1D::parse(elem);

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  std::vector<const TiXmlElement*> parsed = handlePriorityTags(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if (find(parsed.begin(),parsed.end(),child) != parsed.end())
      ;

    else if (!strcasecmp(child->Value(),"gravity")) {
      double g = 0.0;
      utl::getAttribute(child,"g",g);
      if (myPid == 0)
        std::cout <<"\nGravitation constant: "<< g << std::endl;
      klp->setGravity(g);
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      int code = 0;
      if (utl::getAttribute(child,"code",code) && code > 0)
        setPropertyType(code,Property::MATERIAL,mVec.size());
      double E = 1000.0, rho = 1.0, thk = 0.1;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"rho",rho);
      utl::getAttribute(child,"thickness",thk);
      mVec.push_back(new LinIsotropic(E,0.0,rho,true));
      tVec.push_back(thk);
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< rho <<" " << thk << std::endl;
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
        utl::getAttribute(child,"xi",load.xi);
        if (myPid == 0)
          std::cout <<"\n\tPoint: P"<< load.patch
                    <<" xi = "<< load.xi
                    <<" load = "<< load.pload;
        myLoads.push_back(load);
      }
    }

    else if (!strcasecmp(child->Value(),"pressure")) {
      int code = 0;
      double p = 0.0;
      char* fn = NULL;
      utl::getAttribute(child,"code",code);
      if (child->FirstChild())
        p = atof(child->FirstChild()->Value());
      if (myPid == 0)
        std::cout <<"\tPressure code "<< code <<": ";
      std::string function;
      if (utl::getAttribute(child,"function",function))
        fn = strtok(const_cast<char*>(function.c_str())," ");
      myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(fn,p));
      std::cout << std::endl;
      if (code > 0)
        setPropertyType(code,Property::BODYLOAD);
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "expression") {
        std::cout <<"\nAnalytical solution: Expression"<< std::endl;
        mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
                  << type <<" (ignored)"<< std::endl;
    }

    child = child->NextSiblingElement();
  }

  return true;
}


bool SIMLinElBeamC1::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  klp->setMaterial(mVec[propInd]);
  if (tVec[propInd] != 0.0)
    klp->setThickness(tVec[propInd]);

  return true;
}


bool SIMLinElBeamC1::initBodyLoad (size_t patchInd)
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


bool SIMLinElBeamC1::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!this->SIMbase::preprocess(ignored,fixDup))
    return false;

  // Preprocess the nodal point loads
  for (PloadVec::iterator p = myLoads.begin(); p != myLoads.end();)
  {
    int pid = this->getLocalPatchIndex(p->patch);
    if (pid < 1 || myModel[pid-1]->empty())
      p = myLoads.erase(p);
    else if ((p->inod = myModel[pid-1]->evalPoint(&p->xi,&p->xi,p->X)) < 1)
    {
      p = myLoads.erase(p);
      std::cerr <<"  ** SIMLinElBeamC1::preprocess: Load point ("<< p->xi
		<<") on patch #"<< p->patch <<" is not a nodal point"
		<<" (ignored)." << std::endl;
    }
    else
    {
      int ipt = 1 + (int)(p-myLoads.begin());
      if (ipt == 1) std::cout <<'\n';
      std::cout <<"Load point #"<< ipt <<": patch #"<< p->patch
		<<" (u,v)=("<< p->xi <<"), node #"<< p->inod
		<<", X = "<< p->X << std::endl;
      p++;
    }
  }

  return true;
}


bool SIMLinElBeamC1::finalizeAssembly (bool newLHSmatrix)
{
  SystemVector* b = myEqSys->getVector();
  for (size_t i = 0; i < myLoads.size() && b; i++)
    if (!mySam->assembleSystem(*b,&myLoads[i].pload,myLoads[i].inod))
      return false;

  return this->SIMbase::finalizeAssembly(newLHSmatrix);
}


double SIMLinElBeamC1::externalEnergy (const Vectors& psol) const
{
  double energy = this->SIMbase::externalEnergy(psol);

  // External energy from the nodal point loads
  for (size_t i = 0; i < myLoads.size(); i++)
    energy += myLoads[i].pload * psol.front()(myLoads[i].inod);

  return energy;
}
