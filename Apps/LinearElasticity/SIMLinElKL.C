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
#include "Utilities.h"
#include "Functions.h"
#include "Property.h"
#include "AnaSol.h"
#include "Vec3Oper.h"
#include <string.h>

#include "tinyxml.h"


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
    if (mySol) return true;

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
	  mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d));
	}
	else
	  mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta));
      }
      else
	mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz));
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


bool SIMLinElKL::parse(const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"kirchhofflove"))
    return SIMLinEl2D::parse(elem);

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp)
    return false;

  std::vector<const TiXmlElement*> parsed = handlePriorityTags(elem);
  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if (find(parsed.begin(),parsed.end(),child) != parsed.end()) {
      child = child->NextSiblingElement();
      continue;
    }
    if (!strcasecmp(child->Value(),"gravity")) {
      double g=0;
      if (child->Attribute("g"))
        g = atof(child->Attribute("g"));
      klp->setGravity(g);
      if (myPid == 0)
        std::cout <<"\nGravitation constant: "<< g << std::endl;
    } else if (!strcasecmp(child->Value(),"isotropic")) {
      int code=0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      if (code > 0)
        setPropertyType(code,Property::MATERIAL,mVec.size());
      double E=1000.f, nu=0.3f, rho = 1.f, thk=0;
      if (child->Attribute("E"))
        E = atof(child->Attribute("E"));
      if (child->Attribute("rho"))
        rho = atof(child->Attribute("rho"));
      if (child->Attribute("nu"))
        nu = atof(child->Attribute("nu"));
      if (child->Attribute("thickness"))
        thk = atof(child->Attribute("thickness"));
      mVec.push_back(new LinIsotropic(E,nu,rho,true));
      tVec.push_back(thk);
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho <<" " << thk << std::endl;
      if (!mVec.empty())
        klp->setMaterial(mVec.front());
      if (!tVec.empty() && tVec.front() != 0.0)
        klp->setThickness(tVec.front());
    } else if (!strcasecmp(child->Value(),"pointload")) {
      PointLoad load;
      if (child->Attribute("patch"))
        load.patch = atoi(child->Attribute("patch"));
      if (child->Attribute("xi"))
        load.xi[0] = atof(child->Attribute("xi"));
      if (child->Attribute("eta"))
        load.xi[1] = atof(child->Attribute("eta"));
      if (child->FirstChild() && child->FirstChild()->Value())
        load.pload = atof(child->FirstChild()->Value());
        if (myPid == 0)
          std::cout <<"\n\tPoint: P"<< load.patch
                    <<" xi = "<< load.xi[0] <<" "<< load.xi[1]
                    <<" load = "<< load.pload;
        myLoads.push_back(load);
    } else if (!strcasecmp(child->Value(),"pressure")) {
      int code=0;
      double p=0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      if (child->FirstChild() && child->FirstChild()->Value())
        p = atof(child->FirstChild()->Value());
      if (myPid == 0)
        std::cout <<"\tPressure code "<< code <<": ";
      char* function;
      if (child->Attribute("function"))
        function = strdup(child->Attribute("function"));
      else
        function = strdup("");
      char* function2 = function;
      function = strtok(function," ");
      myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(const_cast<char*>(function),p));
      free(function2);
      std::cout << std::endl;
      if (code > 0)
	setPropertyType(code,Property::BODYLOAD);
    } else if (!strcasecmp(child->Value(),"anasol")) {
      if (child->Attribute("type") &&
          !strcasecmp(child->Attribute("type"),"navierplate")) {
        double a=0, b=0, c=0, d=0, t=0, E=0, nu=0, pz=0, xi=0, eta=0;
        if (child->Attribute("a"))
          a = atof(child->Attribute("a"));
        if (child->Attribute("b"))
          b = atof(child->Attribute("b"));
        if (child->Attribute("c"))
          c = atof(child->Attribute("c"));
        if (child->Attribute("d"))
          d = atof(child->Attribute("d"));
        if (child->Attribute("t"))
          t = atof(child->Attribute("t"));
        if (child->Attribute("E"))
          E = atof(child->Attribute("E"));
        if (child->Attribute("nu"))
          nu = atof(child->Attribute("nu"));
        if (child->Attribute("pz"))
          pz = atof(child->Attribute("pz"));
        if (child->Attribute("xi"))
          xi = atof(child->Attribute("xi"));
        if (child->Attribute("eta"))
          eta = atof(child->Attribute("eta"));
        std::cout <<"\nAnalytic solution: NavierPlate a="<< a <<" b="<< b
                  <<" t="<< t <<" E="<< E <<" nu="<< nu <<" pz="<< pz;

        if (xi != 0.f && eta != 0.f) {
          std::cout <<" xi="<< xi <<" eta="<< eta;
          if (c != 0.f && d != 0.f) {
            std::cout <<" c="<< c <<" d="<< d;
            mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d));
          } else
            mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz,xi,eta));
        } else
          mySol = new AnaSol(new NavierPlate(a,b,t,E,nu,pz));
      } else {
        std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
          << (child->Attribute("type")?child->Attribute("type"):"") << std::endl;
      }
    }
    else
      SIMLinEl2D::parse(child);

    child = child->NextSiblingElement();
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
  if (!this->SIMbase::preprocess(ignored,fixDup))
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


bool SIMLinElKL::finalizeAssembly (bool newLHSmatrix)
{
  SystemVector* b = myEqSys->getVector();
  for (size_t i = 0; i < myLoads.size() && b; i++)
    if (!mySam->assembleSystem(*b,&myLoads[i].pload,myLoads[i].inod))
      return false;

  return this->SIMbase::finalizeAssembly(newLHSmatrix);
}


double SIMLinElKL::externalEnergy (const Vectors& psol) const
{
  double energy = this->SIMbase::externalEnergy(psol);

  // External energy from the nodal point loads
  for (size_t i = 0; i < myLoads.size(); i++)
    energy += myLoads[i].pload * psol.front()(myLoads[i].inod);

  return energy;
}
