// $Id$
//==============================================================================
//!
//! \file SIMLinEl3D.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMLinEl3D.h"
#include "LinIsotropic.h"
#include "LinearElasticity.h"
#include "AnalyticSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Property.h"
#include "Tensor.h"
#include "AnaSol.h"
#include "tinyxml.h"


/*!
  \brief Local coordinate system for a cylinder along global z-axis.
*/

class CylinderCS : public LocalSystem
{
public:
  //! \brief Constructor printing a message making user aware of its presense.
  CylinderCS()
  {
    std::cout <<"\nLocal coordinate system: Cylindric"<< std::endl;
  }

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
    double r = hypot(X.x,X.y);
    T(1,1) = X.x/r;
    T(1,2) = X.y/r;
    T(2,1) = -T(1,2);
    T(2,2) = T(1,1);
    T(3,3) = 1.0;
    return T;
  }
};


#ifdef PRINT_CS
#include <fstream>
#endif

/*!
  \brief Local coordinate system for a cylinder along global z-axis,
  closed by a spherical cap.
*/

class CylinderSphereCS : public LocalSystem
{
public:
  //! \brief Constructor printing a message making user aware of its presense.
  CylinderSphereCS(double H = 0.0) : h(H)
  {
    std::cout <<"\nLocal coordinate system: Cylindric with Spherical cap, h="
	      << h << std::endl;
#ifdef PRINT_CS
    sn.open("nodes.dat");
    se.open("elements.dat");
    s1.open("v1.dat");
    s2.open("v2.dat");
    s3.open("v3.dat");
    sn <<"\n*NODES 4\n";
    se <<"\n*ELEMENTS 4\n%NODES #4\n"
       <<"%NO_ID\n%MAP_NODE_INDICES\n%PART_ID 4\n%POINTS\n";
    s1 <<"\n*RESULTS 31\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s2 <<"\n*RESULTS 32\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s3 <<"\n*RESULTS 33\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
  }

  virtual ~CylinderSphereCS()
  {
    std::cout <<"Finalizing VTF-output of local coordinate systems"<< std::endl;
    s1 <<"\n*GLVIEWVECTOR 2\n%NAME \"v1\"\n%STEP 1\n31\n";
    s2 <<"\n*GLVIEWVECTOR 3\n%NAME \"v2\"\n%STEP 1\n32\n";
    s3 <<"\n*GLVIEWVECTOR 4\n%NAME \"v3\"\n%STEP 1\n33\n";
#endif
  }

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
#ifdef PRINT_CS
    sn << X <<'\n';
    static int iel = 0;
    se << ++iel <<'\n';
#endif
    if (patch == 1) // Cylindric system {-z,theta,r}
    {
      T.zero();
      double r = hypot(X.x,X.y);
      T(1,3) = -1.0;
      T(2,1) = -X.y/r;
      T(2,2) =  X.x/r;
      T(3,1) =  T(2,2);
      T(3,2) = -T(2,1);
#ifdef PRINT_CS
      s1 <<"0 0 -1\n";
      s2 << T(2,1) <<" "<< T(2,2) <<" 0\n";
      s3 << T(3,1) <<" "<< T(3,2) <<" 0\n";
#endif
    }
    else // Spherical system {phi,theta,r}
    {
      Vec3 v3(X.x,X.y,X.z-h);
      v3 /= v3.length();
      double theta = atan2(X.y,X.x);
      double phi = acos(v3.z);
      Vec3 v1(cos(theta)*cos(phi),sin(theta)*cos(phi),-sin(phi));
      Vec3 v2(v3,v1);
      for (int i = 1; i <= 3; i++)
      {
	T(1,i) = v1[i-1];
	T(2,i) = v2[i-1];
	T(3,i) = v3[i-1];
      }
#ifdef PRINT_CS
      s1 << v1 <<'\n';
      s2 << v2 <<'\n';
      s3 << v3 <<'\n';
#endif
    }
    return T;
  }

private:
  double h; //!< Height above global origin of the centre of the sphere
#ifdef PRINT_CS
  mutable std::ofstream sn; //!< VTF output stream for CS nodes
  mutable std::ofstream se; //!< VTF output stream for CS point elements
  mutable std::ofstream s1; //!< VTF output stream for vector v1 of local CS
  mutable std::ofstream s2; //!< VTF output stream for vector v2 of local CS
  mutable std::ofstream s3; //!< VTF output stream for vector v3 of local CS
#endif
};


SIMLinEl3D::~SIMLinEl3D ()
{
  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode > 0) myVectors.erase(aCode);

  for (size_t i = 0; i < mVec.size(); i++)
    delete mVec[i];
}


void SIMLinEl3D::clearProperties ()
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


bool SIMLinEl3D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  int nmat = 0;
  int nConstPress = 0;
  int nLinearPress = 0;
  double gx = 0.0, gy = 0.0, gz = 0;

  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    gx = atof(strtok(keyWord+7," "));
    gy = atof(strtok(NULL," "));
    gz = atof(strtok(NULL," "));
    if (myPid == 0)
      std::cout <<"\nGravitation vector: "
		<< gx <<" "<< gy <<" "<< gz << std::endl;
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
      mVec.push_back(new LinIsotropic(E,nu,rho));
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"CONSTANT_PRESSURE",17))
    nConstPress  = atoi(keyWord+17);
  else if (!strncasecmp(keyWord,"LINEAR_PRESSURE",15))
    nLinearPress = atoi(keyWord+15);

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
        mySol = new AnaSol(new Hole(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		<<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0,true));
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
      std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
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
    press.ldim = 2;

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
      if (press.lindx < 1 || press.lindx > 6)
      {
	std::cerr <<" *** SIMLinEl3D::parse: Invalid face index "
		  << (int)press.lindx << std::endl;
	return false;
      }

      if (mySol && mySol->getStressSol())
      {
	std::cout <<"\tTraction on P"<< press.patch
		  <<" F"<< (int)press.lindx << std::endl;
	myTracs[1+i] = new TractionField(*mySol->getStressSol());
      }
      else
      {
	int pdir = atoi(strtok(NULL," "));
	double p = atof(strtok(NULL," "));
	std::cout <<"\tPressure on P"<< press.patch
		  <<" F"<< (int)press.lindx <<" direction "<< pdir <<": ";
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
	  mVec.push_back(new LinIsotropic(E,nu,rho));
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
	  mVec.push_back(new LinIsotropic(E,nu,rho));
	}
    }
  }

  else if (!strncasecmp(keyWord,"LOCAL_SYSTEM",12))
  {
    size_t i = 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    if (!strncasecmp(keyWord+i,"CYLINDRICZ",10))
      this->getIntegrand()->setLocalSystem(new CylinderCS);
    else if (!strncasecmp(keyWord+i,"CYLINDER+SPHERE",15))
    {
      double H = atof(keyWord+i+15);
      this->getIntegrand()->setLocalSystem(new CylinderSphereCS(H));
    }
    else
      std::cerr <<"  ** SIMLinEl3D::parse: Unsupported coordinate system: "
		<< keyWord+i <<" (ignored)"<< std::endl;
  }

  else
    return this->SIM3D::parse(keyWord,is);

  int npres = nConstPress + nLinearPress;
  if (npres > 0)
  {
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

      if (nLinearPress)
      {
	RealFunc* pfl = new ConstTimeFunc(new LinearFunc(p));
	myTracs[code] = new PressureField(pfl,pdir);
      }
      else
	myTracs[code] = new PressureField(p,pdir);
    }
  }

  if (gx != 0.0 || gy != 0.0 || gz != 0.0)
    this->getIntegrand()->setGravity(gx,gy,gz);
  if (nmat > 0 && !mVec.empty())
    this->getIntegrand()->setMaterial(mVec.front());

  return true;
}


bool SIMLinEl3D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIM3D::parse(elem);

  Vec3 g;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity")) {
      utl::getAttribute(child,"x",g.x);
      utl::getAttribute(child,"y",g.y);
      utl::getAttribute(child,"z",g.z);
      if (myPid == 0)
        std::cout <<"\tGravitation vector: "<< g << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      int code = this->parseMaterialSet(child,mVec.size());

      double E = 1000.0, nu = 0.3, rho = 1.0;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);

      mVec.push_back(new LinIsotropic(E,nu,rho));
      if (myPid == 0)
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
    }

    else if (!strcasecmp(child->Value(),"bodyforce")) {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,123);
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
          mySol = new AnaSol(new Hole(a,F0,nu,true));
      }
      else if (type == "lshape") {
        double a = 0.0, F0 = 0.0, nu = 0.0;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"F0",F0);
        utl::getAttribute(child,"nu",nu);
        std::cout <<"\tAnalytical solution: Lshape a="<< a <<" F0="<< F0
                  <<" nu="<< nu << std::endl;
        if (!mySol)
          mySol = new AnaSol(new Lshape(a,F0,nu,true));
      }
      else if (type == "cants") {
        double L = 0.0, H = 0.0, F0 = 0.0;
        utl::getAttribute(child,"L",L);
        utl::getAttribute(child,"H",H);
        utl::getAttribute(child,"F0",F0);
        std::cout <<"\tAnalytical solution: CanTS L="<< L <<" H="<< H
                  <<" F0="<< F0 << std::endl;
        if (!mySol)
          mySol = new AnaSol(new CanTS(L,H,F0,true));
      }
      else if (type == "expression") {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(child,false);
      }
      else
        std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      int code = 0;
      utl::getAttribute(child,"code",code);
      if (code > 0 && mySol && mySol->getStressSol()) {
        std::cout <<"\tNeumann code "<< code
                  <<": Analytical traction"<< std::endl;
        setPropertyType(code,Property::NEUMANN);
        myTracs[code] = new TractionField(*mySol->getStressSol());
      }
    }

    else if (!strcasecmp(child->Value(),"localsystem") && child->FirstChild()) {
      // Caution: When running adaptively, the below will cause a small memory
      // leak because the coordinate system objects are only deleted by the
      // Elasticity destructor (and not in SIMbase::clearProperties).
      if (!strcasecmp(child->FirstChild()->Value(),"cylindricz"))
        this->getIntegrand()->setLocalSystem(new CylinderCS);
      else if (!strcasecmp(child->FirstChild()->Value(),"cylinder+sphere")) {
        double H = 0.0;
        utl::getAttribute(child,"H",H);
        this->getIntegrand()->setLocalSystem(new CylinderSphereCS(H));
      }
      else
        std::cerr <<"  ** SIMLinEl3D::parse: Unsupported coordinate system: "
                  << child->FirstChild()->Value() <<" (ignored)"<< std::endl;
    }

  if (!g.isZero(1.0e-16))
    this->getIntegrand()->setGravity(g.x,g.y,g.z);
  if (!mVec.empty())
    this->getIntegrand()->setMaterial(mVec.front());

  return true;
}


Elasticity* SIMLinEl3D::getIntegrand ()
{
  if (myProblem)
    return dynamic_cast<Elasticity*>(myProblem);

  Elasticity* elp = new LinearElasticity();
  myProblem = elp;

  return elp;
}


bool SIMLinEl3D::preprocess (const std::vector<int>& ignored, bool fixDup)
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

  return this->SIM3D::preprocess(ignored,fixDup);
}


bool SIMLinEl3D::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd]);
  return true;
}


bool SIMLinEl3D::initBodyLoad (size_t patchInd)
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


bool SIMLinEl3D::initNeumann (size_t propInd)
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


std::ostream& SIMLinEl3D::printNorms (const Vectors& norms, std::ostream& os)
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
