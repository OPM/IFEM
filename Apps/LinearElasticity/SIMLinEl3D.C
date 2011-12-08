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
#include <string.h>


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


SIMLinEl3D::SIMLinEl3D (bool checkRHS, int form) : SIM3D(checkRHS)
{
  if (form == SIM::NONE)
    return;
  else
    myProblem = new LinearElasticity();
}


SIMLinEl3D::~SIMLinEl3D ()
{
  for (size_t i = 0; i < mVec.size(); i++)
    delete mVec[i];
}


bool SIMLinEl3D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  int nConstPress = 0;
  int nLinearPress = 0;
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    double gx = atof(strtok(keyWord+7," "));
    double gy = atof(strtok(NULL," "));
    double gz = atof(strtok(NULL," "));
    elp->setGravity(gx,gy,gz);
    if (myPid == 0)
      std::cout <<"\nGravitation vector: "
		<< gx <<" "<< gy <<" "<< gz << std::endl;
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
      mVec.push_back(new LinIsotropic(E,nu,rho));
      if (myPid == 0)
	std::cout <<"\tMaterial code "<< code <<": "
		  << E <<" "<< nu <<" "<< rho << std::endl;
    }
    if (!mVec.empty())
      elp->setMaterial(mVec.front());
  }

  else if (!strncasecmp(keyWord,"CONSTANT_PRESSURE",17))
    nConstPress  = atoi(keyWord+17);
  else if (!strncasecmp(keyWord,"LINEAR_PRESSURE",15))
    nLinearPress = atoi(keyWord+15);

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    int code = -1;
    if (!strncasecmp(cline,"HOLE",4))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      mySol = new AnaSol(new Hole(a,F0,nu,true));
      std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      mySol = new AnaSol(new Lshape(a,F0,nu,true));
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		<<" nu="<< nu << std::endl;
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      mySol = new AnaSol(new CanTS(L,H,F0,true));
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		<<" F0="<< F0 << std::endl;
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::string primary, secondary, variables;
      int lines = atoi(strtok(NULL, " "));
      char* c = strtok(NULL, " ");
      if (c)
        code = atoi(c);
      else
        code = 0;
      STensorFunc* v=NULL;
      for (int i = 0; i < lines; i++) {
        std::string function = utl::readLine(is);
        size_t pos;
        if ((pos = function.find("Variables=")) != std::string::npos) {
          variables = function.substr(pos+10);
          if (variables[variables.size()-1] != ';')
            variables += ";";
        }
        if ((pos = function.find("Stress=")) != std::string::npos) {
          secondary = function.substr(pos+7);
          v = new EvalMultiFunction<STensorFunc,Vec3,SymmTensor>(secondary,
                                                                 6,variables);
        }
      }
      std::cout <<"\nAnalytical solution:" << std::endl;
      if (!variables.empty())
        std::cout << "\t Variables=" << variables << std::endl;
      if (v)
        std::cout << "\t Stress=" << secondary << std::endl;
      mySol = new AnaSol(v);
    }
    else
    {
      std::cerr <<"  ** SIMLinEl3D::parse: Unknown analytical solution "
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

      if (!mVec.empty())
	elp->setMaterial(mVec.front());
    }
  }

  else if (!strncasecmp(keyWord,"LOCAL_SYSTEM",12))
  {
    size_t i = 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    if (!strncasecmp(keyWord+i,"CYLINDRICZ",10))
      elp->setLocalSystem(new CylinderCS);
    else if (!strncasecmp(keyWord+i,"CYLINDER+SPHERE",15))
      elp->setLocalSystem(new CylinderSphereCS(atof(keyWord+i+15)));
    else
      std::cerr <<" *** SIMLinEl3D::parse: Unsupported coordinate system: "
		<< keyWord+i << std::endl;
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

  return true;
}


bool SIMLinEl3D::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd]);
  return true;
}


bool SIMLinEl3D::initNeumann (size_t propInd)
{
  TracFuncMap::const_iterator tit = myTracs.find(propInd);
  if (tit == myTracs.end()) return false;

  Elasticity* elp = dynamic_cast<Elasticity*>(myProblem);
  if (!elp) return false;

  elp->setTraction(tit->second);
  return true;
}
