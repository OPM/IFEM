// $Id$
//==============================================================================
//!
//! \file IntegrandBase.C
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base classes representing FEM integrands.
//!
//==============================================================================

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Utilities.h"


IntegrandBase::~IntegrandBase ()
{
  if (myMats) delete myMats;
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC,
				 const Vec3&, size_t)
{
  return this->initElement(MNPC);
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC)
{
  if (myMats)
    myMats->withLHS = true;

  // Extract the primary solution vectors for this element
  int ierr = 0;
  for (size_t i = 0; i < mySols.size() && i < primsol.size() && ierr == 0; i++)
    if (!primsol[i].empty())
      ierr = utl::gather(MNPC,npv,primsol[i],mySols[i]);

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool IntegrandBase::initElement (const std::vector<int>&,
				 const std::vector<int>&, size_t)
{
  std::cerr <<" *** IntegrandBase::initElement must be reimplemented"
	    <<" for mixed problems."<< std::endl;
  return false;
}


bool IntegrandBase::initElementBou (const std::vector<int>& MNPC)
{
  if (myMats)
    myMats->withLHS = false;

  // Extract current primary solution vector for this element
  int ierr = 0;
  if (!mySols.empty() && !primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),mySols.front());

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElementBou: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool IntegrandBase::initElementBou (const std::vector<int>&,
				    const std::vector<int>&, size_t)
{
  std::cerr <<" *** IntegrandBase::initElementBou must be reimplemented"
	    <<" for mixed problems."<< std::endl;
  return false;
}


bool Integrand::evalSol (Vector&,
			 const Vector&, const Matrix&,
			 const Vec3&, const std::vector<int>&) const
{
  std::cerr <<" *** Integrand::evalSol not implemented."<< std::endl;
  return false;
}


bool Integrand::evalSol (Vector& s,
			 const Vector& N1, const Vector&,
			 const Matrix& dN1dX, const Matrix&, const Vec3& X,
			 const std::vector<int>& MNPC1,
			 const std::vector<int>&) const
{
  return this->evalSol(s,N1,dN1dX,X,MNPC1);
}


bool IntegrandBase::evalPrimSol (Vector&, const VecFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalPrimSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const TensorFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const STensorFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalPrimSol (double&, const RealFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalPrimSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const VecFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


void NormBase::initIntegration (const TimeDomain& time)
{
  myProblem.initIntegration(time);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
			    const Vec3& Xc, size_t nPt)
{
  return myProblem.initElement(MNPC,Xc,nPt);
}


bool NormBase::initElement (const std::vector<int>& MNPC)
{
  return myProblem.initElement(MNPC);
}


bool NormBase::initElement (const std::vector<int>& MNPC1,
			    const std::vector<int>& MNPC2, size_t n1)
{
  return myProblem.initElement(MNPC1,MNPC2,n1);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC)
{
  return myProblem.initElementBou(MNPC);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC1,
			       const std::vector<int>& MNPC2, size_t n1)
{
  return myProblem.initElementBou(MNPC1,MNPC2,n1);
}


ElmNorm& NormBase::getElmNormBuffer (LocalIntegral*& elmInt, const size_t nn)
{
  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (eNorm) return *eNorm;

  static double* data = 0;
  if (!data) data = new double[nn];
  static ElmNorm buf(data);
  memset(data,0,nn*sizeof(double));
  elmInt = &buf;
  return buf;
}
