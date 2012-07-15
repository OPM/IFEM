// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.C
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for finite deformation FEM analysis.
//!
//==============================================================================

#include "NonlinearDriver.h"
#include "SIMbase.h"


bool NonlinearDriver::solutionNorms (const TimeDomain& time, bool energyNorm,
				     double zero_tol, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  size_t nsd = model->getNoSpaceDim();
  size_t iMax[nsd];
  double dMax[nsd];
  double normL2 = model->solutionNorms(solution.front(),dMax,iMax);
  RealArray RF;
  bool haveReac = model->getCurrentReactions(RF,solution.front());

  Vectors gNorm;
  if (energyNorm)
  {
    model->setMode(SIM::RECOVERY);
    model->setQuadratureRule(model->opt.nGauss[1]);
    if (!model->solutionNorms(time,solution,gNorm))
      gNorm.clear();
  }

  if (myPid > 0) return true;

  std::streamsize stdPrec = outPrec > 0 ? std::cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tol;

  std::cout <<"  Primary solution summary: L2-norm            : "
	    << utl::trunc(normL2);

  char D = 'X';
  for (size_t d = 0; d < nsd; d++, D++)
    if (utl::trunc(dMax[d]) != 0.0)
      std::cout <<"\n                            Max " << D
		<<"-displacement : "<< dMax[d] <<" node "<< iMax[d];

  if (haveReac)
  {
    std::cout <<"\n  Total reaction forces: Sum(R) =";
    for (size_t i = 1; i < RF.size(); i++)
      std::cout <<" "<< utl::trunc(RF[i]);
    if (utl::trunc(RF.front()) != 0.0)
      std::cout <<"\n  displacement*reactions: (R,u) = "<< RF.front();
  }

  if (!gNorm.empty())
  {
    const Vector& norm = gNorm.front();
    if (norm.size() > 0)
    {
      std::cout <<"\n  Energy norm:    |u^h| = a(u^h,u^h)^0.5 : "
                << utl::trunc(norm(1));
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t a(u^h,u^h) = "<< utl::trunc(norm(1)*norm(1));
      std::cout.precision(oldPrec);
    }
    if (norm.size() > 1 && utl::trunc(norm(2)) != 0.0)
    {
      std::cout <<"\n  External energy: ((f,u^h)+(t,u^h))^0.5 : "<< norm(2);
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t(f,u)+(t,u) = "<< norm(2)*norm(2);
      std::cout.precision(oldPrec);
    }
    if (norm.size() > 2)
      std::cout <<"\n  Stress norm, L2: (sigma^h,sigma^h)^0.5 : "<< norm(3);
    if (norm.size() > 3)
      std::cout <<"\n  Pressure norm, L2:       (p^h,p^h)^0.5 : "<< norm(4)
                <<"\t(p^h = trace(sigma^h)/3)";
    if (norm.size() > 4)
      std::cout <<"\n  Deviatoric stress norm:  (s^d,s^d)^0.5 : "<< norm(5)
                <<"\t(s^d = sigma^h - p^h*I)";
    if (norm.size() > 5)
      std::cout <<"\n  Stress norm, von Mises: vm(sigma^h)    : "<< norm(6);
    std::cout << std::endl;
  }

  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) std::cout.precision(stdPrec);
  return true;
}
