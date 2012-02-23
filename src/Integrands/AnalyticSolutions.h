// $Id$
//==============================================================================
//!
//! \file AnalyticSolutions.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for linear elasticity problems.
//!
//==============================================================================

#ifndef _ANALYTIC_SOLUTIONS_H
#define _ANALYTIC_SOLUTIONS_H

#include "Function.h"
#include "Tensor.h"


/*!
  \brief Analytic solution for an infinite plate with a hole.
*/

class Hole : public STensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  Hole(double r = 1.0, double f = 1.0, double P = 0.3, bool use3D = false)
    : a(r), F0(f), nu(P), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~Hole() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

private:
  double a;  //!< Hole radius
  double F0; //!< Load factor
  double nu; //!< Poisson's ratio
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
};


/*!
  \brief Analytic solution for the L-shaped domain.
*/

class Lshape : public STensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  Lshape(double r = 1.0, double f = 1.0, double P = 0.3, bool use3D = false);
  //! \brief Empty destructor.
  virtual ~Lshape() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

private:
  double a;  //!< Length parameter
  double F0; //!< Load factor
  double nu; //!< Poisson's ratio
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not

  Tensor T;  //!< Local-to-global stress transformation tensor
};


/*!
  \brief Analytic solution for the cantilever beam with a tip shear load.
*/

class CanTS : public STensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CanTS(double l, double h, double f = 1.0, bool use3D = false)
    : L(l), H(h), F0(f), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~CanTS() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

private:
  double L;  //!< Length
  double H;  //!< Height
  double F0; //!< Load factor
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
};


/*!
  \brief Analytic solution for the cantilever beam with a tip moment load.
*/

class CanTM : public STensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CanTM(double h, double m = 1.0, bool use3D = false)
    : H(h), M0(m), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~CanTM() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

private:
  double H;  //!< Height
  double M0; //!< Load factor
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
};


/*!
  \brief Analytic solution for the curved beam with end shear.
*/

class CurvedBeam : public STensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CurvedBeam(double u0 = 0.1, double Ri = 1.0, double Ro = 2.0,
	     double E = 2.1e7, bool use3D = false);
  //! \brief Empty destructor.
  virtual ~CurvedBeam() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

private:
  double a;  //!< Inner radius
  double b;  //!< Outer radius
  double PN; //!< Load parameter
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
};


/*!
  \brief Analytic solution for the simply supported rectangular thin plate.
*/

class NavierPlate : public STensorFunc
{
public:
  //! \brief Constructor for plate with constant pressure load.
  NavierPlate(double a, double b, double t, double E, double Poiss, double P);
  //! \brief Constructor for plate with partial pressure or point load.
  NavierPlate(double a, double b, double t, double E, double Poiss, double P,
	      double xi_, double eta_, double c = 0.0, double d = 0.0);
  //! \brief Empty destructor.
  virtual ~NavierPlate() {}

protected:
  //! \brief Evaluates the analytic stress resultant tensor at the point \a x.
  virtual SymmTensor evaluate(const Vec3& x) const;

  //! \brief Adds the m'th and n'th terms of the plate solution to the moments.
  void addTerms(std::vector<double>& M, double x, double y, int m, int n) const;

private:
  double alpha; //!< pi/(plate length)
  double beta;  //!< pi/(plate width)
  double pz;    //!< Load parameter
  double nu;    //!< Poisson's ratio

  char   type; //!< Load type parameter
  double xi;   //!< X-position of point/partial load
  double eta;  //!< Y-position of point/partial load
  double c2;   //!< Partial load extension in X-direction
  double d2;   //!< Partial load extension in Y-direction
  int    inc;  //!< Increment in Fourier term summation (1 or 2)
};

#endif
