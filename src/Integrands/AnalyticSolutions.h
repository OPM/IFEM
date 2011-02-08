// $Id: AnalyticSolutions.h,v 1.10 2011-02-08 12:19:37 rho Exp $
//==============================================================================
//!
//! \file AnalyticSolutions.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for some linear elasticity and Poisson problems.
//!
//==============================================================================

#ifndef _ANALYTIC_SOLUTIONS_H
#define _ANALYTIC_SOLUTIONS_H

#include "Function.h"
#include "Tensor.h"


/*!
  \brief Analytic solution for an infinite plate with a hole.
*/

class Hole : public TensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  Hole(double r = 1.0, double f = 1.0, double P = 0.3, bool use3D = false)
    : a(r), F0(f), nu(P), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~Hole() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual Tensor evaluate(const Vec3& x) const;

private:
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
  double a;  //!< Hole radius
  double F0; //!< Load factor
  double nu; //!< Poisson's ratio
};


/*!
  \brief Analytic solution for the L-shaped domain.
*/

class Lshape : public TensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  Lshape(double r = 1.0, double f = 1.0, double P = 0.3, bool use3D = false);
  //! \brief Empty destructor.
  virtual ~Lshape() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual Tensor evaluate(const Vec3& x) const;

private:
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
  double a;  //!< Length parameter
  double F0; //!< Load factor
  double nu; //!< Poisson's ratio

  Tensor T;  //!< Local-to-global stress transformation tensor
};


/*!
  \brief Analytic solution for the cantilever beam with a tip shear load.
*/

class CanTS : public TensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CanTS(double l, double h, double f = 1.0, bool use3D = false)
    : L(l), H(h), F0(f), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~CanTS() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual Tensor evaluate(const Vec3& x) const;

private:
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
  double L;  //!< Length
  double H;  //!< Height
  double F0; //!< Load factor
};


/*!
  \brief Analytic solution for the cantilever beam with a tip moment load.
*/

class CanTM : public TensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CanTM(double h, double m = 1.0, bool use3D = false)
    : H(h), M0(m), is3D(use3D) {}
  //! \brief Empty destructor.
  virtual ~CanTM() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual Tensor evaluate(const Vec3& x) const;

private:
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
  double H;  //!< Height
  double M0; //!< Load factor
};


/*!
  \brief Analytic solution for the curved beam with end shear.
*/

class CurvedBeam : public TensorFunc
{
public:
  //! \brief Constructor with some default parameters.
  CurvedBeam(double u0 = 0.1, double Ri = 1.0, double Ro = 2.0,
	     double E = 2.1e7, bool use3D = false);
  //! \brief Empty destructor.
  virtual ~CurvedBeam() {}

protected:
  //! \brief Evaluates the analytic stress tensor at the point \a x.
  virtual Tensor evaluate(const Vec3& x) const;

private:
  bool is3D; //!< Flag telling whether to return a 3D stress tensor or not
  double a;  //!< Inner radius
  double b;  //!< Outer radius
  double PN; //!< Load parameter

  mutable Tensor T; //!< Local-to-global stress transformation tensor
};


/*!
  \brief Analytic solution for the Poisson equation on a square.
*/

class Square2D : public VecFunc
{
public:
  //! \brief Empty constructor.
  Square2D(double = 1.0) {}
  //! \brief Empty destructor.
  virtual ~Square2D() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Heat source for the Poisson equation on a square.
*/

class Square2DHeat : public RealFunc
{
public:
  //! \brief Empty constructor.
  Square2DHeat(double = 1.0) {}
  //! \brief Empty destructor.
  virtual ~Square2DHeat() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;
};


/*!
  \brief Analytic solution for the L-shape problem for the Poisson equation.
*/

class LshapePoisson : public VecFunc
{
public:
  //! \brief Empty constructor.
  LshapePoisson() {}
  //! \brief Empty destructor.
  virtual ~LshapePoisson() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Analytic solution for the Poisson equation on a cube.
*/

class PoissonCube : public VecFunc
{
public:
  //! \brief Empty Constructor.
  PoissonCube() {}
  //! \brief Empty destructor.
  virtual ~PoissonCube() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Heat source for the Poisson equation on a cube.
*/

class PoissonCubeSource : public RealFunc
{
public:
  //! \brief Empty constructor.
  PoissonCubeSource() {}
  //! \brief Empty destructor.
  virtual ~PoissonCubeSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;
};


/*!
  \brief Analytic solution for the Poisson equation on a line.
*/

class PoissonLine : public VecFunc
{
public:
  //! \brief Empty Constructor.
  PoissonLine(double r = 1.0) : L(r) {}
  //! \brief Empty destructor.
  virtual ~PoissonLine() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  double L; //!< Length parameter
};


/*!
  \brief Heat source for the Poisson equation on a line.
*/

class PoissonLineSource : public RealFunc
{
public:
  //! \brief Empty constructor.
  PoissonLineSource(double r = 1.0) : L(r) {}
  //! \brief Empty destructor.
  virtual ~PoissonLineSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double L; //!< Length parameter
};

#endif
