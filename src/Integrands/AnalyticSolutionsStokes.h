// $Id$
//==============================================================================
//!
//! \file AnalyticSolutionsStokes.h
//!
//! \date Des 10 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Analytical solutions for some Stokes and Navier-Stokes problems.
//!
//==============================================================================

#ifndef _ANALYTIC_SOLUTIONS_STOKES_H
#define _ANALYTIC_SOLUTIONS_STOKES_H

#include "AnaSol.h"
#include "Tensor.h"


/*!
  \brief Analytic solution for Poiseuille flow.
*/

class Poiseuille : public AnaSol
{
 public:
  //! \brief Constructor with some default parameters.
  Poiseuille(double P = 1.0, double diam = 1.0, double len = 1.0,
	     double visc = 1.0);
  //! \brief Empty destructor.
  virtual ~Poiseuille() {}

 private:
  double D;   //!< Diameter of channel
  double L;   //!< Length of channel
  double mu;  //!< Fluid viscosity
  double Pin; //!< Inlet pressure value

  //! \brief Analytical pressure field.
  class Pressure : public RealFunc
  {
    const Poiseuille& params;

  protected:
    virtual double evaluate(const Vec3& x) const;

  public:
    Pressure(const Poiseuille& p) : params(p) {}
    virtual ~Pressure() {}
  };

  //! \brief Analytical velocity field.
  class Velocity : public VecFunc
  {
    const Poiseuille& params;

  protected:
    virtual Vec3 evaluate(const Vec3& x) const;

  public:
    Velocity(const Poiseuille& p) : params(p) {}
    virtual ~Velocity() {}
  };


  //! \brief Analytical pressure gradient field.
  class PressureGrad : public VecFunc
  {
    const Poiseuille& params;

  protected:
    virtual Vec3 evaluate(const Vec3& x) const;

  public:
   PressureGrad(const Poiseuille& p) : params(p) {}
    virtual ~PressureGrad() {}
  };

  //! \brief Analytical velocity gradient field.
  class VelocityGrad : public TensorFunc
  {
    const Poiseuille& params;

  protected:
    virtual Tensor evaluate(const Vec3& x) const;

  public:
    VelocityGrad(const Poiseuille& p) : params(p) {}
    virtual ~VelocityGrad() {}
  };
};


/*!
  \brief Artificial solution for Navier-Stokes equations.
*/

class TestSolution : public AnaSol
{
 public:
  //! \brief Constructor with some default parameters.
  //! \param[in] dens Fluid density
  //! \param[in] visc Fluid viscosity
  TestSolution(double dens = 1.0, double visc = 1.0);
  //! \brief Empty destructor.
  virtual ~TestSolution() {}

 private:
  double rho; //!< Fluid density
  double mu;  //!< Fluid viscosity

  //! \brief Analytical pressure field.
  class Pressure : public RealFunc
  {
  protected:
    virtual double evaluate(const Vec3& x) const;

  public:
    Pressure() {}
    virtual ~Pressure() {}
  };

  //! \brief Analytical velocity field.
  class Velocity : public VecFunc
  {
  protected:
    virtual Vec3 evaluate(const Vec3& x) const;

  public:
    Velocity() {}
    virtual ~Velocity() {}
  };

  //! \brief Analytical pressure gradient field.
  class PressureGrad : public VecFunc
  {
  protected:
    virtual Vec3 evaluate(const Vec3& x) const;

  public:
    PressureGrad() {}
    virtual ~PressureGrad() {}
  };

  //! \brief Analytical velocity gradient field.
  class VelocityGrad : public TensorFunc
  {
  protected:
    virtual Tensor evaluate(const Vec3& x) const;

  public:
    VelocityGrad() {}
    virtual ~VelocityGrad() {}
  };
};

#endif
