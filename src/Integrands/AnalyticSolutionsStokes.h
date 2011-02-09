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
  \brief Analytic solution for Poiseuille flow
*/

class Poiseuille : public AnaSol
{
 public:
  //! \brief Constructor for pressure inlet BC
  Poiseuille(double P = 1.0, double diam = 1.0, double len = 1.0, double visc = 1.0);
  //! \brief Destructor
  virtual ~Poiseuille();

 private:
  double D;    // Diameter of channel
  double L;    // Length of channel
  double mu;   // Fluid viscosity
  double Pin;  // Inlet value

  // Analytical pressure
  class Pressure : public RealFunc
  {
  private:
    Poiseuille& params;

  protected:
    real evaluate(const Vec3& x) const;

  public:
    Pressure(Poiseuille& p) : params(p) {}
    virtual ~Pressure() {}
  };

  // Analytical velocity
  class Velocity : public VecFunc
  {
  private:
    Poiseuille& params;

  protected:
    Vec3 evaluate(const Vec3& x) const;

  public:
    Velocity(Poiseuille& p) : params(p) {}
    virtual ~Velocity() {}
  };


  // Analytical pressure gradient
  class PressureGrad : public VecFunc
  {
  private:
    Poiseuille& params;

  protected:
    Vec3 evaluate(const Vec3& x) const;

  public:
   PressureGrad(Poiseuille& p) : params(p) {}
    virtual ~PressureGrad() {}
  };

  // Analytical velocity gradient
  class VelocityGrad : public TensorFunc
  {
  private:
    Poiseuille& params;

  protected:
    Tensor evaluate(const Vec3& x) const;

  public:
    VelocityGrad(Poiseuille& p) : params(p) {}
    virtual ~VelocityGrad() {}
  };
};


/*!
  \brief Artificial solution for Navier-Stokes equations
*/

class TestSolution : public AnaSol
{
 public:
  //! \brief Constructor 
  //! \param[in] rho Fluid density
  //! \param[in] mu  Fluid viscosity
  TestSolution(double rho = 1.0, double mu = 1.0);
  //! \brief Destructor
  virtual ~TestSolution();

 private:
  double rho;  // Fluid density
  double mu;   // Fluid viscosity

  // Analytical pressure
  class Pressure : public RealFunc
  {
  protected:
    real evaluate(const Vec3& x) const;

  public:
    Pressure() {}
    virtual ~Pressure() {}
  };

  // Analytical velocity
  class Velocity : public VecFunc
  {
  protected:
    Vec3 evaluate(const Vec3& x) const;

  public:
    Velocity() {}
    virtual ~Velocity() {}
  };


  // Analytical pressure gradient
  class PressureGrad : public VecFunc
  {
  protected:
    Vec3 evaluate(const Vec3& x) const;

  public:
    PressureGrad() {}
    virtual ~PressureGrad() {}
  };

  // Analytical velocity gradient
  class VelocityGrad : public TensorFunc
  {
  protected:
    Tensor evaluate(const Vec3& x) const;

  public:
    VelocityGrad() {}
    virtual ~VelocityGrad() {}
  };
};

#endif
