//==============================================================================
//!
//! \file Spalding.h
//!
//! \date Jan 25 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Implementation of Spalding parametrization of a turbulent
//!        boundary layer. Computes the mean velocity uplus parallel
//!        to a solid wall given the distance y and the tangential
//!        velocity component ut.
//!
//==============================================================================

#ifndef SPALDING_H
#define SPALDING_H

#include <cmath>


/*!
  \brief Class representing Spalding parametrization of a turbulent
         boundary layer
*/

class Spalding
{
 public:
  //! \brief The default constructor initializes parameters.
 Spalding(double eps = 1.0e-10, int mit = 400) : rtol(eps), maxit(mit) {}
  
  //! \brief Empty destructor
  virtual ~Spalding() {}

  // Computed tangent velocity
  bool computeTauB(double hb, double CbI, double nu, double ut, double& tauB) const;

 protected:
  //! \brief Function defining Spaling parametrization
  class f_Spalding 
  {
  public:
    //! \brief Constructor defining constants
    f_Spalding() : Chi(0.4), B(5.5) {}
    //! \brief Empty destructur
    virtual ~f_Spalding() {};

  protected:
    double Chi;       //! Constant in Spalding model
    double B;         //! Constant in Spalding model

    //! \brief Spaling function definition
    //! \param[in] up U plus parameter
    //! \param[in] y  Distance y from wall
    //! \param[in] ut Norm of tangential velocity
    //! \param[in] nu Viscosity
    virtual double evaluate(double up) const
    {  
      double val = up + exp(-Chi*B)*(exp(Chi*up)-1-Chi*up-0.5*pow(Chi*up,2.0)-pow(Chi*up,3.0)/6.0);

      return val;
    }

  public:
    //! \brief Operator returning the function value for the given argument.
    double operator()(double up) const 
    { return this->evaluate(up); }
  };

  //! \brief Function defining derivative of Spaling parametrization wrt. tauB
  class drdtauB_Spalding
  {
  public:
    //! \brief Constructor defining constants
    drdtauB_Spalding() : Chi(0.4), B(5.5) {}
    //! \brief Empty destructur
    virtual ~drdtauB_Spalding() {};

  protected:
    double Chi;       //! Constant in Spalding model
    double B;         //! Constant in Spalding model

    //! \brief Spaling residual definition
    //! \param[in] up    U plus parameter
    //! \param[in] tauB  Friction coefficient
    //! \param[in] ut Norm of tangential velocity
    //! \param[in] nu Viscosity
    //! \param[in] hb Grid size (y)
    //! \param[in] CbI Coefficient
    virtual double evaluate(double up, double ut, double nu, double hb, double CbI, double tauB) const
    {  
      double coeff = sqrt(ut/tauB);

      double val = 0.5*hb/(nu*CbI)*coeff;
      val += 0.5*(1.0 + Chi*exp(-Chi*B)*(exp(Chi*up) - 1.0 - Chi*up - 0.5*pow(Chi*up,2.0)))*coeff/tauB;

      return val;
    }

  public:
    //! \brief Operator returning the function value for the given argument.
    double operator()(double up, double ut, double nu, double hb, double CbI, double tauB) const 
    { return this->evaluate(up,ut,nu,hb,CbI,tauB); }
  };


  // Newton-Raphson parameters
  double rtol;                //!< Residual tolerance
  int    maxit;               //!< Maximal number of iterations
  
  // Function definition
  f_Spalding        f;              //!< Spalding parametrization
  drdtauB_Spalding  drdtauB;        //!< Residual of Spalding differentiated wrt. tauB
};

#endif
