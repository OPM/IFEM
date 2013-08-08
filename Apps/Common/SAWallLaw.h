//==============================================================================
//!
//! \file SAWallLaw.h
//!
//! \date Jun 19 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Implementation of Spalding parametrization of a turbulent
//!        boundary layer. Computes the mean velocity uplus parallel
//!        to a solid wall given the distance y and the tangential
//!        velocity component ut.
//!
//==============================================================================

#ifndef SA_WALL_LAW_H
#define SA_WALL_LAW_H

#include <cmath>
#include <iostream>

/*!
  \brief Class representing Spalart-Allmaras parametrization of a turbulent
         boundary layer
*/

class SAWallLaw 
{
 protected:
  // Function defining Spalart-Allmaras parametrization
  class f_SA_wall_law
  {
  protected:
    //! Model constants
    const double B;
    const double a1;
    const double a2;
    const double b1;
    const double b2;
    const double c1;
    const double c2;
    const double c3;
    const double c4;
    
  public:
    //! \brief Constructor defining constants
    f_SA_wall_law(); 
    //! \brief Empty destructur
    virtual ~f_SA_wall_law() {};
    
    //! \brief Spalart-Allmaras function definition
    //! \param[in] yp y plus parameter
    double up(double yp) const
    {
      double yppa1 = yp + a1;
      double yppa2 = yp + a2;
      
      double up =  B + c1*log(pow(yppa1,2.0) + pow(b1,2.0));
      up -= c2*log(pow(yppa2,2.0) + pow(b2,2.0));
      up -= c3*atan2(b1,yppa1);
      up -= c4*atan2(b2,yppa2);  
      
      return up;
    }
    
    double dupdyp(double yp) const
    {
      double yppa1 = yp + a1;
      double yppa2 = yp + a2;

      double value  = 2.0*c1*yppa1/(pow(yppa1,2.0) + pow(b1,2.0));
      value -= 2.0*c2*yppa2/(pow(yppa2,2.0) + pow(b2,2.0));
      value -= c3/(yppa1*(1.0+pow(b1/yppa1,2.0)));
      value -= c4/(yppa2*(1.0+pow(b2/yppa2,2.0)));
      
      return value;
    }
  };
  
 protected:
  // Newton-Raphson parameters
  double rtol;                // Residual tolerance
  int    maxit;               // Maximal number of iterations
  
  // Parameters
  const double K;             // von Karman constant
  
  // Function definition
  f_SA_wall_law SA;           // Spalart-Allmaras parametrization

 public:
  //! \brief The default constructor initializes parameters.
 SAWallLaw(double eps = 1.0e-10, int mit = 400) : rtol(eps), maxit(mit), K(0.41) {}
  
  //! \brief Empty destructor
  virtual ~SAWallLaw() {}

  bool computeYplus(double y, double nu, double utan, double& yplus) const;
  bool computeUstar(double y, double nu, double utan, double& ustar) const;
  bool computeTauB(double y, double nu, double utan, double& tauB) const;
  bool computeNu(double y, double nu, double utan, double& nuwall) const;

 
};    

#endif
