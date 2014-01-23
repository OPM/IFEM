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
  //! \brief Function defining Spalart-Allmaras parametrization
  class f_SA_wall_law
  {
  protected:
    const double B;  //!< Model constant
    const double a1; //!< Model constant
    const double a2; //!< Model constant
    const double b1; //!< Model constant
    const double b2; //!< Model constant
    const double c1; //!< Model constant
    const double c2; //!< Model constant
    const double c3; //!< Model constant
    const double c4; //!< Model constant
    
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
    
    //! \brief Spalart-Allmaras derivative definition
    //! \param[in] yp y plus parameter
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
  double rtol;                //!< Residual tolerance
  int    maxit;               //!< Maximal number of iterations
  
  // Parameters
  const double K;             //!< von Karman constant
  
  // Function definition
  f_SA_wall_law SA;           //!< Spalart-Allmaras parametrization

 public:
  //! \brief The default constructor initializes parameters.
  //! \param eps Relative tolerance in Newton-Raphson loop
  //! \param mit Maximum number of Newton iterations
 SAWallLaw(double eps = 1.0e-10, int mit = 400) : rtol(eps), maxit(mit), K(0.41) {}
  
  //! \brief Empty destructor
  virtual ~SAWallLaw() {}

  //! \brief Compute the y-plus value
  //! \param y The y (wall distance) value
  //! \param nu The viscosity
  //! \param utan The tangential velocity
  //! \param[out] yplus The yplus value
  //! \return True if value was obtained
  bool computeYplus(double y, double nu, double utan, double& yplus) const;

  //! \brief Compute the u-star value
  //! \param y The y (wall distance) value
  //! \param nu The viscosity
  //! \param utan The tangential velocity
  //! \param[out] ustar The ustar value
  //! \return True if value was obtained
  bool computeUstar(double y, double nu, double utan, double& ustar) const;

  //! \brief Compute the tau-b value
  //! \param y The y (wall distance) value
  //! \param nu The viscosity
  //! \param utan The tangential velocity
  //! \param[out] tauB The tau-b value
  //! \return True if value was obtained
  bool computeTauB(double y, double nu, double utan, double& tauB) const;

  //! \brief Compute the wall-viscosity value
  //! \param y The y (wall distance) value
  //! \param nu The viscosity
  //! \param utan The tangential velocity
  //! \param[out] nuwall The nuwall value
  //! \return True if value was obtained
  bool computeNu(double y, double nu, double utan, double& nuwall) const;
};    

#endif
