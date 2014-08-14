// $Id$
//==============================================================================
//!
//! \file BDF.h
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper functions for BDF based time stepping.
//!
//==============================================================================

#ifndef _BDF_H_
#define _BDF_H_

#include <vector>
#include <iostream>

namespace TimeIntegration //! Utilities for time integration.
{
  /*!
    \brief Helper class for BDF methods.
    \details Keeps track of coefficients, extrapolation and startup.
  */

  class BDF
  {
  public:
    //! \brief The constructor initializes the coefficients.
    //! \param[in] order The order of the BDF scheme
    //! \param[in] step_ Initial step position
    BDF(int order = 0, int step_ = 0);
    
    //! \brief Destructor
    virtual ~BDF() {}
    
    //! \brief Sets the order of the scheme.
    virtual void setOrder(int order);

    //! \brief Returns the order to be used for current time step.
    virtual int getOrder() const { return step < 2 ? coefs1.size()-1 : coefs.size()-1; }
    //! \brief Returns the order of the scheme.
    virtual int getActualOrder() const { return coefs.size()-1; }
    //! \brief Returns the degree for the time derivative approximation
    virtual int getDegree() const { return 1; }

    //! \brief Advances the time stepping scheme.
    virtual void advanceStep() { ++step; }
    //! \brief Advances the time stepping scheme.
    virtual void advanceStep(double dt, double dtn);

    //! \brief Returns the BDF coefficients.
    virtual const std::vector<double>& getCoefs() const;

    //! \brief Indexing operator returning the idx'th coefficient.
    double operator[](int idx) const { return this->getCoefs()[idx]; }

    //! \brief Extrapolates values.
    //! \param[in] values The values to extrapolate
    //! \return The extrapolated value
    virtual double extrapolate(const double* values) const;

  protected:
    std::vector<double> coefs;  //!< The BDF coefficients
    std::vector<double> coefs1; //!< BDF coefficients for first time step
    int                 step;   //!< Time step counter
  };

  /*!
    \brief Helper class for BDF methods for 2nd order problems
    \details Keeps track of coefficients, extrapolation and startup.
  */
  class BDFD2
  {
  public:
    //! \brief The constructor initializes the coefficients.
    //! \param[in] order The order of the BDF scheme
    BDFD2(int order = 2);

    //! \brief Destructor
    virtual ~BDFD2() {}
    
    //! \brief Sets the order of the scheme.
    void setOrder(int order);

    //! \brief Returns the order to be used for current time step.
    int getOrder() const  { return step < 3 ? 1 : coefs.size()-2; }
    //! \brief Returns the order of the scheme.
    int getActualOrder() const { return coefs.size()-2; }
    //! \brief Returns the degree for the time derivative approximation
    virtual int getDegree() const { return 2; }

    //! \brief Advances the time stepping scheme.
    void advanceStep(double dt, double dtn) { ++step; } 

    //! \brief Returns the BDF coefficients.
    const std::vector<double>& getCoefs() const;

    //! \brief Indexing operator returning the idx'th coefficient.
    double operator[](int idx) const { return this->getCoefs()[idx]; }

    //! \brief Extrapolates values.
    //! \param[in] values The values to extrapolate
    //! \return The extrapolated value
    double extrapolate(const double* values) const;

  protected:
    std::vector<double> coefs;  //!< The BDF coefficients
    std::vector<double> coefs1; //!< BDF coefficients for first time step
    std::vector<double> coefs2; //!< BDF coefficients for second time step
    int                 step;   //!< Time step counter
  };
}

#endif
