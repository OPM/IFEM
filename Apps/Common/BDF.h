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

#ifndef BDF_H_
#define BDF_H_

#include <vector>


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
    //! \param[in] order The order of the BDF scheme (1,2 currently)
    BDF(int order);

    //! \brief Returns order to be used for current time step.
    int getOrder() const { return step < 2 ? 1 : coefs.size()-1; }

    //! \brief Returns order of the schemes.
    int getActualOrder() const { return coefs.size()-1; }

    //! \brief Set the order
    void setOrder(int order);

    //! \brief Advances the time stepping scheme.
    void advanceStep() { step++; }

    //! \brief Advances the time stepping scheme.
    void advanceStep(double dt, double dtn);

    //! \brief Returns the BDF coefficients.
    const std::vector<double>& getCoefs() const;

    //! \brief Indexing operator returning the idx'th coefficient.
    double operator[](int idx) const { return this->getCoefs()[idx]; }

    //! \brief Extrapolates values.
    //! \param[in] values The values to extrapolate based on
    //! \return The extrapolated value
    double extrapolate(const double* values) const;

  protected:
    std::vector<double> coefs;  //!< The BDF coefficients
    std::vector<double> coefs1; //!< BDF coefficients for first time step
    int                 step;   //!< Time step counter
  };
}

#endif
