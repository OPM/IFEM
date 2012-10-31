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

#include "MatVec.h"
/*!
  \brief Helper class for BDF methods.
  \details Keeps track of coefficients, extrapolation and startup
!*/

class BDF {
  public:
    //! \brief Constructor
    //! \param[in] order The order of the BDF scheme (1,2 currently)
    BDF(int order);

    //! \brief Empty destructor
    virtual ~BDF() {}

    //! \brief Return order to be used for current timestep
    int getOrder() const { return step < 2?1:coefs.size()-1; }

    //! \brief Return order of the scheme
    int getActualOrder() const { return coefs.size()-1; }

    //! \brief Advance time scheme
    void advanceStep() { step++; }

    //! \brief Returns the BDF coefficients
    const std::vector<double>& getCoefs() const;

    double operator[](int idx) const
    {
      return getCoefs()[idx];
    }

    //! \brief Extrapolate value
    //! \param[in] values The values to extrapolate based on
    //! \return The extrapolated value
    double extrapolate(const double* values) const;
  protected:
    std::vector<double> coefs;
    std::vector<double> coefs1;
    int step;
};

#endif
