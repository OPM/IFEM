// $Id$
//==============================================================================
//!
//! \file BDF.h
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper functions for Backwards Difference based time stepping.
//!
//==============================================================================

#ifndef _BDF_H_
#define _BDF_H_

#include <vector>


namespace TimeIntegration //! Utilities for time integration.
{
  /*!
    \brief Helper class for Backwards Difference Formula integration schemes.
    \details Keeps track of coefficients, startup and extrapolation.
  */

  class BDF
  {
  public:
    //! \brief Default constructor.
    //! \param[in] order The order of the BDF scheme
    explicit BDF(int order = 0);
    //! \brief Empty destructor.
    virtual ~BDF() {}

    //! \brief Initializes the coefficients for the specified \a order.
    void setOrder(int order);

    //! \brief Returns the order to be used for current time step.
    int getOrder() const;
    //! \brief Returns the order of the scheme.
    int getActualOrder() const { return coefs.size() - degree; }
    //! \brief Returns the degree of the time derivative approximation.
    short int getDegree() const { return degree; }

    //! \brief Advances the time stepping scheme.
    bool advanceStep(double dt = 0.0, double dtn = 0.0);

    //! \brief Returns the BDF coefficients.
    virtual const std::vector<double>& getCoefs() const;

    //! \brief Indexing operator returning the idx'th coefficient.
    double operator[](int idx) const { return this->getCoefs()[idx]; }

    //! \brief Extrapolates values.
    //! \param[in] values The values to be extrapolated
    //! \return The extrapolated value
    template<class T> T extrapolate(const std::vector<T>& values) const
    {
      const T& v0 = values.front();
      if (step > 1 && this->getActualOrder() == 2) // second order
        return 2.0*v0 - values[1];
      else // first order
        return v0;
    }

  protected:
    short int           degree; //!< Degree of the time derivative approximation
    int                 step;   //!< Time step counter
    std::vector<double> coefs;  //!< The BDF coefficients
    std::vector<double> coefs1; //!< BDF coefficients for first time step
  };


  /*!
    \brief Helper class for Backwards Difference schemes for 2nd order problems.
  */

  class BDFD2 : public BDF
  {
  public:
    //! \brief Default constructor.
    //! \param[in] order The order of the BDF scheme
    //! \param[in] step_ Initial step position
    explicit BDFD2(int order = 2, int step_ = 0);

    //! \brief Returns the BDF coefficients.
    const std::vector<double>& getCoefs() const override;

  protected:
    std::vector<double> coefs2; //!< BDF coefficients for second time step
  };
}

#endif
