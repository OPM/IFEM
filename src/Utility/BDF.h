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


namespace TimeIntegration //! Utilities for time integration.
{
  /*!
    \brief Helper class for BDF schemes.
    \details Keeps track of coefficients, startup and extrapolation.
  */

  class BDF
  {
  public:
    //! \brief Default constructor.
    //! \param[in] order The order of the BDF scheme
    BDF(int order = 0) : step(0), coefs1(1,1.0) { this->setOrder(order); }
    //! \brief Empty destructor.
    virtual ~BDF() {}

    //! \brief Initializes the coefficients for the specified \a order.
    virtual void setOrder(int order);

    //! \brief Returns the order to be used for current time step.
    int getOrder() const;
    //! \brief Returns the order of the scheme.
    int getActualOrder() const { return coefs.size() - this->getDegree(); }

    //! \brief Returns the degree for the time derivative approximation.
    virtual int getDegree() const { return 1; }

    //! \brief Advances the time stepping scheme.
    virtual void advanceStep(double dt = 0.0, double dtn = 0.0);

    //! \brief Returns the BDF coefficients.
    virtual const std::vector<double>& getCoefs() const;

    //! \brief Indexing operator returning the idx'th coefficient.
    double operator[](int idx) const { return this->getCoefs()[idx]; }

    //! \brief Extrapolates values.
    //! \param[in] values The values to be extrapolated
    //! \return The extrapolated value
    template<class T>
    T extrapolate(const std::vector<T>& values) const
    {
      if (step > 1 && this->getActualOrder() == 2) // second order
        return 2.0*values[0] - values[1];
      else // first order
        return values[0];
    }

  protected:
    int                 step;   //!< Time step counter
    std::vector<double> coefs;  //!< The BDF coefficients
    std::vector<double> coefs1; //!< BDF coefficients for first time step
  };

  /*!
    \brief Helper class for BDF schemes for 2nd order problems.
  */

  class BDFD2 : public BDF
  {
  public:
    //! \brief Default constructor.
    //! \param[in] order The order of the BDF scheme
    //! \param[in] step_ Initial step position
    BDFD2(int order = 2, int step_ = 0) { this->setOrder(order); step = step_; }
    //! \brief Empty destructor.
    virtual ~BDFD2() {}

    //! \brief Initializes the coefficients for the specified \a order.
    virtual void setOrder(int order);

    //! \brief Returns the degree for the time derivative approximation.
    virtual int getDegree() const { return 2; }

    //! \brief Advances the time stepping scheme.
    virtual void advanceStep(double = 0.0, double = 0.0) { ++step; }

    //! \brief Returns the BDF coefficients.
    virtual const std::vector<double>& getCoefs() const;

  protected:
    std::vector<double> coefs2; //!< BDF coefficients for second time step
  };
}

#endif
