// $Id$
//==============================================================================
//!
//! \file BasisFunctionCache.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function cache.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_CACHE_H
#define _BASIS_FUNCTION_CACHE_H

#include "BasisFunctionVals.h"
#include "ASMenums.h"

#include <array>
#include <memory>

class Integrand;


/*!
  \brief Base class for basis function caches.
*/

template<size_t Dim> class BasisFunctionCache
{
public:
  //! \brief Default constructor.
  BasisFunctionCache();

  //! \brief Copy-constructor.
  BasisFunctionCache(const BasisFunctionCache& rhs);

  //! \brief Empty destructor.
  virtual ~BasisFunctionCache() = default;

  //! \brief Clears the basis function cache.
  void clear();

  //! \brief Initialize the basis function cache according to policy.
  //! \param nd Number of derivatives
  bool init(int nd);

  //! \brief Obtain basis function values/derivatives in an integration point.
  //! \param el Element of integration point (0-indexed)
  //! \param gp Integratin point on element (0-indexed)
  //! \param reduced If true, returns values for reduced integration scheme
  const BasisFunctionVals& getVals(size_t el, size_t gp, bool reduced = false);

  //! \brief Defines the integrand to be used.
  //! \param itg Integrand to use
  void setIntegrand(const Integrand* itg) { integrand = itg; }

  //! \brief Returns number of integration points.
  //! \param reduced True to return reduced quadrature
  const std::array<int,Dim>& nGauss(bool reduced = false)
  { return reduced ? reducedQ->ng : mainQ->ng; }

  //! \brief Return integration scheme weights.
  //! \param reduced True to return reduced quadrature
  const std::array<const double*,Dim>& weight(bool reduced = false) const
  { return reduced ? reducedQ->wg : mainQ->wg; }

  //! \brief Return integration scheme nodes.
  //! \param reduced True to return reduced quadrature
  const std::array<const double*,Dim>& coord(bool reduced = false) const
  { return reduced ? reducedQ->xg : mainQ->xg; }

  //! \brief Returns whether or not a reduced quadrature is enabled.
  bool hasReduced() const { return !reducedQ->gpar[0].empty(); }

  //! \brief Obtain a single integration point parameter.
  //! \param dir Direction of for integration point
  //! \param el Element number in given direction
  //! \param gp Integration point in given direction
  //! \param reduced True to return parameter for reduced quadrature
  virtual double getParam(int dir, size_t el, size_t gp, bool reduced = false) const;

  int basis = 1; //!< Basis to use

  //! \brief Called if application changes number of threads.
  void resizeThreadBuffers();

protected:
  //! \brief Template struct holding information about a quadrature.
  struct Quadrature {
    std::array<int,Dim> ng; //!< Number of integration point in each dimension
    std::array<const double*,Dim> xg; //!< Integration scheme nodes
    std::array<const double*,Dim> wg; //!< Integration scheme weights
    std::array<RealArray,Dim> gpar; //!< Global integration point parameters

    //! \brief Clears out configured quadrature.
    void reset()
    {
      ng.fill(0);
      xg.fill(nullptr);
      wg.fill(nullptr);
      gpar.fill({});
    }
  };

  //! \brief Class-specified initialization.
  virtual bool internalInit() = 0;

  //! \brief Calculates basis function info in a single integration point.
  //! \param el Element of integration point (0-indexed)
  //! \param gp Integratin point on element (0-indexed)
  //! \param reduced If true, returns values for reduced integration scheme
  virtual BasisFunctionVals calculatePt(size_t el, size_t gp,
                                        bool reduced = false) const = 0;

  //! \brief Calculates basis function info in all integration points.
  virtual void calculateAll() = 0;

  //! \brief Obtain global integration point index.
  //! \param el Element of integration point (0-indexed)
  //! \param gp Integration point on element (0-indexed)
  //! \param reduced If true return index for reduced integration scheme
  virtual size_t index (size_t el, size_t gp, bool reduced) const;

  //! \brief Obtain structured integration point indices.
  //! \param gp Global integration point
  //! \param reduced If true, returns values for reduced integration scheme
  std::array<size_t,Dim> gpIndex(size_t gp, bool reduced) const;

  std::vector<BasisFunctionVals> values; //!< Cache for main quadrature
  std::vector<BasisFunctionVals> valuesRed; //!< Cache for reduced quadrature
  const Integrand* integrand = nullptr; //!< Integrand to use
  int nderiv = 0; //!< Number of derivatives

  std::shared_ptr<Quadrature> mainQ; //!< Main quadrature information
  std::shared_ptr<Quadrature> reducedQ; //!< Reduced quadrature information

  size_t nTotal = 0; //!< Total number of main integration points
  size_t nTotalRed = 0; //!< Total number of reduced integration points
};

#endif
