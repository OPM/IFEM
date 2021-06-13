// $Id$
//==============================================================================
//!
//! \file SIMgeneric.h
//!
//! \date Aug 28 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generic SIM class with some added functionalities.
//!
//==============================================================================

#ifndef _SIM_GENERIC_H_
#define _SIM_GENERIC_H_

#include "SIMoutput.h"


/*!
  \brief Generic SIM class with some added functionalities.
  \details This class extends the SIMbase class with some added functionalities
  of generic character, which can be used to access the FE data and structures
  in a more flexible way.
*/

class SIMgeneric : public SIMoutput
{
protected:
  //! \brief Default constructor.
  explicit SIMgeneric(IntegrandBase* itg = nullptr) : SIMoutput(itg) {}

public:
  //! \brief Empty destructor.
  virtual ~SIMgeneric() {}

  //! \brief Creates a model with the default geometry (line, plane, cube).
  //! \return Pointer to the (first) spline patch of the model
  ASMbase* createDefaultModel();

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] par Parameters of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index containing the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, const double* par,
                     int deriv = 0, int patch = 1) const;

  //! \brief Evaluates the mapping of the geometry at the given point.
  //! \param[in] xi Dimensionless parameters in range [0,1] of the point
  //! \param[out] X The Cartesian coordinates of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[in] patch 1-based patch index containing the evaluation point
  //! \param[in] global If \e true, return global number, otherwise patch-local
  //! \return Patch-local or global node number of node that matches the point
  int evalPoint(const double* xi, Vec3& X, double* param = nullptr,
                int patch = 1, bool global = false) const;

  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameters of the point in the knot-span domain
  //! \param[in] patch 1-based patch index containing the point
  //! \param[in] global If \e true, return global number, otherwise patch-local
  //! \return Patch-local or global number of the element containing the point
  int findElementContaining(const double* param,
                            int patch = 1, bool global = false) const;

  //! \brief Calculates surface traction resultants.
  virtual bool calcBouForces(Vectors&, const Vectors&) { return false; }

  //! \brief Returns the norm index for the a VCP-recovered quantity.
  size_t getVCPindex(size_t idx = 1) const;
  //! \brief Returns the norm index for the integrated volume (3D) or area (2D).
  virtual size_t getVolumeIndex() const { return 0; }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] gNorm Global norm values
  //! \param[in] w Total number of characters in the norm labels
  virtual void printNorms(const Vectors& gNorm, size_t w = 36) const;

  //! \brief Prints a norm group to the log stream (app-specific).
  virtual void printNormGroup(const Vector&, const Vector&,
                              const std::string&) const {}

  //! \brief Returns the reference norm to base mesh adaptation upon.
  //! \param[in] gNorm Global norm values
  //! \param[in] adaptor 0-based norm group index to be used for mesh adaptation
  virtual double getReferenceNorm(const Vectors& gNorm, size_t adaptor) const;

  //! \brief Returns the global effectivity index.
  //! \param[in] gNorm Global norm values
  //! \param[in] idx 0-based norm group index
  //! \param[in] inorm 1-based norm index within the specified norm group
  virtual double getEffectivityIndex(const Vectors& gNorm,
                                     size_t idx, size_t inorm) const;

  //! \brief Interface for static condensation of the linear equation system.
  virtual bool staticCondensation(Matrix&, Vector&) { return false; }
  //! \brief Interface for recovery of internal DOFs from static condensation.
  virtual bool recoverInternals(const Vector&, Vector&) { return false; }
  //! \brief Interface for returning the superelement file name, if any.
  virtual std::string getSupelName() const { return ""; }
  //! \brief Interface for dumping supernode coordinates to file.
  virtual void dumpSupernodes(std::ostream&) const {}

protected:
  //! \brief Reverts the square-root operation on some norm quantities.
  //! \param gNorm Global norm values
  //! \param eNorm Element norm values
  bool revertSqrt(Vectors& gNorm, Matrix* eNorm);
};

#endif
