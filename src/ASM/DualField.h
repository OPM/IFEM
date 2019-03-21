// $Id$
//==============================================================================
//!
//! \file DualField.h
//!
//! \date Oct 22 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of dual fields for goal-oriented error estimation.
//!
//==============================================================================

#ifndef _DUAL_FIELD_H
#define _DUAL_FIELD_H

#include "Function.h"


/*!
  \brief Class representing a dual field for goal-oriented error estimation.
*/

class DualRealFunc : public RealFunc
{
public:
  //! \brief Constructor for 2D problems (in XY-plane).
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  DualRealFunc(const Vec3& o, const Vec3& n,
               double d = 1.0, double w = 0.0);
  //! \brief Constructor for 3D problems.
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] XZp Point in the local XZ-plane
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  DualRealFunc(const Vec3& o, const Vec3& n, const Vec3& XZp,
               double d = 1.0, double w = 0.0);
  //! \brief Empty destructor.
  virtual ~DualRealFunc() {}

  //! \brief Returns the offset from the cross section origin \ b X0.
  double ecc(const Vec3& X, int i) const { return X(i) - X0(i); }
  //! \brief Returns the function value at the given point.
  //! \param[in] X Spatial evaluation point
  //! \param[in] ignoreDomain If \e true, the function is evaluated without
  //! considering the boundaries of the function domain. Otherwise, it is
  //! assumed identically zero outside these boundaries.
  double value(const Vec3& X, bool ignoreDomain = false) const;
  //! \brief Checks if the point \b X is within the function domain.
  virtual bool inDomain(const Vec3& X) const;

protected:
  //! \brief Evaluates the dual field function.
  virtual double evaluate(const Vec3& X) const { return this->value(X); }

private:
  Vec3      X0; //!< Global coordinates of the cross section origin
  Vec3  normal; //!< Outward-directed normal vector of the cross section
  Vec3 tangent; //!< Vector defining the local y-direction of the cross section
  double depth; //!< Depth of the the dual function domain
  double width; //!< Width of the the dual function domain (0=infinite)
};


/*!
  \brief Class representing a dual vector field.
*/

class DualVecFunc : public VecFunc
{
  int comp; //!< Which section force component to extract for [1,6]

public:
  //! \brief Constructor for 2D problems (in XY-plane).
  //! \param[in] c Sectional force component index to do extraction for
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  DualVecFunc(int c, const Vec3& o, const Vec3& n,
              double d = 1.0, double w = 0.0);
  //! \brief Constructor for 3D problems.
  //! \param[in] c Sectional force component index to do extraction for
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] XZp Point in the local XZ-plane
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  DualVecFunc(int c, const Vec3& o, const Vec3& n, const Vec3& XZp,
              double d = 1.0, double w = 0.0);
  //! \brief Empty destructor.
  virtual ~DualVecFunc() {}

  //! \brief Checks if the point \b X is within the function domain.
  virtual bool inDomain(const Vec3& X) const { return W.inDomain(X); }

protected:
  //! \brief Evaluates the dual field function.
  //! \param[in] X The spatial point to evaluate the function at
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  DualRealFunc W; //!< The scalar dual field
};

#endif
