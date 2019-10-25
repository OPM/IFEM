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
#include "Point.h"
#include <set>

class ASMbase;


/*!
  \brief Class representing a dual field for goal-oriented error estimation.
*/

class DualRealFunc : public RealFunc
{
public:
  //! \brief Constructor for 3D problems.
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] XZp Point in the local XZ-plane
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  //! \param[in] p The affected patch (null=all)
  DualRealFunc(const Vec3& o, const Vec3& n, const Vec3& XZp,
               double d, double w = 0.0, ASMbase* p = nullptr);
  //! \brief Constructor for 2D problems (in XY-plane).
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  //! \param[in] p The affected patch (null=all)
  DualRealFunc(const Vec3& o, const Vec3& n,
               double d, double w = 0.0, ASMbase* p = nullptr);
  //! \brief Constructor for point extraction.
  //! \param[in] o Point to extract the point quantity at
  //! \param[in] d Lower-left and upper-rigth corner of function domain
  //! \param[in] p The affected patch (null=all)
  //! \param[in] eps Parameter tolerance for searching elements around a point
  DualRealFunc(const utl::Point& o, const Vec3Pair& d,
               ASMbase* p = nullptr, double eps = 0.0);
  //! \brief Empty destructor.
  virtual ~DualRealFunc() {}

  //! \brief Returns the local X (normal) direction of the cross section.
  const Vec3& x() const { return normal; }
  //! \brief Returns the local Y (tangent) direction of the cross section.
  const Vec3& y() const { return tangent; }
  //! \brief Returns the local Z (2nd tangent) direction of the cross section.
  Vec3 z() const { return Vec3(normal,tangent); }

  //! \brief Returns the offset from the cross section origin \ b X0.
  double ecc(const Vec3& X, int i) const { return X(i) - X0(i); }
  //! \brief Returns the function value at the given point.
  //! \param[in] X Spatial evaluation point
  //! \param[in] ignoreDomain If \e true, the function is evaluated without
  //! considering the boundaries of the function domain. Otherwise, it is
  //! assumed identically zero outside these boundaries.
  double value(const Vec3& X, bool ignoreDomain = false) const;
  //! \brief Returns \e true, if the dual field is for extracting point values.
  bool isPointExtraction() const { return depth <= 0.0; }
  //! \brief Checks if the point \b X is within the function domain.
  virtual bool inDomain(const Vec3& X) const;
  //! \brief Returns \e true if current patch is affected by this function.
  virtual bool initPatch(size_t idx);

protected:
  //! \brief Evaluates the dual field function.
  virtual double evaluate(const Vec3& X) const { return this->value(X); }

private:
  utl::Point X0; //!< Coordinates of cross section origin or reference point

  Vec3  normal; //!< Outward-directed normal vector of the cross section
  Vec3 tangent; //!< Vector defining the local y-direction of the cross section
  double depth; //!< Depth of the the dual function domain (0=point extraction)
  double width; //!< Width of the the dual function domain (0=infinite)
  Vec3     Xll; //!< Lower-left corner of box domain for point extraction
  Vec3     Xur; //!< Upper-right corner of box domain for point extraction

  ASMbase*      patch; //!< The patch that this field is defined on
  std::set<int> Delem; //!< List of elements in domain
};


/*!
  \brief Class representing a dual vector field.
*/

class DualVecFunc : public VecFunc
{
  int comp; //!< Which section force component to extract for [1,6]

public:
  //! \brief Constructor for 3D problems.
  //! \param[in] c Sectional force component index to do extraction for
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] XZp Point in the local XZ-plane
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  //! \param[in] p The affected patch (null=all)
  DualVecFunc(int c, const Vec3& o, const Vec3& n, const Vec3& XZp,
              double d, double w = 0.0, ASMbase* p = nullptr);
  //! \brief Constructor for 2D problems (in XY-plane).
  //! \param[in] c Sectional force component index to do extraction for
  //! \param[in] o Origin of local cross section coordinate system
  //! \param[in] n Cross section normal
  //! \param[in] d Depth of dual function domain
  //! \param[in] w Width of dual function domain (0=infinite)
  //! \param[in] p The affected patch (null=all)
  DualVecFunc(int c, const Vec3& o, const Vec3& n,
              double d, double w = 0.0, ASMbase* p = nullptr);
  //! \brief Constructor for stress extraction.
  //! \param[in] c Stress component index to do extraction for
  //! \param[in] o Point to extract the point quantity at
  //! \param[in] d Lower-left and upper-rigth corner of function domain
  //! \param[in] p The affected patch (null=all)
  //! \param[in] eps Parameter tolerance for searching elements around a point
  DualVecFunc(int c, const utl::Point& o, const Vec3Pair& d,
              ASMbase* p = nullptr, double eps = 0.0);
  //! \brief Empty destructor.
  virtual ~DualVecFunc() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return W.isPointExtraction() ? 3:2; }

  //! \brief Checks if the point \b X is within the function domain.
  virtual bool inDomain(const Vec3& X) const { return W.inDomain(X); }
  //! \brief Returns \e true if current patch is affected by this function.
  virtual bool initPatch(size_t idx) { return W.initPatch(idx); }

  //! \brief Returns a const reference to the scalar dual field.
  const RealFunc& getW() const { return W; }

  //! \brief Returns a representative scalar equivalent of the function value.
  virtual double getScalarValue(const Vec3& X) const { return W.value(X); }

protected:
  //! \brief Evaluates the dual field function.
  //! \param[in] X The spatial point to evaluate the function at
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  DualRealFunc W; //!< The scalar dual field
};

#endif
