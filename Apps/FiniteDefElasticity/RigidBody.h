// $Id$
//==============================================================================
//!
//! \file RigidBody.h
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of rigid bodies in contact analysis.
//!
//==============================================================================

#ifndef _RIGID_BODY_H
#define _RIGID_BODY_H

#include "MatVec.h"
#include "Tensor.h"
#include "Vec3.h"
#include <map>

class ElementBlock;


/*!
  \brief Base class representing a general rigid body for contact analysis.
*/

class RigidBody
{
protected:
  //! \brief The constructor is protected to allow sub-classes only.
  //! \param[in] n Number of space dimensions
  //! \param[in] np Number of internal points
  RigidBody(unsigned char n, size_t np);

public:
  //! \brief Empty destructor.
  virtual ~RigidBody() {}

  //! \brief Initializes the internal points of the rigid body.
  void initPoints(const std::vector<Vec3>& p);
  //! \brief Initializes the global node numbers of the internal points.
  void initNodes(const std::vector<int>& nodes);
  //! \brief Renumbers the global node numbers of the internal points.
  void renumberNodes(const std::map<int,int>& old2new);

  //! \brief Creates a tesselated geometry of the body for visualization.
  virtual ElementBlock* tesselate() const = 0;

  //! \brief Prints out the rigid body definition to the given stream.
  virtual void print(std::ostream& os) const;

  //! \brief Updates the position of the rigid body.
  //! \param[in] displ Total displacement vector of the whole FE model
  virtual bool update(const RealArray& displ);

  //! \brief Returns the current position of the rigid body.
  Vec3 getPosition() const;

  //! \brief Evaluates the gap function at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] N Weight functions of the internal nodal points
  //! \return The gap between the rigid surface and a corresponding slave
  virtual double evalGap(const Vec3& X, Vec3& normal, RealArray& N) const = 0;

  //! \brief Evaluate the geometric stiffness at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] Tg Geometric stiffness matrix
  //! \param[out] N Weight functions of the internal nodal points
  virtual void geometricStiffness(const Vec3& X, Vec3& normal,
				  Matrix& Tg, RealArray& N) const = 0;

  //! \brief Returns the number of internal nodes.
  size_t getNoNodes() const { return MLGN.size(); }

  //! \brief Returns the global node number for the given internal point.
  //! \param[in] inod 1-based node index local to current body
  int getNodeID(size_t inod) const { return MLGN[inod-1]; }

protected:
  std::vector<Vec3> X0;   //!< Initial coordinates of the internal points
  std::vector<Vec3> Xn;   //!< Updated coordinates of the internal points
  std::vector<int>  MLGN; //!< Matrix of Local to Global Node numbers
  unsigned char     nsd;  //!< Number of space dimensions

public:
  double eps;    //!< Penalty parameter for the contact constraint of this body
  int    code;   //!< Property code associated with this contact body
  int    gBlock; //!< Geometry block ID on the VTF file for this contact body
};


/*!
  \brief Class representing a rigid sphere for contact analysis.
*/

class RigidSphere : public RigidBody
{
public:
  //! \brief Default constructor.
  RigidSphere(double r = 1.0) : RigidBody(3,1), R(r) {}
  //! \brief Empty destructor.
  virtual ~RigidSphere() {}

  //! \brief Prints out the rigid body definition to the given stream.
  virtual void print(std::ostream& os) const;

  //! \brief Creates a tesselated geometry of the body for visualization.
  virtual ElementBlock* tesselate() const;

  //! \brief Evaluates the gap function at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] N Weight functions of the internal nodal points
  //! \return The gap between the rigid surface and a corresponding slave
  virtual double evalGap(const Vec3& X, Vec3& normal, RealArray& N) const;

  //! \brief Evaluate the geometric stiffness at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] Tg Geometric stiffness matrix
  //! \param[out] N Weight functions of the internal nodal points
  virtual void geometricStiffness(const Vec3& X, Vec3& normal,
				  Matrix& Tg, RealArray& N) const;

private:
  double R; //!< Sphere radius
};


/*!
  \brief Class representing a rigid cylinder for contact analysis.
*/

class RigidCylinder : public RigidBody
{
public:
  //! \brief Default constructor.
  RigidCylinder(double r = 1.0, unsigned char n = 3);
  //! \brief Empty destructor.
  virtual ~RigidCylinder() {}

  //! \brief Prints out the rigid body definition to the given stream.
  virtual void print(std::ostream& os) const;

  //! \brief Creates a tesselated geometry of the body for visualization.
  virtual ElementBlock* tesselate() const;

  //! \brief Updates the position of the rigid body.
  //! \param[in] displ Total displacement vector of the whole FE model
  virtual bool update(const RealArray& displ);

  //! \brief Evaluates the gap function at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] N Weight functions of the internal nodal points
  //! \return The gap between the rigid surface and a corresponding slave
  virtual double evalGap(const Vec3& X, Vec3& normal, RealArray& N) const;

  //! \brief Evaluate the geometric stiffness at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] Tg Geometric stiffness matrix
  //! \param[out] N Weight functions of the internal nodal points
  virtual void geometricStiffness(const Vec3& X, Vec3& normal,
				  Matrix& Tg, RealArray& N) const;

private:
  Vec3 dirC; //!< Cylinder axis direction vector
  double R;  //!< Cylinder radius
  double L;  //!< Cylinder length
};


/*!
  \brief Class representing a rigid plane for contact analysis.
*/

class RigidPlane : public RigidBody
{
public:
  //! \brief Default constructor.
  RigidPlane(unsigned char n = 3) : RigidBody(n,n), Area(0.0) {}
  //! \brief Empty destructor.
  virtual ~RigidPlane() {}

  //! \brief Prints out the rigid body definition to the given stream.
  virtual void print(std::ostream& os) const;

  //! \brief Creates a tesselated geometry of the body for visualization.
  virtual ElementBlock* tesselate() const;

  //! \brief Updates the position of the rigid body.
  //! \param[in] displ Total displacement vector of the whole FE model
  virtual bool update(const RealArray& displ);

  //! \brief Evaluates the gap function at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] N Weight functions of the internal nodal points
  //! \return The gap between the rigid surface and a corresponding slave
  virtual double evalGap(const Vec3& X, Vec3& normal, RealArray& N) const;

  //! \brief Evaluate the geometric stiffness at an integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] normal Outward-direction normal vector of the rigid surface
  //! \param[out] Tg Geometric stiffness matrix
  //! \param[out] N Weight functions of the internal nodal points
  virtual void geometricStiffness(const Vec3& X, Vec3& normal,
				  Matrix& Tg, RealArray& N) const;

private:
  Vec3   nVec; //!< Plane normal vector
  double Area; //!< Area of triangle spanned by the three points
};

#endif
