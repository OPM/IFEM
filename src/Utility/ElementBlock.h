// $Id$
//==============================================================================
//!
//! \file ElementBlock.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a standard FE grid block.
//!
//==============================================================================

#ifndef _ELEMENT_BLOCK_H
#define _ELEMENT_BLOCK_H

#include "Point.h"
#include "MatVec.h"
#include <array>


/*!
  \brief Class for storage of a standard FE grid block.
*/

class ElementBlock
{
public:
  //! \brief The constructor defines the number of nodes per element \a nenod.
  explicit ElementBlock(size_t nenod = 8);
  //! \brief Empty destructor.
  virtual ~ElementBlock() = default;

  //! \brief Reallocates the internal arrays to fit a structured grid.
  //! \param[in] nI Number of nodes in I-direction
  //! \param[in] nJ Number of nodes in J-direction
  //! \param[in] nK Number of nodes in K-direction
  void resize(size_t nI, size_t nJ = 1, size_t nK = 1);

  //! \brief Reallocates the internal arrays to fit an unstructured grid.
  //! \param[in] nEl Number of elements
  //! \param[in] nPts Number of nodes
  //! \param[in] nMNPC Total size of element connectivity array
  void unStructResize(size_t nEl, size_t nPts, size_t nMNPC = 0);

  //! \brief Defines the \a j'th coordinate of node \a i.
  bool setCoor(size_t i, size_t j, Real x);
  //! \brief Defines the coordinates of node \a i.
  bool setCoor(size_t i, const Vec3& X);
  //! \brief Defines the coordinates of node \a i.
  bool setCoor(size_t i, Real x, Real y, Real z)
  { return this->setCoor(i,Vec3(x,y,z)); }
  //! \brief Defines the parameter values for node \a i.
  bool setParams(size_t i, Real u, Real v, Real w = Real(0));

  //! \brief Defines the global number of element node \a i.
  bool setNode(size_t i, int nodeNumb);
  //! \brief Marks the end of current element for unstructured grids.
  bool endOfElm(size_t& i);

  //! \brief Adds a line element to the grid, assuming \ref nen is equal to two.
  bool addLine(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2);
  //! \brief Adds a line element to the grid, assuming \ref nen is equal to two.
  //! \param[in] i1 Index of existing node to use as start point
  //! \param[in] X2 Coordinates of new node to use as end point
  //! \param[in] elmId External element ID (generate if negative)
  size_t addLine(size_t i1, const Vec3& X2, int elmId = -1);

  //! \brief Assigns an external id to an element.
  void setElmId(size_t i, int iel) { MINEX[i-1] = iel; }
  //! \brief Returns the external id of an element.
  int getElmId(size_t i) const { return MINEX[i-1]; }
  //! \brief Returns the internal index of an element in case of mixed types.
  size_t getElmIndex(size_t i) const { return i < elmIdx.size() ? elmIdx[i]:i; }

  //! \brief Returns the total number of nodes in the block.
  size_t getNoNodes() const { return coord.size(); }
  //! \brief Returns the total number of elements in the block.
  size_t getNoElms() const { return MINEX.size(); }
  //! \brief Returns the number of nodes per element.
  size_t getNoElmNodes() const { return nen; }
  //! \brief Sets the number of nodes per element.
  void setNoElmNodes(size_t nenod) { nen = nenod; }

  //! \brief Merges another element block into this one.
  void merge(const ElementBlock* other,
             std::vector<int>& nodeNums, bool uniqNodes = true);
  //! \brief Merges another element block into this one.
  void merge(const ElementBlock& other, bool uniqNodes = true);

  //! \brief Applies a transformation matrix from local to global system.
  void transform(const Matrix& Tlg);

  //! \brief Returns the beginning of the coordinate array.
  std::vector<Vec3>::const_iterator begin_XYZ() const { return coord.begin(); }
  //! \brief Returns the end of the coordinate array.
  std::vector<Vec3>::const_iterator end_XYZ() const { return coord.end(); }

  //! \brief Returns the coordinate of a given node.
  const Vec3& getCoord(size_t i) const { return coord[i]; }
  //! \brief Returns a pointer to the parameter values of a given node.
  const Real* getParam(size_t i) const;

  //! \brief Returns a pointer to the element connectivity array.
  const int* getElements() const { return MMNPC.data(); }
  //! \brief Get element connectivity array for elements with \a nenod nodes.
  bool getElements(std::vector<int>& mnpc, size_t nenod) const;

  //! \brief Returns the coordinates of the center of the given elemment.
  utl::Point getCenter(size_t i) const;

protected:
  using Prm3 = std::array<Real,3>; //!< Convenience type

  std::vector<Vec3> coord; //!< Vector of nodal coordinates
  std::vector<Prm3> param; //!< Vector of parameter values of the nodal points
  std::vector<int>  MMNPC; //!< Matrix of Matrices of Nodal Point Correspondance
  std::vector<int>  MINEX; //!< Matrix of Internal to External element numbers
  size_t            nen;   //!< Number of Element Nodes

private:
  std::vector<size_t> elmIdx; //!< Internal element order in case of mixed types

public:
  static double eps; //!< Element shrinkage factor
};


/*!
  \brief Class for single-element plane geometries.
  \details This class is used to create visualization of rigid contact planes.
*/

class PlaneBlock : public ElementBlock
{
public:
  //! \brief Constructor defining a plane from three points.
  //! \param[in] X0 First corner
  //! \param[in] X1 Second corner
  //! \param[in] X2 Third corner
  PlaneBlock(const Vec3& X0, const Vec3& X1, const Vec3& X2);
};


/*!
  \brief Class for single-element cube geometries.
  \details This class is used to create small cube shapes in the model,
  mainly for visualization of result points, etc.
*/

class CubeBlock : public ElementBlock
{
public:
  //! \brief The constructor defines a cube centred at specified point.
  //! \param[in] X0 Center of the cube
  //! \param[in] dX Length of each cube edge
  CubeBlock(const Vec3& X0, double dX);
};


/*!
  \brief Class for sphere geometries.
  \details This class is used to create sphere shapes in the model,
  for visualization of single-point elements, rigid contact objects, etc.
*/

class SphereBlock : public ElementBlock
{
public:
  //! \brief The constructor defines a sphere centred at specified point.
  //! \param[in] X0 Center of the sphere
  //! \param[in] R Sphere diameter
  //! \param[in] nTheta Number of elements around equator
  //! \param[in] nPhi Number of elements from pole to pole
  SphereBlock(const Vec3& X0, double R, size_t nTheta = 180, size_t nPhi = 60);
};


/*!
  \brief Class for cylinder geometries.
  \details This class is used to create cylinder shapes in the model,
  for visualization of rigid contect objects, etc.
*/

class CylinderBlock : public ElementBlock
{
public:
  //! \brief The constructor defines a sphere centred at specified point.
  //! \param[in] X0 First end point
  //! \param[in] X1 Second end point
  //! \param[in] R Cylinder diameter
  //! \param[in] nTheta Number of elements in circular direction
  CylinderBlock(const Vec3& X0, const Vec3& X1,
                double R, size_t nTheta = 180);
};

#endif
