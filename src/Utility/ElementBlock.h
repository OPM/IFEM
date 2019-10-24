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
#include <array>
#include <vector>


/*!
  \brief Class for storage of a standard FE grid block.
*/

class ElementBlock
{
public:
  //! \brief The constructor defines the number of nodes per element \a nenod.
  explicit ElementBlock(size_t nenod = 8);
  //! \brief Empty destructor.
  virtual ~ElementBlock() {}

  //! \brief Reallocates the internal arrays to fit a structured grid.
  //! \param[in] nI Number of nodes in I-direction
  //! \param[in] nJ Number of nodes in J-direction
  //! \param[in] nK Number of nodes in K-direction
  void resize(size_t nI, size_t nJ = 1, size_t nK = 1);

  //! \brief Reallocates the internal arrays to fit an unstructured grid.
  //! \param[in] nEl  Number of elements
  //! \param[in] nPts Number of nodes
  void unStructResize(size_t nEl, size_t nPts);

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

  //! \brief Adds a line element to the grid, assuming \a nen is equal to two.
  bool addLine(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2);

  //! \brief Assigns an external id to an element.
  void setElmId(size_t i, int iel) { MINEX[i-1] = iel; }
  //! \brief Returns the external id of an element.
  int getElmId(size_t i) const { return MINEX[i-1]; }

  //! \brief Returns the total number of nodes in the block.
  size_t getNoNodes() const { return coord.size(); }
  //! \brief Returns the total number of elements in the block.
  size_t getNoElms() const { return MINEX.size(); }
  //! \brief Returns the number of nodes per element.
  size_t getNoElmNodes() const { return nen; }
  //! \brief Sets the number of nodes per element.
  void setNoElmNodes(size_t nenod) { nen = nenod; }

  //! \brief Merges another element block into this one.
  void merge(const ElementBlock* other, std::vector<int>& nodeNums);
  //! \brief Merges another element block into this one.
  void merge(const ElementBlock& other);

  //! \brief Returns the beginning of the coordinate array.
  std::vector<Vec3>::const_iterator begin_XYZ() const { return coord.begin(); }
  //! \brief Returns the end of the coordinate array.
  std::vector<Vec3>::const_iterator end_XYZ() const { return coord.end(); }

  //! \brief Returns the coordinate of a given node.
  const Vec3& getCoord(size_t i) const { return coord[i]; }
  //! \brief Returns a pointer to the parameter values of a given node.
  const Real* getParam(size_t i) const { return param[i].data(); }

  //! \brief Returns a pointer to the element connectivity array.
  const int* getElements() const { return MMNPC.data(); }

  //! \brief Returns the coordinates of the center of the given elemment.
  utl::Point getCenter(size_t i) const;

private:
  typedef std::array<Real,3> Prm3; //!< Convenience type

  std::vector<Vec3> coord; //!< Vector of nodal coordinates
  std::vector<Prm3> param; //!< Vector of parameter values of the nodal points
  std::vector<int>  MMNPC; //!< Matrix of Matrices of Nodal Point Correspondance
  std::vector<int>  MINEX; //!< Matrix of Internal to External element numbers
  size_t            nen;   //!< Number of Element Nodes
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
  //! \brief Empty destructor.
  virtual ~CubeBlock() {}
};

#endif
