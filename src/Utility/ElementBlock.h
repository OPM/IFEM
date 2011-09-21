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

#ifndef _ELEMENTBLOCK_H
#define _ELEMENTBLOCK_H

#include "Vec3.h"
#include <vector>

/*!
  \brief Class for storage of a standard FE grid block.
*/

class ElementBlock
{
public:
  //! The constructor defines the number of nodes per element \a nenod.
  ElementBlock(size_t nenod = 8);

  //! \brief Reallocates the internal arrays to fit a structured element block.
  //! \param[in] nI Number of element in I-direction
  //! \param[in] nJ Number of element in J-direction
  //! \param[in] nK Number of element in K-direction
  void resize(size_t nI, size_t nJ = 1, size_t nK = 1);

  //! \brief Reallocates the internal arrays to fit an unstructured element block.
  //! \param[in] nEl  Number of elements
  //! \param[in] nPts Number of nodes
  void unStructResize(size_t nEl, size_t nPts);

  //! \brief Defines the coordinates of node \a i
  bool setCoor(size_t i, real x, real y, real z);
  //! \brief Defines the \a j'th coordinate of node \a i
  bool setCoor(size_t i, size_t j, real x);

  //! \brief Defines the global number of element node \a i
  bool setNode(size_t i, int nodeNumb);

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

  //! \brief Merges another element block into this one.
  void merge(const ElementBlock* other, std::vector<int>& nodeNums);

  //! \brief Returns the beginning of the coord array.
  std::vector<Vec3>::const_iterator begin_XYZ() const { return coord.begin(); }
  //! \brief Returns the end of the coord array.
  std::vector<Vec3>::const_iterator end_XYZ() const { return coord.end(); }

  //! \brief Returns a pointer to the element connectivity array.
  const int* getElements() const { return &MMNPC.front(); }

private:
  std::vector<Vec3> coord; //!< Vector of nodal coordinates
  std::vector<int>  MMNPC; //!< Matrix of Matrices of Nodal Point Correspondance
  std::vector<int>  MINEX; //!< Matrix of Internal to External element numbers
  size_t            nen;   //!< Number of Element Nodes
};

#endif
