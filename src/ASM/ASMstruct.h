// $Id$
//==============================================================================
//!
//! \file ASMstruct.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for structured spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_STRUCT_H
#define _ASM_STRUCT_H

#include "ASMbase.h"

namespace Go {
  class GeomObject;
}

/*!
  \brief Base class for structured spline-based FE assembly drivers.
  \details This class contains methods common for structured spline patches.
*/

class ASMstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMstruct(const ASMstruct& patch, unsigned char n_f);

public:
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == NULL; }

  //! \brief Resets the global element and node counters.
  static void resetNumbering(int n = 0);

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const = 0;

  //! \brief Returns the number of elements in each parameter direction.
  //! \param[out] n1 Number of elements in first (u) direction
  //! \param[out] n2 Number of elements in second (v) direction
  //! \param[out] n3 Number of elements in third (w) direction
  virtual bool getNoStructElms(int& n1, int& n2, int& n3) const = 0;

  using ASMbase::evalSolution;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integr Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const IntegrandBase& integr) const = 0;

protected:
  //! \brief Adds extraordinary nodes associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  //!
  //! \details This method is only a helper that is used by the addXElms methods
  //! of the dimension-specific sub-classes.
  bool addXNodes(unsigned short int dim, size_t nXn, IntVec& nodes);

protected:
  Go::GeomObject* geo; //!< Pointer to the actual spline geometry object

  //! Auxilliary node number map used when establishing Dirichlet constraints
  static std::map<int,int> xNode;

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter
};

#endif
