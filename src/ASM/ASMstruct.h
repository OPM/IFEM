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
  //! \brief The destructor frees the dynamically allocated spline objects.
  virtual ~ASMstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geomB == nullptr; }

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

  //! \brief Evaluates the basis functions at the specified point.
  //! \param[in] u First parameter value of evaluation point
  //! \param[in] v Second parameter value of evaluation point
  //! \param[in] w Third parameter value of evaluation point
  //! \param[out] N Basis function values
  virtual void evaluateBasis(double u, double v, double w, Vector& N) const = 0;

  using ASMbase::evalSolution;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integr Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const IntegrandBase& integr) const = 0;

  //! \brief Integrates a spatial dirac-delta function over a patch.
  //! \param integr Object with problem-specific data and methods
  //! \param glInt The integrated quantity
  //! \param[in] u Parameters of the non-zero point of the dirac-delta function
  //! \param[in] pval Function value at the specified point
  virtual bool diracPoint(Integrand& integr, GlobalIntegral& glInt,
                          const double* u, const Vec3& pval);

  //! \brief Checks if a separate projection basis is used for this patch.
  virtual bool separateProjectionBasis() const;

protected:
  //! \brief Adds extraordinary nodes associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  //!
  //! \details This method is only a helper that is used by the addXElms methods
  //! of the dimension-specific sub-classes.
  bool addXNodes(unsigned short int dim, size_t nXn, IntVec& nodes);

  //! \brief Performs a sanity check on the thread groups.
  //! \param[in] nodes The nodes to santiy check
  //! \param[in] group The group to check for
  //! \param[in] ignoreGlobalLM If \e true, ignore global lagrange multipliers
  //! \return \e true if the groups pass checks, otherwise \e false
  //!
  //! \details This checks that no nodes exist on several threads.
  bool checkThreadGroups(const std::vector<std::set<int>>& nodes,
                         int group, bool ignoreGlobalLM);

  //! \brief Computes the element border parameters.
  //! \param[in] iel 1-based element index
  //! \param[out] u Parameter values of the element borders
  virtual void getElementBorders(int iel, double* u) const = 0;

protected:
  Go::GeomObject* geomB; //!< Pointer to spline object of the geometry basis
  Go::GeomObject* projB; //!< Pointer to spline object of the projection basis
  Go::GeomObject* altProjB; //!< Pointer to spline object of the alternative projection basis
};

#endif
