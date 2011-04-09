// $Id$
//==============================================================================
//!
//! \file ASMs2DmxLag.h
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D Lagrange mixed FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_MX_LAG_H
#define _ASM_S2D_MX_LAG_H

#include "ASMs2DLag.h"
#include "ASMmxBase.h"


/*!
  \brief Driver for assembly of structured 2D Lagrange mixed FE models.
  \details This class implements a two-field mixed formulation with Lagrangian
  basis functions. The geometry and the first field are of equal order and
  is one order higher than the second field.
*/

class ASMs2DmxLag : public ASMs2DLag, private ASMmxBase
{
public:
  //! \brief Constructor creating an instance by reading the given file.
  ASMs2DmxLag(const char* fName = 0, unsigned char n_s = 2,
	      unsigned char n_f1 = 2, unsigned char n_f2 = 1);
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs2DmxLag(std::istream& is, unsigned char n_s = 2,
	      unsigned char n_f1 = 2, unsigned char n_f2 = 1);
  //! \brief Default constructor creating an empty patch.
  ASMs2DmxLag(unsigned char n_s = 2,
	      unsigned char n_f1 = 2, unsigned char n_f2 = 1);
  //! \brief Empty destructor.
  virtual ~ASMs2DmxLag() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  virtual void clear();

  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis = 0) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int* sysMadof);

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  virtual bool connectPatch(int edge, ASMs2D& neighbor, int nedge,
			    bool revers = false);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  virtual void closeEdges(int dir, int = 0, int = 1);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());


  // Post-processing methods
  // =======================

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
			   const IntVec& nodes) const;

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const int*) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int*) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char = 0) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[in] basis Which basis to return size parameters for
  virtual bool getSize(int& n1, int& n2, int basis = 0) const;

private:
  size_t nx2; //!< Number of nodes in 1st parameter direction for second basis
  size_t ny2; //!< Number of nodes in 2nd parameter direction for second basis
};

#endif
