// $Id$
//==============================================================================
//!
//! \file SAMpatchPara.h
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for multi-patch models.
//!
//==============================================================================

#ifndef _SAM_PATCH_PARA_H
#define _SAM_PATCH_PARA_H

#include "SAMpatch.h"
#ifdef PARALLEL_PETSC
#include "petscksp.h"
#endif


/*!
  \brief This is a sub-class of SAMpatch with support for parallel processing.
*/

class SAMpatchPara : public SAMpatch
{
public:
  //! \brief The constructor initializes the \a l2gn array.
  SAMpatchPara(const IntVec& l2gn_mp);
  //! \brief The destructor destroys the index set arrays.
  virtual ~SAMpatchPara();

  //! \brief Returns the number of equations (free DOFs) on this processor.
  virtual int getNoEquations() const { return nleq; }

  //! \brief Computes number of couplings for each local dof
  //! in a distributed matrix.
  //! \param[in] ifirst Global number of first local DOF
  //! \param[in] ilast Global number of last local DOF pluss 1
  //! \param[out] d_nnz Number of diagonal couplings for each local DOF
  //! \param[out] o_nnz Number of off-diagonal couplings for each local DOF
  //! \return \e false if number of couplings is not computed, otherwise \e true
  virtual bool getNoDofCouplings(int ifirst, int ilast,
				 IntVec& d_nnz, IntVec& o_nnz) const;

  //! \brief Initializes the system load vector prior to the element assembly.
  //! \param sysRHS The system right-hand-side load vector to be initialized
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e false if no free DOFs in the system, otherwise \e true
  virtual bool initForAssembly(SystemVector& sysRHS,
			       Vector* reactionForces = 0) const;

  //! \brief Adds element stiffness contributions to the system load vector.
  //! \param sysRHS  The right-hand-side system load vector
  //! \param[in] eK  The element stiffness matrix
  //! \param[in] iel Identifier for the element that \a eK belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are added into the right-hand-side system load vector.
  virtual bool assembleSystem(SystemVector& sysRHS,
			      const Matrix& eK, int iel = 0,
			      Vector* reactionForces = 0) const;

  //! \brief Adds an element load vector into the system load vector.
  //! \param sysRHS  The right-hand-side system load vector
  //! \param[in] eS  The element load vector
  //! \param[in] iel Identifier for the element that \a eS belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assembleSystem(SystemVector& sysRHS,
			      const RealArray& eS, int iel = 0,
			      Vector* reactionForces = 0) const;

  //! \brief Finds the matrix of equation numbers for an element.
  //! \param[out] meen Matrix of element equation numbers
  //! \param[in] iel Identifier for the element to get the equation numbers for
  //! \param[in] nedof Number of degrees of freedom in the element
  //! (used for internal consistency checking, unless zero)
  virtual bool getElmEqns(IntVec& meen, int iel, int nedof = 0) const;

  //! \brief Updates the multi-point constraint array \a TTCC.
  virtual bool updateConstraintEqs(const std::vector<ASMbase*>& model,
				   const Vector* prevSol = 0);

  //! \brief Expands a solution vector from equation-ordering to DOF-ordering.
  //! \param[in] solVec Solution vector, length = NEQ
  //! \param[out] displ Displacement vector, length = NDOF = 3*NNOD
  //! \return \e false if the length of \a solVec is invalid, otherwise \e true
  //!
  //! \details The size of the solution vector that comes out of the linear
  //! equation solver equals the number of free DOFs in the system (=NEQ).
  //! That is, all fixed or constrained (slave) DOFs are not present.
  //! Before we can compute derived element quantities we therefore need to
  //! extract the resulting displacement values also for the constrained DOFs.
  virtual bool expandSolution(const SystemVector& solVec, Vector& displ) const;

  //! \brief Computes the dot-product of two vectors.
  //! \param[in] x First vector in dot-product
  //! \param[in] y Second vector in dot-product
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of dot-product
  //!
  //! \details Both vectors are defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual real dot(const Vector& x, const Vector& y, char dofType = 'D') const;

  //! \brief Computes the L2-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of L2-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual real normL2(const Vector& x, char dofType = 'D') const;

  //! \brief Computes the inf-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of inf-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual real normInf(const Vector& x, size_t& comp, char dofType = 'D') const;

protected:
  //! \brief Initializes the multi-point constraint arrays
  //! \a MPMCEQ, \a MMCEQ and \a TTCC.
  virtual bool initConstraintEqs(const std::vector<ASMbase*>& model);

  //! \brief Initializes the DOF-to-equation connectivity array \a MEQN.
  virtual bool initSystemEquations();

private:
  // Parameters for parallel computing
  int    nProc;      //!< Number of processes
  int    nleq;       //!< Number of equations for this processor
  int    nnodGlob;   //!< Number of global nodes;
  IntVec ghostNodes; //!< Indices for ghost nodes
  IntVec l2gn;       //!< Local-to-global node numbers for this processor
#ifdef PARALLEL_PETSC
  IS     iglob;      //!< Index set for global numbering
  IS     iloc;       //!< Index set for local numbering
#endif
};

#endif
