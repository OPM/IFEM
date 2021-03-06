// $Id$
//==============================================================================
//!
//! \file SAMpatchPETSc.h
//!
//! \date Aug 31 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of global norms for distributed models.
//!
//==============================================================================

#ifndef _SAM_PATCH_PETSC_H
#define _SAM_PATCH_PETSC_H

#include "SAMpatch.h"
#include "PETScSupport.h"
#include <map>

class ProcessAdm;


/*!
  \brief This is a sub-class of SAMpatch with support for parallel processing.
*/

class SAMpatchPETSc : public SAMpatch
{
public:
  //! \brief The constructor initializes the \a l2gn array.
  //! \param[in] g2ln Global-to-local node numbers for this processor
  //! \param[in] padm Parallel process administrator
  SAMpatchPETSc(const std::map<int,int>& g2ln, const ProcessAdm& padm);
  //! \brief The destructor destroys the index set arrays.
  virtual ~SAMpatchPETSc();

  //! \brief Computes the dot-product of two vectors.
  //! \param[in] x First vector in dot-product
  //! \param[in] y Second vector in dot-product
  //! \param[in] dofType Only consider DOFs of this type
  //! \return Value of dot-product
  //!
  //! \details Both vectors are defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real dot(const Vector& x, const Vector& y, char dofType = 'D') const;

  //! \brief Computes the L2-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param[in] dofType Only consider DOFs of this type
  //! \return Value of L2-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normL2(const Vector& x, char dofType = 'D') const;

  //! \brief Computes the inf-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider DOFs of this type
  //! \return Value of inf-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normInf(const Vector& x, size_t& comp, char dofType = 'D') const;

  //! \brief Expands a solution vector from equation-ordering to DOF-ordering.
  //! \param[in] solVec Solution vector, length = NEQ
  //! \param[out] dofVec Degrees of freedom vector, length = NDOF
  //! \param[in] scaleSD Scaling factor for specified (slave) DOFs
  //! \return \e false if the length of \a solVec is invalid, otherwise \e true
  //!
  //! \details The size of the solution vector that comes out of the linear
  //! equation solver equals the number of free DOFs in the system (=NEQ).
  //! That is, all fixed or constrained (slave) DOFs are not present.
  //! Before we can compute derived element quantities we therefore need to
  //! extract the resulting solution values also for the constrained DOFs.
  virtual bool expandSolution(const SystemVector& solVec,
                              Vector& dofVec, Real scaleSD) const;

private:
  // Parameters for parallel computing
  int    nProc;      //!< Number of processes
  IntVec ghostNodes; //!< Indices for the ghost nodes
  IntVec l2gn;       //!< Local-to-global node numbers for this processor

  const ProcessAdm& adm; //!< Parallel process administrator

#ifdef HAVE_MPI
  mutable IS glob2LocEq = nullptr; //!< Index set for global-to-local equations

  //! \brief Struct holding dof vector info
  struct DofIS {
    IS local;        //!< Local index set for dof type
    IS global;       //!< Global index set for dof type
    int nDofs;       //!< Number of DOFs on this process
    VecScatter ctx;  //!< Scatterer
    bool scatterCreated = false; //!< True if the scatterer is created
  };
  mutable std::map<char,DofIS> dofIS; //!< Map of dof type scatter info

  //! \brief Returns a parallel index set for a given dofType.
  DofIS& getIS(char dofType) const;
#endif
};

#endif
