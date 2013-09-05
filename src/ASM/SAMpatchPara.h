// $Id$
//==============================================================================
//!
//! \file SAMpatchPara.h
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for distributed models.
//!
//==============================================================================

#ifndef _SAM_PATCH_PARA_H
#define _SAM_PATCH_PARA_H

#include "SAMpatch.h"
#ifdef HAS_PETSC
#include "petscksp.h"
#include "petscsys.h"
#else
typedef int PetscInt; //!< To avoid compilation failures
typedef double PetscReal; //!< To avoid compilation failures
#endif
#include <map>

typedef std::vector<PetscInt>    PetscIntVec;  //!< PETSc integer vector
typedef std::vector<PetscReal>   PetscRealVec; //!< PETSc real vector
typedef std::vector<PetscIntVec> PetscIntMat;  //!< PETSc integer matrix


/*!
  \brief This is a sub-class of SAMpatch with support for parallel processing.
*/

class SAMpatchPara : public SAMpatch
{
  typedef std::vector<IntVec> IntMat; //!< General integer matrix

public:
  //! \brief The constructor initializes the \a l2gn array.
  //! \param[in] g2ln Global-to-local node numbers for this processor
  SAMpatchPara(const std::map<int,int>& g2ln);
  //! \brief The destructor destroys the index set arrays.
  virtual ~SAMpatchPara();

  //! \brief Allocates the dynamic arrays and populates them with data.
  //! \param[in] model All spline patches in the model
  //! \param[in] numNod Total number of unique nodes in the model
  virtual bool init(const std::vector<ASMbase*>& model, int numNod = 0);

  //! \brief Returns the number of equations (free DOFs) on this processor.
  virtual int getNoEquations() const { return nleq; }
  //! \brief Returns min global equation number for this process.
  virtual int getMinEqNumber() const { return ieqmin; }
  //! \brief Returns max global equation number for this process.
  virtual int getMaxEqNumber() const { return ieqmax; }

  //! \brief Computes number of couplings for each local dof
  //! in a distributed matrix.
  //! \param[in] ifirst Global number of first local DOF
  //! \param[in] ilast Global number of last local DOF pluss 1
  //! \param[out] d_nnz Number of diagonal couplings for each local DOF
  //! \param[out] o_nnz Number of off-diagonal couplings for each local DOF
  //! \return \e false if number of couplings is not computed, otherwise \e true
  virtual bool getNoDofCouplings(int ifirst, int ilast,
				 IntVec& d_nnz, IntVec& o_nnz) const;

  //! \brief Computes number of couplings for each local dof
  //! in a distributed block matrix.
  //! \param[in] ifirst First equation number
  //! \param[in] ilast  Last equation number
  //! \param[in] ncomps Number of nodal components per block
  //! \param[out] d_nnz Number of diagonal couplings for each local DOF
  //! \param[out] o_nnz Number of off-diagonal couplings for each local DOF
  //! \return \e false if number of couplings is not computed, otherwise \e true
  virtual bool getNoDofCouplings(int ifirst, int ilast,
				 const IntVec& ncomps,
				 std::vector<IntMat>& d_nnz,
				 std::vector<IntMat>& o_nnz) const;

  //! \brief Initializes the system load vector prior to the element assembly.
  //! \param sysRHS The system right-hand-side load vector to be initialized
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e false if no free DOFs in the system, otherwise \e true
  virtual bool initForAssembly(SystemVector& sysRHS,
			       Vector* reactionForces = NULL) const;

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
			      Vector* reactionForces = NULL) const;

  //! \brief Adds an element load vector into the system load vector.
  //! \param sysRHS  The right-hand-side system load vector
  //! \param[in] eS  The element load vector
  //! \param[in] iel Identifier for the element that \a eS belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assembleSystem(SystemVector& sysRHS,
			      const RealArray& eS, int iel = 0,
			      Vector* reactionForces = NULL) const;

  //! \brief Finds the matrix of equation numbers for an element.
  //! \param[out] meen Matrix of element equation numbers
  //! \param[in] iel Identifier for the element to get the equation numbers for
  //! \param[in] nedof Number of degrees of freedom in the element
  //! (used for internal consistency checking, unless zero)
  virtual bool getElmEqns(IntVec& meen, int iel, int nedof = 0) const;
  //! \brief Finds the matrix of equation numbers for an element.
  //! \param[out] meen Matrix of element equation numbers
  //! \param[in] iel Identifier for the element to get the equation numbers for
  //! \param[in] f1 The first field to extract equation numbers for
  //! \param[in] f2 The last field to extract equation numbers for
  //! \param[in] globalEq If \e true, extract global equation numbers
  //! \param[in] nedof Number of degrees of freedom in the element
  //! (used for internal consistency checking, unless zero)
  virtual bool getElmEqns(IntVec& meen, int iel, int f1, int f2,
			  bool globalEq = false, int nedof = 0) const;

  //! \brief Expands a solution vector from equation-ordering to DOF-ordering.
  //! \param[in] solVec Solution vector, length = NEQ
  //! \param[out] displ Displacement vector, length = NDOF = 3*NNOD
  //! \param[in] scaleSD Scaling factor for specified (slave) DOFs
  //! \return \e false if the length of \a solVec is invalid, otherwise \e true
  virtual bool expandSolution(const SystemVector& solVec, Vector& displ,
			      Real scaleSD = 1.0) const;

  //! \brief Computes the dot-product of two vectors.
  //! \param[in] x First vector in dot-product
  //! \param[in] y Second vector in dot-product
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of dot-product
  //!
  //! \details Both vectors are defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real dot(const Vector& x, const Vector& y, char dofType = 'D') const;

  //! \brief Computes the L2-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of L2-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normL2(const Vector& x, char dofType = 'D') const;

  //! \brief Computes the inf-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of inf-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normInf(const Vector& x, size_t& comp, char dofType = 'D') const;

  //! \brief Gets the mapping for ordering the local DOFs in line ordering
  //!  Assumes initial XYZ ordering.
  //! \param[out] order Index set with mapping from old to new ordering
  //! \param[in]  dir   Direction of line ordering (ex. 213)
  //! \param[in]  nf    Number of fields
  bool getDirOrdering(PetscIntVec& order, int dir, int nf = 0) const;

  //! \brief Splits the local dofs into a number of unique subdomains.
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  bool getLocalSubdomains(PetscIntMat& locSubds,
			  int nx = 1, int ny = 1, int nz = 1) const;

  //! \brief Splits the local dofs into a number of unique subdomains.
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in]  overlap Overlap between subdomains
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  bool getSubdomains(PetscIntMat& locSubds, int overlap = 1,
		     int nx = 1, int ny = 1, int nz = 1) const;

  //! \brief Splits the local dofs corresponding to given fields
  //! into a number of unique subdomains.
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in]  f1, f2 First and last field to extract dofs from
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  bool getLocalSubdomainsBlock(PetscIntMat& locSubds, int f1 = 0, int f2 = 0,
			       int nx = 1, int ny = 1, int nz = 1) const;

  //! \brief Splits the local dofs corresponding to given fields
  //! into a number of unique subdomains.
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in]  f1, f2 First and last field to extract dofs from
  //! \param[in]  overlap Overlap between subdomains
  //! \param[in]  nx, ny, nz Number of paritionings in each direction
  bool getSubdomainsBlock(PetscIntMat& locSubds, int f1 = 0, int f2 = 0,
			  int overlap = 1, int nx = 1, int ny = 1, int nz = 1) const;

  //! \brief Extracts the node coordinates of all local nodes.
  //! \param[out] coords Coordinates of local nodes stored as [x0,y0,z0,x1,y1,z1,...]
  bool getLocalNodeCoordinates(PetscRealVec& coords) const;

  //! \brief Get number of space dimensions
  size_t getNoSpaceDim() const;

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
  IntVec ghostNodes; //!< Indices for the ghost nodes
  IntVec l2gn;       //!< Local-to-global node numbers for this processor
  int    ieqmin;     //!< Minium equation number
  int    ieqmax;     //!< Maximun equation number
#ifdef PARALLEL_PETSC
  IS     iglob;      //!< Index set for global numbering
  IS     iloc;       //!< Index set for local numbering
#endif

  // For domain decomposition preconditioner
  std::vector<ASMbase*> patch; //!< The spline patches

  //! \brief Splits the local dofs into a number of unique subdomains (1D).
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains1D(const IntVec& nx,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    PetscIntMat& locSubds) const;

  //! \brief Splits the local dofs into a number of unique subdomains (2D).
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains2D(const IntVec& nx, const IntVec& ny,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    PetscIntMat& locSubds) const;

  //! \brief Splits the local dofs into a number of unique subdomains (3D).
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getLocalSubdomains3D(const IntVec& nx, const IntVec& ny, const IntVec& nz,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    PetscIntMat& locSubds) const;

  //! \brief Splits the dofs into a number of overlapping subdomains (1D).
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains1D(const IntVec& nx, int overlap,
		       PetscIntMat& locSubds) const;

  //! \brief Splits the local dofs into a number of unique subdomains (2D).
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains2D(const IntVec& nx, const IntVec& ny, int overlap,
		       PetscIntMat& locSubds) const;

  //! \brief Splits the local dofs into a number of unique subdomains (3D).
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[out] locSubds Global node number for each of the subdomains
  bool getSubdomains3D(const IntVec& nx, const IntVec& ny, const IntVec& nz,
		       int overlap, PetscIntMat& locSubds) const;

  //! \brief Splits the local dofs into a number of unique subdomains (1D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getLocalSubdomains1D(PetscIntMat& locSubds, const IntVec& nx,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    int f1, int f2) const;

  //! \brief Splits the local dofs into a number of unique subdomains (2D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getLocalSubdomains2D(PetscIntMat& locSubds,
			    const IntVec& nx, const IntVec& ny,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    int f1, int f2) const;

  //! \brief Splits the local dofs into a number of unique subdomains (3D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] minNodeId Minimum node number for each patch
  //! \param[in] maxNodeId Maximum node number for each patch
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getLocalSubdomains3D(PetscIntMat& locSubds,
			    const IntVec& nx, const IntVec& ny, const IntVec& nz,
			    const IntVec& minNodeId, const IntVec& maxNodeId,
			    int f1, int f2) const;

  //! \brief Splits the dofs into a number of overlapping subdomains (1D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getSubdomains1D(PetscIntMat& locSubds, const IntVec& nx,
		       int overlap, int f1, int f2) const;

  //! \brief Splits the local dofs into a number of unique subdomains (2D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx, ny Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getSubdomains2D(PetscIntMat& locSubds,
		       const IntVec& nx, const IntVec& ny,
		       int overlap, int f1, int f2) const;

  //! \brief Splits the local dofs into a number of unique subdomains (3D).
  //! \param[out] locSubds Global node number for each of the subdomains
  //! \param[in] nx, ny, nz Split each patch in npart smaller subdomains
  //! \param[in] overlap Overlap for subdomains
  //! \param[in] f1 First field to compute subdomains for
  //! \param[in] f2 Last field to compute subdomains for
  bool getSubdomains3D(PetscIntMat& locSubds,
		       const IntVec& nx, const IntVec& ny, const IntVec& nz,
		       int overlap, int f1, int f2) const;

  //! \brief Compute min and max global node number for each patch
  //! \param[out] minNodeId Min node number on all patches
  //! \param[out] maxNodeId Max node number on all patches
  bool getMinMaxNode(IntVec& minNodeId, IntVec& maxNodeId) const;
};


#endif
