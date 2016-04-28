// $Id$
//==============================================================================
//!
//! \file SAM.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices.
//!
//==============================================================================

#ifndef _SAM_H
#define _SAM_H

#include "MatVec.h"
#include <set>

class SystemMatrix;
class SystemVector;

typedef std::vector<int> IntVec; //!< General integer vector
typedef std::set<int>    IntSet; //!< General integer set


/*!
  \brief This class contains data and functions for the assembly of FE matrices.
  \details The names and meanings of (most of) the data members of this class
  are adopted from Kolbein Bell's pionering work on the field.
  See his reports on the SAM library for a thorough elaboration.

  The class does not contain methods for initializing the data members.
  That has to be done by deriving sub-classes specific to the solution methods.
*/

class SAM
{
public:
  //! \brief The constructor initializes an empty object.
  SAM();
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~SAM();

  //! \brief Prints out the key data to the given stream.
  void print(std::ostream& os) const;

  //! \brief Returns the number of elements in the model.
  int getNoElms() const { return nel; }
  //! \brief Returns the number of FE nodes in the model.
  //! \param[in] dofType Only consider nodes of this type (default All)
  int getNoNodes(char dofType = 'A') const;
  //! \brief Returns the total number of DOFs in the model.
  int getNoDOFs() const { return ndof; }
  //! \brief Returns the number of equations (free DOFs) in the model.
  virtual int getNoEquations() const { return neq; }
  //! \brief Returns the equations numbers for a given dofType.
  //! \param nodeType Type of node
  //! \param ldof DOFs in nodes (0 = all)
  IntSet getEquations(char nodeType, int dof=0) const;
  //! \brief Returns the Matrix of Accumulated DOFs.
  const int* getMADOF() const { return madof; }
  //! \brief Returns the Matrix of EQuation Numbers.
  const int* getMEQN() const { return meqn; }

  //! \brief Returns max number of DOF couplings in the model.
  int getMaxDofCouplings() const;

  //! \brief Computes number of couplings for each free DOF in the model.
  //! \param[out] nnz Number of couplings (non-zeroes) for each DOF
  virtual bool getNoDofCouplings(IntVec& nnz) const;
  //! \brief Interface to computation of DOF couplings for distributed matrices.
  virtual bool getNoDofCouplings(int, int, IntVec&, IntVec&) const
  { return false; }

  //! \brief Computes the sparse structure (DOF couplings) in the system matrix.
  //! \param[out] irow start index for each row in jcol
  //! \param[out] jcol column indices for non-zero entries
  bool getDofCouplings(IntVec& irow, IntVec& jcol) const;
  //! \brief Finds the set of free DOFs coupled to each free DOF.
  bool getDofCouplings(std::vector<IntSet>& dofc) const;

  //! \brief Initializes the system matrices prior to the element assembly.
  //! \param sysK   The system left-hand-side matrix to be initialized
  //! \param sysRHS The system right-hand-side load vector to be initialized
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \param[in] dontLockSP If \e true, do not lock the matrix sparsity pattern
  //! \return \e false if no free DOFs in the system, otherwise \e true
  //!
  //! \details This method must be called once before the first call to
  //! \a assembleSystem for a given load case or time step.
  bool initForAssembly(SystemMatrix& sysK, SystemVector& sysRHS,
		       Vector* reactionForces = nullptr,
		       bool dontLockSP = false) const;

  //! \brief Initializes a system matrix prior to the element assembly.
  //! \param sysM The system left-hand-side matrix to be initialized
  //! \return \e false if no free DOFs in the system, otherwise \e true
  //!
  //! \details This method must be called once before the first call to
  //! \a assembleSystem for a given load case or time step.
  bool initForAssembly(SystemMatrix& sysM) const;

  //! \brief Initializes the system load vector prior to the element assembly.
  //! \param sysRHS The system right-hand-side load vector to be initialized
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e false if no free DOFs in the system, otherwise \e true
  virtual bool initForAssembly(SystemVector& sysRHS,
			       Vector* reactionForces = nullptr) const;

  //! \brief Adds an element stiffness matrix into the system stiffness matrix.
  //! \param sysK    The left-hand-side system stiffness matrix
  //! \param sysRHS  The right-hand-side system load vector
  //! \param[in] eK  The element stiffness matrix
  //! \param[in] iel Identifier for the element that \a eK belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the right-hand-side system load vector.
  bool assembleSystem(SystemMatrix& sysK, SystemVector& sysRHS,
		      const Matrix& eK, int iel = 0,
		      Vector* reactionForces = nullptr) const;

  //! \brief Adds an element matrix into the corresponding system matrix.
  //! \param sysM    The left-hand-side system matrix
  //! \param[in] eM  The element matrix
  //! \param[in] iel Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  bool assembleSystem(SystemMatrix& sysM,
		      const Matrix& eM, int iel = 0) const;

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
			      Vector* reactionForces = nullptr) const;

  //! \brief Adds an element load vector into the system load vector.
  //! \param sysRHS  The right-hand-side system load vector
  //! \param[in] eS  The element load vector
  //! \param[in] iel Identifier for the element that \a eS belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assembleSystem(SystemVector& sysRHS,
			      const RealArray& eS, int iel = 0,
			      Vector* reactionForces = nullptr) const;

  //! \brief Adds a node load vector into the system load vector.
  //! \param sysRHS The right-hand-side system load vector
  //! \param[in] nS The nodal load vector
  //! \param[in] inod Identifier for the node that \a nS belongs to
  //! \param reactionForces Pointer to vector of nodal reaction forces
  //! \return \e true on successful assembly, otherwise \e false
  bool assembleSystem(SystemVector& sysRHS,
                      const Real* nS, int inod = 0,
                      Vector* reactionForces = nullptr) const;

  //! \brief Adds a node load vector into the system load vector.
  //! \param sysRHS The right-hand-side system load vector
  //! \param[in] S The nodal load vector
  //! \param[in] dof Node and local dof number for the load value \a S
  //! \return \e true on successful assembly, otherwise \e false
  bool assembleSystem(SystemVector& sysRHS, Real S,
                      const std::pair<int,int>& dof) const;

  //! \brief Adds a global load vector into the system load vector.
  //! \param sysRHS The right-hand-side system load vector
  //! \param[in] S  The global load vector
  virtual void addToRHS(SystemVector& sysRHS, const RealArray& S) const;

  //! \brief Finds the matrix of nodal point correspondance for an element.
  //! \param[out] mnpc Matrix of nodal point correspondance
  //! \param[in] iel Identifier for the element to get the node numbers for
  bool getElmNodes(IntVec& mnpc, int iel) const;
  //! \brief Finds the matrix of equation numbers for an element.
  //! \param[out] meen Matrix of element equation numbers
  //! \param[in] iel Identifier for the element to get the equation numbers for
  //! \param[in] nedof Number of degrees of freedom in the element
  //! (used for internal consistency checking, unless zero)
  virtual bool getElmEqns(IntVec& meen, int iel, int nedof = 0) const;
  //! \brief Returns the number equations for an element.
  //! \param[in] iel Identifier for the element to get number of equations for
  size_t getNoElmEqns(int iel) const;

  //! \brief Finds the matrix of equation numbers for a node.
  //! \param[out] mnen Matrix of node equation numbers
  //! \param[in] inod Identifier for the node to get the equation numbers for
  bool getNodeEqns(IntVec& mnen, int inod) const;

  //! \brief Returns the DOF classification of a given node.
  //! \param[in] inod Identifier for the node to get the classification for
  char getNodeType(int inod) const
  { return inod-- > 0 && inod < (int)nodeType.size() ? nodeType[inod] : ' '; }

  //! \brief Returns the first and last DOFs for a node.
  //! \param[in] inod Identifier for the node to get the DOF numbers for
  std::pair<int,int> getNodeDOFs(int inod) const;

  //! \brief Returns the internal node number and local index for a global DOF.
  //! \param[in] idof Global DOF-number in the range [1,NDOF]
  //! \return first = internal node number in the range [1,NNOD]
  //! \return second = local DOF number in the range [1,NNDOF]
  std::pair<int,int> getNodeAndLocalDof(int idof) const;

  //! \brief Finds the equation number corresponding to a local nodal DOF.
  //! \param[in] inod Identifier for the node to get the equation number for
  //! \param[in] ldof Local index of the DOF within node \a inod
  //! \return Equation number, or zero if the DOF is fixed or constrained
  int getEquation(int inod, int ldof) const;

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
  virtual bool expandSolution(const SystemVector& solVec, Vector& dofVec,
			      Real scaleSD = 1.0) const;

  //! \brief Expands a solution vector from equation-ordering to DOF-ordering.
  //! \param[in] solVec Solution vector, length = NEQ
  //! \param[out] dofVec Degrees of freedom vector, length = NDOF
  //! \return \e false if the length of \a solVec is invalid, otherwise \e true
  //!
  //! \details This version is typically used to expand eigenvectors.
  bool expandVector(const Vector& solVec, Vector& dofVec) const;

  //! \brief Applies the non-homogenous Dirichlet BCs to the given vector.
  //! \param dofVec Degrees of freedom vector, length = NDOF
  //!
  //! \details This method is typically used with explicit time integration.
  bool applyDirichlet(Vector& dofVec) const;

  //! \brief Computes the dot-product of two vectors of length NDOF.
  //! \param[in] x The first vector of the dot-product
  //! \param[in] y The second vector of the dot-product
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  virtual Real dot(const Vector& x, const Vector& y, char dofType = 'D') const;
  //! \brief Computes the l2-norm of a vector of length NDOF.
  //! \param[in] x The vector to compute the norm of
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  Real norm2(const Vector& x, char dofType = 'D') const
  { return sqrt(this->dot(x,x,dofType)); }
  //! \brief Computes the L2-norm of a vector of length NDOF.
  //! \param[in] x The vector to compute the norm of
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  virtual Real normL2(const Vector& x, char dofType = 'D') const;
  //! \brief Computes the L_infinity-norm of a vector of length NDOF.
  //! \param[in] x The vector to compute the norm of
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  virtual Real normInf(const Vector& x, size_t& comp, char dofType = 'D') const;

  //! \brief Computes the energy norm contributions from nodal reaction forces.
  //! \param[in] u The (incremental) nodal displacement vector
  //! \param[in] rf Compressed reaction force vector for the entire model
  virtual Real normReact(const Vector& u, const Vector& rf) const;
  //! \brief Returns the total reaction force in the given coordinate direction.
  //! \param[in] dir 1-based coordinate direction index
  //! \param[in] rf Compressed reaction force vector for the entire model
  virtual Real getReaction(int dir, const Vector& rf) const;
  //! \brief Returns a vector of reaction forces for a given node.
  //! \param[in] inod Identifier for the node to get the reaction forces for
  //! \param[in] rf Compressed reaction force vector for the entire model
  //! \param[out] nrf Nodal reaction forces
  //! \return \e true if the specified node has reaction forces
  bool getNodalReactions(int inod, const Vector& rf, Vector& nrf) const;

protected:
  //! \brief Initializes the DOF-to-equation connectivity array \a MEQN.
  virtual bool initSystemEquations();

  //! \brief Adds a scalar value into a system right hand-side vector.
  //! \param RHS The right-hand-side system load vector
  //! \param[in] value The scalar value to add in
  //! \param[in] ieq Global equation number for the scalar value
  //!
  //! \details If \a ieq is zero, this is a fixed DOF and the value is ignored.
  //! If \a ieq is less than zero, this is a constrained DOF and the value is
  //! added to the master DOFs of the governing constraint equation of this DOF.
  void assembleRHS(Real* RHS, Real value, int ieq) const;

  //! \brief Assembles reaction forces for the fixed and prescribed DOFs.
  //! \param reac The vector of reaction forces
  //! \param[in] eS  The element load vector
  //! \param[in] iel Identifier for the element that \a eS belongs to
  void assembleReactions(Vector& reac, const RealArray& eS, int iel) const;

  //! \brief Expands a solution vector from equation-ordering to DOF-ordering.
  //! \param[in] solVec Pointer to solution vector, length = NEQ
  //! \param[out] dofVec Degrees of freedom vector, length = NDOF
  //! \param[in] scaleSD Scaling factor for specified (slave) DOFs
  bool expandVector(const Real* solVec, Vector& dofVec, Real scaleSD) const;

private:
  int mpar[50]; //!< Matrix of parameters

protected:
  // The following parameters are pointers to specific locations in MPAR
  int& nnod;   //!< Number of nodes
  int& nel;    //!< Number of elements
  int& ndof;   //!< Number of DOFs
  int& nspdof; //!< Number of specified DOFs
  int& nceq;   //!< Number of constraint equations
  int& neq;    //!< Number of system equations
  int& nmmnpc; //!< Number of elements in MMNPC
  int& nmmceq; //!< Number of elements in MMCEQ

  // The standard SAM arrays (see K. Bell's reports for detailed explanation).
  // We are using plane C-pointers for these items such that they more easily
  // can be passed directly as arguments to FORTRAN subroutines.
  int*  mpmnpc; //!< Matrix of pointers to MNPCs in MMNPC
  int*  mmnpc;  //!< Matrix of matrices of nodal point correspondances
  int*  madof;  //!< Matrix of accumulated DOFs
  int*  msc;    //!< Matrix of status codes
  int*  mpmceq; //!< Matrix of pointers to MCEQs in MMCEQ
  int*  mmceq;  //!< Matrix of matrices of constraint equation definitions
  Real* ttcc;   //!< Table of tables of constraint equation coefficients
  int*  minex;  //!< Matrix of internal to external node numbers
  int*  meqn;   //!< Matrix of equation numbers

  std::vector<char> nodeType; //!< Nodal DOF classification

  friend class DenseMatrix;
  friend class SPRMatrix;
  friend class SparseMatrix;
  friend class PETScMatrix;
  friend class PETScBlockMatrix;
};

#endif
