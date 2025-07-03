// $Id$
//==============================================================================
//!
//! \file GlbL2projector.h
//!
//! \date Oct 16 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General integrand for L2-projection of secondary solutions.
//!
//==============================================================================

#ifndef _GLB_L2_PROJECTOR_H
#define _GLB_L2_PROJECTOR_H

#include "Integrand.h"
#include "LinAlgenums.h"
#include "MatVec.h"

class ASMbase;
class IntegrandBase;
class FunctionBase;
class LinSolParams;
class SparseMatrix;
class StdVector;
class ProcessAdm;

typedef std::vector<int>           IntVec;      //!< Vector of integers
typedef std::vector<size_t>        uIntVec;     //!< Vector of unsigned integers
typedef std::vector<FunctionBase*> FunctionVec; //!< Vector of functions


/*!
  \brief Abstract class for evaluating integrand or function.
*/

class L2Integrand
{
public:
  //! \brief The constructor initializes the patch reference.
  //! \param[in] patch ASM holding geometry to evaluate on
  L2Integrand(const ASMbase& patch) : m_patch(patch) {}

  //! \brief Evaluates the entity in a set of points.
  //! \param[out] sField Matrix with results
  //! \param[in] gpar Points to evaluate in
  virtual bool evaluate(Matrix& sField, const RealArray* gpar) const = 0;

  //! \brief Returns dimension of entity to evaluate.
  virtual size_t dim() const = 0;

protected:
  const ASMbase& m_patch; //!< Reference to ASM holding geometry
};


/*!
  \brief Evaluation class for secondary solutions of an integrand.
*/

class L2ProbIntegrand : public L2Integrand
{
public:
  //! \brief The constructor initializes the integrand reference.
  //! \param[in] patch ASM holding geometry to evaluate on
  //! \param[in] itg Integrand to evaluate
  L2ProbIntegrand(const ASMbase& patch, const IntegrandBase& itg);

  //! \brief Evaluates the secondary solutions in a set of points.
  //! \param[out] sField Matrix with results
  //! \param[in] gpar Points to evaluate in
  bool evaluate(Matrix& sField, const RealArray* gpar) const override;

  //! \brief Returns number of secondary solutions.
  size_t dim() const override;

private:
  const IntegrandBase& m_itg; //!< Reference to integrand
};


/*!
  \brief Evaluation class for functions.
*/

class L2FuncIntegrand : public L2Integrand
{
public:
  //! \brief The constructor initializes the function reference.
  //! \param[in] patch ASM holding geometry to evaluate on
  //! \param[in] func Function to evaluate
  L2FuncIntegrand(const ASMbase& patch, const FunctionBase& func);

  //! \brief Evaluates the function in a set of points.
  //! \param[out] sField Matrix with results
  //! \param[in] gpar Points to evaluate in
  bool evaluate(Matrix& sField, const RealArray* gpar) const override;

  //! \brief Returns number of function components.
  size_t dim() const override;

private:
  const FunctionBase& m_func; //!< Reference to function
};


/*!
  \brief General integrand for L2-projection of secondary solutions.
*/

class GlbL2 : public Integrand
{
public:
  static LinAlg::MatrixType MatrixType;   //!< Matrix type for projection
  static LinSolParams*      SolverParams; //!< Linear solver params projection

  //! \brief The constructor initializes the projection matrices.
  //! \param[in] p The main problem integrand
  //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
  GlbL2(IntegrandBase* p, size_t n);
  //! \brief Alternative constructor for projection of an explicit function.
  //! \param[in] f The function to do L2-projection on
  //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
  GlbL2(FunctionBase* f, size_t n);
  //! \brief Alternative constructor for projection of explicit functions.
  //! \param[in] f The functions to do L2-projection on
  //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
  GlbL2(const FunctionVec& f, size_t n);
  //! \brief The destructor frees the system matrix and system vector.
  virtual ~GlbL2();

  //! \brief Returns current solution mode.
  virtual SIM::SolutionMode getMode(bool) const { return SIM::RECOVERY; }
  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const;

  using Integrand::getLocalIntegral;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt Local integral for element
  virtual bool initElement(const IntVec& MNPC, const FiniteElement& fe,
                           const Vec3& X0, size_t nPt,
                           LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  virtual bool initElement(const IntVec& MNPC1,
                           const MxFiniteElement& fe,
                           const uIntVec& elem_sizes,
                           const uIntVec& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Dummy implementation.
  virtual bool initElement(const IntVec& MNPC1,
                           const uIntVec& elem_sizes,
                           const uIntVec& basis_sizes,
                           LocalIntegral& elmInt) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElement(const IntVec&, LocalIntegral&) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElementBou(const IntVec&, LocalIntegral&) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElementBou(const IntVec&, const uIntVec&, const uIntVec&,
                              LocalIntegral&) { return false; }

  using Integrand::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt,
                       const FiniteElement& fe, const Vec3& X) const;

  using Integrand::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const Vec3& X) const;

  //! \brief Pre-computes the sparsity pattern of the projection matrix \b A.
  //! \param[in] MMNPC Matrix of matrices of nodal point correspondances
  //! \param[in] nel Number of elements
  void preAssemble(const std::vector<IntVec>& MMNPC, size_t nel);

  //! \brief Solves the projection equation system and evaluates nodal values.
  //! \param[out] sField Nodal/control-point values of the projected results.
  bool solve(Matrix& sField);
  //! \brief Solves the projection equation system and evaluates nodal values.
  //! \param[out] sField Nodal/control-point values of the projected results.
  bool solve(const std::vector<Matrix*>& sField);

private:
  //! \brief Allocates the system L2-projection matrices.
  void allocate(size_t n);

public:
  mutable SparseMatrix* pA; //!< Left-hand-side matrix of the L2-projection
  mutable StdVector*    pB; //!< Right-hand-side vectors of the L2-projection

private:
  IntegrandBase* problem; //!< The main problem integrand
  FunctionVec  functions; //!< Explicit functions to L2-project
  size_t            nrhs; //!< Number of right-hand-size vectors
#ifdef HAS_PETSC
  ProcessAdm* adm; //!< Process administrator for PETSc
#endif
};

#endif
