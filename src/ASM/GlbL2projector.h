// $Id$
//==============================================================================
//!
//! \file GlbL2projector.h
//!
//! \date Oct 16 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General integrands for global L2-projection.
//!
//==============================================================================

#ifndef _GLB_L2_PROJECTOR_H
#define _GLB_L2_PROJECTOR_H

#include "LinAlgenums.h"
#include "MatVec.h"

class ASMbase;
class IntegrandBase;
class FunctionBase;
class LinSolParams;
class ProcessAdm;


namespace GlbL2 //! Global control parameters for L2-projection
{
  extern LinAlg::MatrixType MatrixType;
  extern LinSolParams*      SolverParams;
  extern const ProcessAdm*  Adm;
}


/*!
  \brief Abstract class for evaluating integrand or function.
*/

class L2Integrand
{
public:
  //! \brief The constructor initializes the patch reference.
  //! \param[in] patch ASM holding geometry to evaluate on
  L2Integrand(const ASMbase& patch) : m_patch(patch) { myTime = 0.0; }
  //! \brief Empty destructor.
  virtual ~L2Integrand() {}

  //! \brief Evaluates the entity in a set of points.
  //! \param[out] sField Matrix with results
  //! \param[in] gpar Points to evaluate in
  virtual bool evaluate(Matrix& sField, const RealArray* gpar) const = 0;

  //! \brief Returns dimension of entity to evaluate.
  virtual size_t dim() const = 0;

  //! \brief Returns current time/load parameter for 2ndary solution evaluation.
  double getTimeLevel() const { return myTime; }

protected:
  const ASMbase& m_patch; //!< Reference to ASM holding geometry
  double         myTime;  //!< Current time or load parameter
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
  //! \param[in] time Current time
  L2FuncIntegrand(const ASMbase& patch,
                  const FunctionBase& func,
                  double time = 0.0);

  //! \brief Evaluates the function in a set of points.
  //! \param[out] sField Matrix with results
  //! \param[in] gpar Points to evaluate in
  bool evaluate(Matrix& sField, const RealArray* gpar) const override;

  //! \brief Returns number of function components.
  size_t dim() const override;

private:
  const FunctionBase& m_func; //!< Reference to function
};

#endif
