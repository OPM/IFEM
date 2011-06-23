// $Id$
//==============================================================================
//!
//! \file IntegrandBase.h
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base classes representing FEM integrands.
//!
//==============================================================================

#ifndef _INTEGRAND_BASE_H
#define _INTEGRAND_BASE_H

#include "Integrand.h"
#include "SIMenums.h"
#include "Function.h"

class NormBase;
class ElmMats;
class ElmNorm;
class AnaSol;
class VTF;


/*!
  \brief Base class representing a system level integrated quantity.
*/

class IntegrandBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  IntegrandBase() : npv(0) {}

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~IntegrandBase();

  //! \brief Prints out the problem definition to the given output stream.
  virtual void print(std::ostream&) const {}


  // Global initialization interface
  // ===============================

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode) {}
  //! \brief Initializes the integrand for a new integration loop.
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain.
  virtual void initIntegration(const TimeDomain&) {}
  //! \brief Initializes the integrand for a new result point loop.
  //! \details This method is invoked once before starting the evaluation of
  //! the secondary solution at all result sampling points, after the converged
  //! primary solution has been found.
  virtual void initResultPoints(double) {}


  // Element-level initialization interface
  // ======================================

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  virtual bool initElement(const std::vector<int>& MNPC);
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1);


  // Solution field evaluation interface
  // ===================================

  //! \brief Evaluates the analytical primary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (vector field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalPrimSol(Vector& s,
			   const VecFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const TensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (symmetric tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical primary solution at a result point.
  //! \param[out] s The solution field value at current point
  //! \param[in] asol The analytical solution field (scalar field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalPrimSol(double& s,
			   const RealFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (vector field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const VecFunc& asol, const Vec3& X) const;


  // Various service methods
  // =======================

  //! \brief Returns whether a mixed formulation is used.
  virtual bool mixedFormulation() const { return false; }

  //! \brief Writes surface tractions/fluxes for a given time step to VTF-file.
  virtual bool writeGlvT(VTF*, int, int&) const { return true; }

  //! \brief Returns whether there are any traction/flux values to write to VTF.
  virtual bool hasTractionValues() const { return false; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const { return 0; }

  //! \brief Returns the number of primary/secondary solution field components.
  virtual size_t getNoFields(int = 2) const { return 0; }
  //! \brief Returns the name of a primary solution field component.
  virtual const char* getField1Name(size_t, const char* = 0) const { return 0; }
  //! \brief Returns the name of a secondary solution field component.
  virtual const char* getField2Name(size_t, const char* = 0) const { return 0; }

  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const { return primsol.size(); }

  //! \brief Accesses the primary solution vector of current patch.
  Vector& getSolution(size_t n = 0) { return primsol[n]; }

protected:
  Vectors  primsol; //!< Primary solution vectors for current patch
  ElmMats* myMats;  //!< Local element matrices
  Vectors  mySols;  //!< Local element solution vectors

  unsigned short int npv; //!< Number of primary solution variables per node
};


/*!
  \brief Base class representing a system level norm quantity.
*/

class NormBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  NormBase(IntegrandBase& p) : myProblem(p) {}

public:
  //! \brief Empty destructor.
  virtual ~NormBase() {}

  //! \brief Initializes the integrand for a new integration loop.
  virtual void initIntegration(const TimeDomain& time);

  //! \brief Initializes current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC);

  //! \brief Initializes current element for numerical integration (mixed).
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1);

  //! \brief Initializes current element for boundary integration.
  virtual bool initElementBou(const std::vector<int>& MNPC);

  //! \brief Initializes current element for boundary integration (mixed).
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1);

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Returns a 1-based index of the external energy norm.
  virtual size_t indExt() const { return 0; }

  //! \brief Returns the number of field components.
  virtual size_t getNoFields() const { return 0; }

protected:
  //! \brief Returns the element norm object to use in the integration.
  //! \param elmInt The local integral object to receive norm contributions
  //! \param[in] nn Size of static norm buffer
  //!
  //! \details If \a elmInt is NULL or cannot be casted to an ElmNorm pointer,
  //! a local static buffer is used instead.
  static ElmNorm& getElmNormBuffer(LocalIntegral*& elmInt, const size_t nn = 4);

protected:
  IntegrandBase& myProblem; //!< The problem-specific data
};

#endif
