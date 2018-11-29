// $Id$
//==============================================================================
//!
//! \file ASMu2Dnurbs.h
//!
//! \date Nov 28 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D NURBS FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_NURBS_H
#define _ASM_U2D_NURBS_H

#include "ASMu2D.h"


/*!
  \brief Driver for assembly of unstructured 2D NURBS FE models.
  \details This class contains methods common for 2D LR-NURBS patches.
*/

class ASMu2Dnurbs : public ASMu2D
{
public:
  //! \brief Default constructor.
  ASMu2Dnurbs(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Copy constructor.
  ASMu2Dnurbs(const ASMu2Dnurbs& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu2Dnurbs() {}

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);

protected:
  //! \brief Evaluates the basis functions and derivatives of an element.
  virtual bool evaluateBasis(int iel, FiniteElement& fe, int derivs) const;

  //! \brief Evaluate basis functions in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisPtsSf& bas, int iel,
                            const LR::LRSplineSurface* spline) const;

  //! \brief Evaluate basis functions and first derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf& bas, int iel,
                            const LR::LRSplineSurface* spline) const;
  //! \brief Evaluate basis functions and two derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf2& bas, int iel) const;

  //! \brief Evaluate basis functions and two derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf3& bas, int iel) const;

  //! \brief Converts current tensor spline object to LR-spline.
  virtual LR::LRSplineSurface* createLRfromTensor();

private:
  bool noNurbs; //!< If \e true, we read a spline and thus forward to ASMu2D
};

#endif
