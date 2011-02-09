//==============================================================================
//!
//! \file AnaSol.h
//!
//! \date Des 7 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Analytical solution fields (primary and secondary)
//!
//==============================================================================

#ifndef ANA_SOL_H
#define ANA_SOL_H

#include "Function.h"


/*!
  \brief Data type for analytical scalar solution fields (primary and secondary)
*/

class AnaSol
{
 protected:
  bool scalarSol;          //!< If scalar solution field defined
  bool vectorSol;          //!< If a vector solution field is defined
  

  // Scalar analytical solutions
  RealFunc* scalSol;       //<! Primary scalar solution field
  VecFunc*  scalSecSol;    //<! Secondary scalar solution fields

  // Vector solutions
  VecFunc*      vecSol;    //<! Primary Vec3 solution fields
  TensorFunc*   vecSecSol; //<! Secondary vector solution fields
  //GradientFunc* vecGrad;   //<! Gradient of solution fields

 public:
  //! \brief Constructor
  AnaSol();

  //! \brief Constructor
  //! \param[in] sp Primary scalar solution field
  //! \param[in] ss Secondary scalar solution field
  //! \param[in] vp Primary vector solution field
  //! \param[in] vs Secondary vector solution field
  AnaSol(RealFunc* sp, VecFunc* ss, VecFunc* vp, TensorFunc* vs);

  //! \brief Destructor
  virtual ~AnaSol();

  //! \brief Returns \e true if a scalar solution is defined
  virtual bool hasScalarSol() const { return scalarSol; }

  //! \brief Returns \e true if a vector solution is defined
  virtual bool hasVectorSol() const { return vectorSol;  }

  //! \brief Returns the scalar solution if any
  virtual RealFunc* getScalarSol()    const { return scalSol;    }

  //! \brief Returns the secondary scalar solution if any
  virtual VecFunc*  getScalarSecSol() const { return scalSecSol; }
  
  //! \brief Returns the vector solution if any
  virtual VecFunc*    getVectorSol()    const { return vecSol;    }

  //! \brief Returns the secondary vector solution if any
  virtual TensorFunc* getVectorSecSol() const { return vecSecSol; }
};

#endif
