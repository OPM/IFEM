//==============================================================================
//!
//! \file SplineField3D.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for spline based finite element field in 3D
//!
//==============================================================================

#ifndef _SPLINE_FIELD_3D_H
#define _SPLINE_FIELD_3D_H

namespace Go {
  class SplineVolume;
  class Point;
  void volume_ratder(double const eder[],int idim,int ider,double gder[]);
}

#include "SplineField.h"
#include "GoTools/trivariate/SplineVolume.h"

/* 
   \brief Base class for spline-based finite element fields in 3D.

   \details This class implements the functions required to evaluate
   a 3D spline field at a given point in parametrical or physical 
   coordinates.
*/


class SplineField3D : public SplineField
{
 public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] geometry Spline geometry
  //! \param[in] name Name of spline field
  SplineField3D(Go::SplineVolume *geometry, char* name = NULL); 

  //! \brief Empty destructor
  virtual ~SplineField3D();

  // Methods to compute field values
  //================================================

  //! \brief Computes the value in a given node/control point
  //! \param[in] node Node number 
  double valueNode(int node) const;

  //! \brief Computes the value at a given local coordinate
  //! \param[in] fe Finite element definition
  double valueFE(const FiniteElement& fe) const;
  
  //! \brief Computed the value at a given global coordinate
  //! \param[in] x Global/physical coordinate for point
  double valueCoor(const Vec3& x) const;

  //! \brief Computes the gradient for a given local coordinate
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Vector& grad) const;

  //! \brief Computes the gradient for a given global/physical coordinate
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  bool gradCoor(const Vec3& x, Vector& grad) const;

 protected:
  Go::SplineVolume *vol;    //!< Spline volume geometry description

};
   


#endif
