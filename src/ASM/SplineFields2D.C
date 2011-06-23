// $Id$
//==============================================================================
//!
//! \file SplineFields2D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "SplineFields2D.h"
#include "FiniteElement.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineUtils.h"
#ifndef __BORLANDC__
#include "boost/lambda/lambda.hpp"
#endif


SplineFields2D::SplineFields2D(Go::SplineSurface *geometry, char* name)
  : Fields(2,name),  surf(geometry) 
{
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  nno = n1*n2;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  nelm = (n1-p1+1)*(n2-p2+1);

  // Number of fields set in fill
  nf = 0;
}


SplineFields2D::~SplineFields2D()
{
  // Set geometry to NULL; should not be deallocated here
  surf = NULL;
} 


bool SplineFields2D::valueNode(int node, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields2D::valueFE(const FiniteElement& fe, Vector& vals) const
{
  if (!surf) return false;

//   Go::Point pt;
//   surf->pointValue(pt,values,fe.u,fe.v);
//   vals.resize(pt.size());
//   for (int i = 0;i < pt.size();i++)
//     vals[i] = pt[i];

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const int unum = surf->numCoefs_u();
  const int vnum = surf->numCoefs_v();

  const int dim = surf->dimension();
  const bool rational = surf->rational();
  int kdim = rational ? dim + 1 : dim;

  const int ncomp = values.size()/(unum*vnum);
  Go::Point result(ncomp);

  static Go::ScratchVect<double, 10> Bu(uorder);
  static Go::ScratchVect<double, 10> Bv(vorder);
  static Go::ScratchVect<double, 4>  tempVal(ncomp);
  static Go::ScratchVect<double, 4>  tempResult(ncomp);

  Bu.resize(uorder);
  Bv.resize(vorder);
  tempVal.resize(ncomp);

  surf->basis_u().computeBasisValues(fe.u,Bu.begin());
  surf->basis_v().computeBasisValues(fe.v,Bv.begin());
  
  const int uleft = surf->basis_u().lastKnotInterval();
  const int vleft = surf->basis_v().lastKnotInterval();

  //compute the tensor product value
  const int val_start_ix = (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * ncomp;
  
  register double* vtemp;
  register const double* val_ptr = &values[val_start_ix]; 
  std::fill(result.begin(), result.end(), double(0));

  if (rational) {
      double wtemp;
      double w = 0.0;

      const int w_start_ix  = (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * kdim + dim;
      register const double* w_ptr  = &(surf->rcoefs_begin()[w_start_ix]);
      
      for (register double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
        register const double bval_v = *bval_v_ptr;
	std::fill(tempVal.begin(), tempVal.end(), 0);
	wtemp = 0.0;
        for (register double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
          register const double bval_u = *bval_u_ptr;
          for (vtemp = tempVal.begin(); vtemp != tempVal.end(); ++vtemp) {
            *vtemp += bval_u * (*val_ptr++)*(*w_ptr);
          }
          wtemp += bval_u * (*w_ptr);
          w_ptr += kdim;
        }

        vtemp = tempVal.begin();
        for (register double* v = result.begin(); v != result.end(); ++v) {
          *v += (*vtemp++) * bval_v;
        }
	w += wtemp * bval_v;

        val_ptr += ncomp * (unum - uorder);
        w_ptr   += kdim * (unum - uorder);
      }

      const double w_inv = double(1) / w;
#ifdef __BORLANDC__ //C++Builder does not support boost lambda
      std::transform(result.begin(), result.end(), result.begin(), ScaleBy(w_inv));
#else
      std::transform(result.begin(), result.end(), result.begin(), boost::lambda::_1 * w_inv);
#endif
    }
  else {
    for (register double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
      register const double bval_v = *bval_v_ptr;
      std::fill(tempVal.begin(), tempVal.end(), 0);
      for (register double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	register const double bval_u = *bval_u_ptr;
	for (vtemp = tempVal.begin(); vtemp != tempVal.end(); ++vtemp) {
	  *vtemp += bval_u * (*val_ptr++);
	}
      }
      vtemp = tempVal.begin();
      for (register double* v = result.begin(); v != result.end(); ++v) {
	*v += (*vtemp++) * bval_v;
      }
      val_ptr += ncomp * (unum - uorder);
    }
  } 

  vals.resize(result.size());
  for (int i = 0;i < result.size();i++)
    vals[i] = result[i];

  return true;
}



bool SplineFields2D::valueCoor(const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields2D::gradFE(const FiniteElement& fe, Matrix& grad) const
{
  if (!surf) return false;

//   // Derivatives of solution in parametric space
//   std::vector<Go::Point> pts(3);
//   surf->pointValue(pts,values,fe.u,fe.v,1);

//   // Gradient of field wrt parametric coordinates
//   Matrix gradXi(nf,2);
//   for (size_t i = 1;i <= 2;i++) 
//     for (size_t n = 0;n < nf;n++) 
//       gradXi(n+1,i) = pts[i][n];
  
//   // Derivatives of coordinate mapping 
//   std::vector<Go::Point> xpts(3);
//   surf->point(xpts,fe.u,fe.v,1);

//   // Jacobian matrix
//   Matrix Jac(2,2);
//   for (size_t i = 1;i <= 2;i++)
//     for (size_t n = 0;n < 2;n++)
//       Jac(n+1,i) = xpts[i][n];
  
//   // Invert Jacobian matrix
//   Jac.inverse();

//   // Gradient of solution in given parametric point
//   // df/dX = df/dXi * Jac^-1 = df/dXi * [dX/dXi]^-1
//   grad.multiply(gradXi,Jac);

  // Gradient of field 
  Matrix gradXi(nf,2);

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const int unum = surf->numCoefs_u();
  const int vnum = surf->numCoefs_v();

  const int derivs = 1;
  const bool u_from_right = true;
  const bool v_from_right = true;
  const double resolution = 1.0e-12;

  const int totpts = (derivs + 1)*(derivs + 2)/2;
  const int ncomp = values.size()/(unum*vnum);
  
  // Dimension
  const int dim       = surf->dimension();
  const bool rational = surf->rational();
  
  // Take care of the rational case
  const std::vector<double>::iterator co = rational ? surf->rcoefs_begin() : surf->coefs_begin();
  int kdim = dim   + (rational ? 1 : 0);
  int ndim = ncomp + (rational ? 1 : 0); 

  // Make temporary storage for the basis values and a temporary
  // computation cache.
  Go::ScratchVect<double, 30> b0(uorder*(derivs+1));
  Go::ScratchVect<double, 30> b1(vorder*(derivs+1));
  Go::ScratchVect<double, 30> temp(ndim * totpts);
  Go::ScratchVect<double, 30> restemp(ndim * totpts);
  std::fill(restemp.begin(), restemp.end(), 0.0);

  // Compute the basis values and get some data about the spline spaces
  if (u_from_right) {
    surf->basis_u().computeBasisValues(fe.u, &b0[0], derivs, resolution);
  } else {
    surf->basis_u().computeBasisValuesLeft(fe.u, &b0[0], derivs, resolution);
  }
  int uleft  = surf->basis_u().lastKnotInterval();

  if (v_from_right) {
    surf->basis_v().computeBasisValues(fe.v, &b1[0], derivs, resolution);
  } else {
    surf->basis_v().computeBasisValuesLeft(fe.v, &b1[0], derivs, resolution);
  }
  int vleft = surf->basis_v().lastKnotInterval();

  // Compute the tensor product value
  int coefind = uleft-uorder+1 + unum*(vleft-vorder+1);
  int derivs_plus1=derivs+1;
  double weight = 1.0;
  for (int jj = 0; jj < vorder; ++jj) {
    int jjd=jj*(derivs_plus1);
    std::fill(temp.begin(), temp.end(), 0.0);
    
    for (int ii = 0; ii < uorder; ++ii) {
      int iid=ii*(derivs_plus1);
      
      if (rational) {
	weight = co[coefind*kdim + dim];
	int temp_ind = ncomp;
        for (int vder = 0; vder < derivs_plus1; ++vder) {
          for (int uder = 0; uder < vder+1; ++uder) {
            temp[temp_ind] += b0[iid+vder - uder]*weight;
            temp_ind += ndim;
          }
        }
      }

      const double *val_p = &values[coefind*ncomp];
      for (int dd = 0; dd < ncomp; ++dd,++val_p) {
        int temp_ind=dd;
        for (int vder = 0; vder < derivs_plus1; ++vder) {
          for (int uder = 0; uder < vder+1; ++uder) {
            temp[temp_ind] += b0[iid+vder - uder]*(*val_p)*(weight);
            temp_ind += ndim;
          }
        }
      }

      coefind += 1;
    }
    
    for (int dd = 0; dd < ndim; ++dd) {
      int dercount = 0;
      for (int vder = 0; vder < derivs_plus1; ++vder) {
        for (int uder = 0; uder < vder + 1; ++uder) {
          restemp[dercount*ndim + dd] 
            += temp[dercount*ndim + dd]*b1[uder + jjd];
          ++dercount;
        }
      }
    }
    
    coefind += unum - uorder;
  }

  // Copy from restemp to result
  if (rational) {
    std::vector<double> restemp2(totpts*ncomp);
    Go::surface_ratder(&restemp[0], ncomp, derivs, &restemp2[0]);
    for (int i = 1; i < totpts; ++i) {
      for (int dd = 0; dd < ncomp; ++dd) 
	gradXi(dd+1,i) = restemp2[i*ncomp + dd];
    }
  } else {
    double* restemp_it=restemp.begin() + ncomp;
    for (int i = 1; i < totpts; ++i) {
      for (int dd = 0; dd < ncomp; ++dd) {
	gradXi(dd+1,i) = *restemp_it;
        ++restemp_it;
      }
    }
  }  

  // Derivatives of coordinate mapping 
  std::vector<Go::Point> xpts(3);
  surf->point(xpts,fe.u,fe.v,1);

  // Jacobian matrix
  Matrix Jac(2,2);
  for (size_t i = 1;i <= 2;i++)
    for (size_t n = 0;n < 2;n++)
      Jac(n+1,i) = xpts[i][n];
  
  // Invert Jacobian matrix
  Jac.inverse();

  // Gradient of solution in given parametric point
  // df/dX = df/dXi * Jac^-1 = df/dXi * [dX/dXi]^-1
  grad.multiply(gradXi,Jac);

  return true;
}


bool SplineFields2D::gradCoor(const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
