// $Id$
//==============================================================================
//!
//! \file SplineField2D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element scalar field in 2D.
//!
//==============================================================================

#include "SplineField2D.h"
#include "FiniteElement.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineUtils.h"


SplineField2D::SplineField2D(Go::SplineSurface *geometry, char* name)
  : Field(2,name), surf(geometry) 
{
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  nno = n1*n2;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  nelm = (n1-p1+1)*(n2-p2+1);
}


SplineField2D::~SplineField2D()
{
  // Set geometry pointer to NULL, should not be deallocated here
  surf = NULL;
}


double SplineField2D::valueNode(int node) const
{
  // Not implemented yet
  return 0.0;
}


double SplineField2D::valueFE(const FiniteElement& fe) const
{
//   Go::Point pt;
//   surf->pointValue(pt,values,fe.u,fe.v);
//   return pt[0];

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const int unum = surf->numCoefs_u();
  //const int vnum = surf->numCoefs_v();

  const int dim = surf->dimension();
  const bool rational = surf->rational();
  const int kdim = rational ? dim + 1 : dim;

  //const int ncomp = values.size()/(unum*vnum);
  double val = 0.0;
  double tempVal;
 
  static Go::ScratchVect<double, 10> Bu(uorder);
  static Go::ScratchVect<double, 10> Bv(vorder);

  Bu.resize(uorder);
  Bv.resize(vorder);

  surf->basis_u().computeBasisValues(fe.u,Bu.begin());
  surf->basis_v().computeBasisValues(fe.v,Bv.begin());
  
  const int uleft = surf->basis_u().lastKnotInterval();
  const int vleft = surf->basis_v().lastKnotInterval();

  // compute the tensor product value
  const int val_start_ix = (uleft - uorder + 1 + unum * (vleft - vorder + 1));

  register const double* val_ptr = &values[val_start_ix]; 

  if (rational) {
      double w = 0.0;
      const int w_start_ix  = (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * kdim + dim;
      register const double* w_ptr  = &(surf->rcoefs_begin()[w_start_ix]);
      
      for (register double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
        register const double bval_v = *bval_v_ptr;
	tempVal = 0.0;
        for (register double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
          register const double bval_u = *bval_u_ptr;
	  tempVal += bval_u * (*val_ptr++);

          w += bval_u * (*w_ptr);
          w_ptr += kdim;
        }

	val += tempVal * bval_v;
	val_ptr += unum - uorder;

        w_ptr   += kdim * (unum - uorder);
      }

      val /= w;
  }
  else {
    for (register double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
      register const double bval_v = *bval_v_ptr;
      tempVal = 0.0;
      for (register double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	register const double bval_u = *bval_u_ptr;
	tempVal += bval_u * (*val_ptr++);
      }

      val += tempVal * bval_v;
      val_ptr += unum - uorder;
    }
  } 
  
  return val;
}


double SplineField2D::valueCoor(const Vec3& x) const
{
  // Not implemented yet
  return 0.0;
}


bool SplineField2D::gradFE(const FiniteElement& fe, Vector& grad) const
{
  if (!surf) return false;

//   // Derivatives of solution in parametric space
//   std::vector<Go::Point> pts(3);
//   surf->pointValue(pts,values,fe.u,fe.v,1);

//   // Gradient of field wrt parametric coordinates
//   Vector gradXi(2);
//   gradXi(1) = pts[1][0];
//   gradXi(2) = pts[2][0];

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
//   // df/dX = Jac^-T * df/dXi = [dX/dXi]^-T * df/dXi
//   Jac.multiply(gradXi,grad,true,false);
  
  // Gradient of field wrt parametric coordinates
  Vector gradXi(2);

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const int unum = surf->numCoefs_u();

  const int derivs = 1;
  const int totpts = 3;
  const bool u_from_right = true;
  const bool v_from_right = true;
  const double resolution = 1.0e-12;

  // Dimension
  const int dim       = surf->dimension();
  const bool rational = surf->rational();
  
  // Take care of the rational case
  const std::vector<double>::iterator co = rational ? surf->rcoefs_begin() : surf->coefs_begin();
  int kdim = dim   + (rational ? 1 : 0);
  int ndim = 1 + (rational ? 1 : 0); 
  
  // Make temporary storage for the basis values and a temporary
  // computation cache.
  Go::ScratchVect<double, 30> b0(uorder*(derivs+1));
  Go::ScratchVect<double, 30> b1(vorder*(derivs+1));
  Go::ScratchVect<double, 30> temp(ndim * totpts);
  Go::ScratchVect<double, 30> restemp(ndim * totpts);
  std::fill(restemp.begin(), restemp.end(), 0.0);

  // Compute the basis values and get some data about the spline spaces
  if (u_from_right) 
    surf->basis_u().computeBasisValues(fe.u, &b0[0], derivs, resolution);
  else 
    surf->basis_u().computeBasisValuesLeft(fe.u, &b0[0], derivs, resolution);  
  int uleft  = surf->basis_u().lastKnotInterval();

  if (v_from_right) 
    surf->basis_v().computeBasisValues(fe.v, &b1[0], derivs, resolution);
  else 
    surf->basis_v().computeBasisValuesLeft(fe.v, &b1[0], derivs, resolution);  
  int vleft = surf->basis_v().lastKnotInterval();

  // Compute the tensor product value
  int coefind = uleft-uorder+1 + unum*(vleft-vorder+1);
  int derivs_plus1=derivs+1;
  for (int jj = 0; jj < vorder; ++jj) {
    int jjd=jj*(derivs_plus1);
    std::fill(temp.begin(), temp.end(), 0.0);
    
    for (int ii = 0; ii < uorder; ++ii) {
      int iid=ii*(derivs_plus1);
      const double *val_p = &values[coefind];
      const double *co_p  = &co[coefind*kdim] + dim;
      int temp_ind = 0;
      for (int vder = 0; vder < derivs_plus1; ++vder) 
	for (int uder = 0; uder < vder+1; ++uder) {
	  temp[temp_ind] += b0[iid+vder - uder]*(*val_p);
	  temp_ind += ndim;
	}
      
      for (int dd = dim;dd < kdim;dd++,co_p += kdim) {
	int temp_ind = 1;
        for (int vder = 0; vder < derivs_plus1; ++vder) {
          for (int uder = 0; uder < vder+1; ++uder) {
            temp[temp_ind] += b0[iid+vder - uder]*(*co_p);
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
    std::vector<double> restemp2(totpts);
    Go::surface_ratder(&restemp[0],1, derivs, &restemp2[0]);
    for (int i = 1; i < totpts; ++i) 
      gradXi(i) = restemp2[i];
  } else {
    double* restemp_it=restemp.begin();
    ++restemp_it;
    for (int i = 1; i < totpts; ++i) {
      gradXi(i) = *restemp_it;
      ++restemp_it;
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
  // df/dX = Jac^-T * df/dXi = [dX/dXi]^-T * df/dXi
  Jac.multiply(gradXi,grad,true,false);
  
  return true;
}


bool SplineField2D::gradCoor(const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
