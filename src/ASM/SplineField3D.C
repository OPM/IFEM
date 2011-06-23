// $Id$
//==============================================================================
//!
//! \file SplineField3D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element scalar field in 3D.
//!
//==============================================================================

#include "SplineField3D.h"
#include "FiniteElement.h"
#include "Vec3.h"

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineUtils.h"

namespace Go {
  void volume_ratder(double const eder[],int idim,int ider,double gder[]);
}


SplineField3D::SplineField3D(Go::SplineVolume *geometry, char* name)
  : Field(3,name), vol(geometry) 
{
  const int n1 = vol->numCoefs(0);
  const int n2 = vol->numCoefs(1);
  const int n3 = vol->numCoefs(2);
  nno = n1*n2*n3;

  const int p1 = vol->order(0);
  const int p2 = vol->order(1);
  const int p3 = vol->order(2);
  nelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);
}


SplineField3D::~SplineField3D()
{
  // Set geometry pointer to NULL, should not be deallocated here
  vol = NULL;
}


double SplineField3D::valueNode(int node) const
{
  // Not implemented yet
  return 0.0;
}


double SplineField3D::valueFE(const FiniteElement& fe) const
{
//   Go::Point pt;
//   vol->pointValue(pt,values,fe.u,fe.v,fe.w);
//   return pt[0];

  double val = 0.0;

  const int uorder = vol->order(0);
  const int vorder = vol->order(1);
  const int worder = vol->order(2);
  const int unum = vol->numCoefs(0);
  const int vnum = vol->numCoefs(1);

  const int dim = vol->dimension();
  const bool rational = vol->rational();
  int kdim = rational ? dim + 1 : dim;
  
  static Go::ScratchVect<double, 10> Bu(uorder);
  static Go::ScratchVect<double, 10> Bv(vorder);
  static Go::ScratchVect<double, 10> Bw(worder);
  double temp, temp2, tempResult = 0.0;
  
  Bu.resize(uorder);
  Bv.resize(vorder);
  Bw.resize(worder);
  
  // compute tbe basis values and get some data about the spline spaces
  vol->basis(0).computeBasisValues(fe.u, Bu.begin());
  vol->basis(1).computeBasisValues(fe.v, Bv.begin());
  vol->basis(2).computeBasisValues(fe.w, Bw.begin());
  const int uleft = vol->basis(0).lastKnotInterval();
  const int vleft = vol->basis(1).lastKnotInterval();
  const int wleft = vol->basis(2).lastKnotInterval();

  // compute the tensor product value
  const int val_start_ix = uleft - uorder + 1 + unum * (vleft - vorder + 1 + vnum * (wleft - worder + 1));
  const double* val_ptr = &values[val_start_ix];
  
  double w = 1.0;
  if (rational) {
    double tempw, tempw2;

    w = 0.0;
    const int w_start_ix = (uleft - uorder + 1 + unum * (vleft - vorder + 1 + vnum * (wleft - worder + 1))) * kdim + dim;
    const double *w_ptr = &(vol->rcoefs_begin()[w_start_ix]);
    
    for (double* bval_w_ptr = Bw.begin(); bval_w_ptr != Bw.end(); ++bval_w_ptr) {
      const double bval_w = *bval_w_ptr;
      temp = tempw = 0.0;
      for (double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	const double bval_v = *bval_v_ptr;
	temp2 = tempw2 = 0.0;
	for (double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	  const double bval_u = *bval_u_ptr;

	  temp2  += bval_u * (*val_ptr++) * (*w_ptr);
	  tempw2 += bval_u*(*w_ptr);
	  w_ptr += kdim;
	}

	temp  += temp2 * bval_v;
	tempw += tempw2 * bval_v;

	val_ptr += unum - uorder;
	w_ptr   += kdim * (unum - uorder);
      }

      tempResult += temp * bval_w;
      w += tempw * bval_w;

      val_ptr += unum * (vnum - vorder);
      w_ptr   += kdim  * unum * (vnum - vorder);
    }
  }
  else {
    for (double* bval_w_ptr = Bw.begin(); bval_w_ptr != Bw.end(); ++bval_w_ptr) {
      const double bval_w = *bval_w_ptr;
      temp = 0.0;
      for (double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	const double bval_v = *bval_v_ptr;
	temp2 = 0.0;
	for (double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	  const double bval_u = *bval_u_ptr;
	  temp2 += bval_u * (*val_ptr++);
	}

	temp += temp2 * bval_v;
	val_ptr += unum - uorder;
      }

      tempResult += temp * bval_w;
      val_ptr += unum * (vnum - vorder);
    }
  }
  
  val = tempResult/w;
  return val;
}


double SplineField3D::valueCoor(const Vec3& x) const
{
  // Not implemented yet
  return 0.0;
}


bool SplineField3D::gradFE(const FiniteElement& fe, Vector& grad) const
{
  if (!vol) return false;

//   // Derivatives of solution in parametric space
//   std::vector<Go::Point> pts;
//   vol->pointValue(pts,values,fe.u,fe.v,fe.w,1);

//   // Gradient of field wrt parametric coordinates
//   Vector gradXi(3);
//   gradXi(1) = pts[1][0];
//   gradXi(2) = pts[2][0];
//   gradXi(3) = pts[3][0];

//   // Derivatives of coordinate mapping 
//   std::vector<Go::Point> xpts(3);
//   vol->point(xpts,fe.u,fe.v,fe.w,1);

//   // Jacobian matrix
//   Matrix Jac(3,3);
//   for (size_t i = 1;i <= 3;i++)
//     for (size_t n = 0;n < 3;n++)
//       Jac(n+1,i) = xpts[i][n];

//   // Invert Jacobian matrix
//   Jac.inverse();

//   // Gradient of solution in given parametric point
//   // df/dX = Jac^-T * df/dXi = [dX/dXi]^-T * df/dXi
//   Jac.multiply(gradXi,grad,true,false);

  // Gradient of field wrt parametric coordinates
  Vector gradXi(3);

  int derivs = 1;
  int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
  
  const int unum = vol->numCoefs(0);
  const int vnum = vol->numCoefs(1);
  const int wnum = vol->numCoefs(2);

  const int uorder = vol->order(0);
  const int vorder = vol->order(1);
  const int worder = vol->order(2);

  const bool u_from_right = true;
  const bool v_from_right = true;
  const bool w_from_right = true;

  const double resolution = 1.0e-12;
  
  const int ncomp     = values.size()/(unum*vnum*wnum);
  const int dim       = vol->dimension();
  const bool rational = vol->rational();
  
  // Take care of the rational case
  const std::vector<double>::iterator co = rational ? vol->rcoefs_begin() : vol->coefs_begin();
  int kdim = dim + (rational ? 1 : 0);
  int ndim = ncomp + (rational ? 1 : 0);
  
  // Make temporary storage for the basis values and a temporary
  // computation cache.
  Go::ScratchVect<double, 30> b0(uorder * (derivs+1));
  Go::ScratchVect<double, 30> b1(vorder * (derivs+1));
  Go::ScratchVect<double, 30> b2(worder * (derivs+1));
  Go::ScratchVect<double, 60> temp(ndim * totpts);
  Go::ScratchVect<double, 60> temp2(ndim * totpts);
  Go::ScratchVect<double, 60> restemp(ndim * totpts);
  std::fill(restemp.begin(), restemp.end(), 0.0);
  
  // Compute the basis values and get some data about the spline spaces
  if (u_from_right) 
    vol->basis(0).computeBasisValues(fe.u, &b0[0], derivs, resolution);
  else 
    vol->basis(1).computeBasisValuesLeft(fe.u, &b0[0], derivs, resolution);
  
  int uleft = vol->basis(0).lastKnotInterval();
  if (v_from_right) 
    vol->basis(1).computeBasisValues(fe.v, &b1[0], derivs, resolution);
  else 
    vol->basis(1).computeBasisValuesLeft(fe.v, &b1[0], derivs, resolution);
  
  int vleft = vol->basis(1).lastKnotInterval();
  if (w_from_right) 
    vol->basis(2).computeBasisValues(fe.w, &b2[0], derivs, resolution);
  else 
    vol->basis(2).computeBasisValuesLeft(fe.w, &b2[0], derivs, resolution);
  
  int wleft = vol->basis(2).lastKnotInterval();
  
  // Compute the tensor product value
  int coefind = uleft-uorder+1 + unum*(vleft-vorder+1 + vnum*(wleft-worder+1));
  int derivs_plus1=derivs+1;
  
  double weight = 1.0;
  for (int k = 0; k < worder; ++k) {
    int kd=k*derivs_plus1;
    std::fill(temp.begin(), temp.end(), 0.0);
      
    for (int j = 0; j < vorder; ++j) {
      int jd=j*derivs_plus1;
      std::fill(temp2.begin(), temp2.end(), 0.0);
      
      for (int i = 0; i < uorder; ++i) {
	int id=i*derivs_plus1;

	if (rational) {
	  int temp_ind = 1;
	  weight  = co[coefind*kdim+dim];
	  for (int wder = 0; wder <= derivs; ++wder) {
	    for (int vder = 0; vder <= wder; ++vder) {
	      for (int uder = 0; uder <= vder; ++uder) {
		temp2[temp_ind] += b0[id + wder - vder]*weight;
		temp_ind+=ndim;
	      }
	    }
	  }
	}
	
	int temp_ind = 0;
	const double *val_p = &values[coefind];
	for (int wder = 0; wder <= derivs; ++wder) {
	  for (int vder = 0; vder <= wder; ++vder) {
	    for (int uder = 0; uder <= vder; ++uder) {
	      temp2[temp_ind] += b0[id + wder - vder]*(*val_p)*weight;
	      temp_ind+=kdim;
	    }
	  }
	}
	
	coefind += 1;
      }
      
      for (int d = 0; d < ndim; ++d) {
	int temp_ind = d;
	
	for (int wder = 0; wder <= derivs; ++wder) {
	  for (int vder = 0; vder <= wder; ++vder) {
	    for (int uder = 0; uder <= vder; ++uder) {
	      temp[temp_ind]
		+= temp2[temp_ind]*b1[jd + vder - uder];
	      temp_ind += ndim;
	      }
	  }
	}
      }
      coefind += unum - uorder;
    }
    
    for (int d = 0; d < ndim; ++d) {
      int temp_ind = d;
      
      for (int wder = 0; wder <= derivs; ++wder) {
	for (int vder = 0; vder <= wder; ++vder) {
	  for (int uder = 0; uder <= vder; ++uder) {
	    restemp[temp_ind]
	      += temp[temp_ind]*b2[kd + uder];
	    temp_ind += ndim;
	  }
	}
      }
    }
    coefind += unum * (vnum - vorder);
  }
  
  // Copy from restemp to result
  if (rational) {
    std::vector<double> restemp2(totpts);
    Go::volume_ratder(&restemp[0], dim, derivs, &restemp2[0]);
    for (int i = 1; i < totpts; ++i) 
      gradXi(i) = restemp2[i];
  } 
  else {
    double* restemp_it=restemp.begin();
    for (int i = 1; i < totpts; ++i) {
      gradXi(i) = *restemp_it;
      ++restemp_it;
    }
  }  
  
  // Derivatives of coordinate mapping 
  std::vector<Go::Point> xpts(3);
  vol->point(xpts,fe.u,fe.v,fe.w,1);

  // Jacobian matrix
  Matrix Jac(3,3);
  for (size_t i = 1;i <= 3;i++)
    for (size_t n = 0;n < 3;n++)
      Jac(n+1,i) = xpts[i][n];

  // Invert Jacobian matrix
  Jac.inverse();

  // Gradient of solution in given parametric point
  // df/dX = Jac^-T * df/dXi = [dX/dXi]^-T * df/dXi
  Jac.multiply(gradXi,grad,true,false);

  return true;
}


bool SplineField3D::gradCoor(const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
