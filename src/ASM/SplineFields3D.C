// $Id$
//==============================================================================
//!
//! \file SplineFields3D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "SplineFields3D.h"
#include "FiniteElement.h"
#include "Vec3.h"

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineUtils.h"
#include "boost/lambda/lambda.hpp"

namespace Go {
  void volume_ratder(double const eder[],int idim,int ider,double gder[]);
}


SplineFields3D::SplineFields3D(Go::SplineVolume *bf, char* name)
  : Fields(3,name), basis(bf), vol(bf) 
{
  const int n1 = basis->numCoefs(0);
  const int n2 = basis->numCoefs(1);
  const int n3 = basis->numCoefs(2);
  nno = n1*n2*n3;

  const int p1 = basis->order(0);
  const int p2 = basis->order(1);
  const int p3 = basis->order(2);
  nelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);

  // Number of fields set in fill
  nf = 0;
}


SplineFields3D::SplineFields3D(Go::SplineVolume *bf, 
			       Go::SplineVolume *geometry, 
			       char* name)
  : Fields(3,name), basis(bf), vol(geometry) 
{
  const int n1 = basis->numCoefs(0);
  const int n2 = basis->numCoefs(1);
  const int n3 = basis->numCoefs(2);
  nno = n1*n2*n3;

  const int p1 = basis->order(0);
  const int p2 = basis->order(1);
  const int p3 = basis->order(2);
  nelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);

  // Number of fields set in fill
  nf = 0;
}



SplineFields3D::~SplineFields3D()
{
  // Set geometry to NULL; should not be deallocated here
  basis = NULL;
} 


bool SplineFields3D::valueNode(int node, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields3D::valueFE(const FiniteElement& fe, Vector& vals) const
{
  if (!basis) return false;

  Go::Point pt;

//   basis->pointValue(pt,values,fe.u,fe.v,fe.w);
//   vals.resize(pt.size());
//   for (int i = 0;i < pt.size();i++)
//     vals[i] = pt[i];

  const int uorder = basis->order(0);
  const int vorder = basis->order(1);
  const int worder = basis->order(2);
  const int unum = basis->numCoefs(0);
  const int vnum = basis->numCoefs(1);
  const int wnum = basis->numCoefs(2);

  const int dim = basis->dimension();
  const bool rational = basis->rational();
  int kdim = rational ? dim + 1 : dim;
  
  const int ncomp = values.size()/(unum*vnum*wnum);
  
  static Go::ScratchVect<double, 10> Bu(uorder);
  static Go::ScratchVect<double, 10> Bv(vorder);
  static Go::ScratchVect<double, 10> Bw(worder);
  static Go::ScratchVect<double, 4> tempPt(ncomp);
  static Go::ScratchVect<double, 4> tempPt2(ncomp);
  static Go::ScratchVect<double, 4> tempResult(ncomp);
  
  Bu.resize(uorder);
  Bv.resize(vorder);
  Bw.resize(worder);
  tempPt.resize(ncomp);
  tempPt2.resize(ncomp);
  tempResult.resize(ncomp);
  
  // compute tbe basis values and get some data about the spline spaces
  basis->basis(0).computeBasisValues(fe.u, Bu.begin());
  basis->basis(1).computeBasisValues(fe.v, Bv.begin());
  basis->basis(2).computeBasisValues(fe.w, Bw.begin());
  const int uleft = basis->basis(0).lastKnotInterval();
  const int vleft = basis->basis(1).lastKnotInterval();
  const int wleft = basis->basis(2).lastKnotInterval();
  
  // compute the tensor product value
  const int val_start_ix =  (uleft - uorder + 1 + unum * (vleft - vorder + 1 + vnum * (wleft - worder + 1))) * ncomp;
  
  double* vtemp;
  const double* val_ptr = &values[val_start_ix];
  std::fill(pt.begin(), pt.end(), double(0));
  
  double w = 1.0;
  if (rational) {
    double tempw, tempw2;
    w = 0.0;

    const int w_start_ix = (uleft - uorder + 1 + unum * (vleft - vorder + 1 + vnum * (wleft - worder + 1))) * kdim + dim;
    const double *w_ptr = &(basis->rcoefs_begin()[w_start_ix]);
    
    for (double* bval_w_ptr = Bw.begin(); bval_w_ptr != Bw.end(); ++bval_w_ptr) {
      const double bval_w = *bval_w_ptr;
      std::fill(tempPt.begin(), tempPt.end(), 0);
      tempw = 0.0;
      for (double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	const double bval_v = *bval_v_ptr;
	std::fill(tempPt2.begin(), tempPt2.end(), 0);
	tempw2 = 0.0;
	for (double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	  const double bval_u = *bval_u_ptr;
	  for (vtemp = tempPt2.begin(); vtemp != tempPt2.end(); ++vtemp) {
	    *vtemp += bval_u * (*val_ptr++)*(*w_ptr);
	  }
	  tempw2 += bval_u*(*w_ptr);
	  w_ptr += kdim;
	}

	vtemp = tempPt2.begin();
	for (double* v = tempPt.begin(); v != tempPt.end(); ++v) {
	  *v += (*vtemp++) * bval_v;
	}

	tempw += tempw2*bval_v;
	val_ptr += ncomp * (unum - uorder);
	w_ptr   += kdim * (unum - uorder);
      }

      vtemp = tempPt.begin();
      for (double* v = tempResult.begin(); v != tempResult.end(); ++v) {
	*v += (*vtemp++) * bval_w;
      }
      w += tempw * bval_w;
      val_ptr += ncomp * unum * (vnum - vorder);
      w_ptr   += kdim  * unum * (vnum - vorder);
    }
  }
  else {
    for (double* bval_w_ptr = Bw.begin(); bval_w_ptr != Bw.end(); ++bval_w_ptr) {
      const double bval_w = *bval_w_ptr;
      std::fill(tempPt.begin(), tempPt.end(), 0);
      for (double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	const double bval_v = *bval_v_ptr;
	std::fill(tempPt2.begin(), tempPt2.end(), 0);
	for (double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	  const double bval_u = *bval_u_ptr;
	  for (vtemp = tempPt2.begin(); vtemp != tempPt2.end(); ++vtemp) {
	    *vtemp += bval_u * (*val_ptr++);
	  }
	}
	vtemp = tempPt2.begin();
	for (double* v = tempPt.begin(); v != tempPt.end(); ++v) {
	  *v += (*vtemp++) * bval_v;
	}
	val_ptr += ncomp * (unum - uorder);
      }
      vtemp = tempPt.begin();
      for (double* v = tempResult.begin(); v != tempResult.end(); ++v) {
	*v += (*vtemp++) * bval_w;
      }
      val_ptr += ncomp * unum * (vnum - vorder);
    }
  }
  
  std::copy(tempResult.begin(), tempResult.begin() + dim, pt.begin());
  if (rational) {
    const double w_inv = double(1) / w;
    std::transform(pt.begin(), pt.end(), pt.begin(), boost::lambda::_1 * w_inv);
  }    
  
  vals.resize(pt.size());
  for (int i = 0;i < pt.size();i++)
    vals[i] = pt[i];

  return true;
}



bool SplineFields3D::valueCoor(const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields3D::gradFE(const FiniteElement& fe, Matrix& grad) const
{
  if (!basis) return false;
  
//   // Derivatives of solution in parametric space
//   std::vector<Go::Point> pts;
//   basis->pointValue(pts,values,fe.u,fe.v,fe.w,1);

//  // Gradient of field wrt parametric coordinates 
//   Matrix gradXi(nf,3);
//   for (size_t i = 1;i <= 3;i++) 
//     for (size_t n = 0;n < nf;n++) 
//       gradXi(n+1,i) = pts[i][n];

//   // Derivatives of coordinate mapping 
//   std::vector<Go::Point> xpts(3);
//   basis->point(xpts,fe.u,fe.v,fe.w,1);

//   // Jacobian matrix
//   Matrix Jac(3,3);
//   for (size_t i = 1;i <= 3;i++)
//     for (size_t n = 0;n < 3;n++)
//       Jac(n+1,i) = xpts[i][n];
//   Jac.inverse();
  
//   // Gradient of solution in given parametric point
//   // df/dX = df/dXi * Jac^-1 = df/dXi * [dX/dXi]^-1
//   grad.multiply(gradXi,Jac);

  // Derivatives of solution in parametric space
  std::vector<Go::Point> pts(4);

  int derivs = 1;
  int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
  
  const int unum = basis->numCoefs(0);
  const int vnum = basis->numCoefs(1);
  const int wnum = basis->numCoefs(2);

  const int uorder = basis->order(0);
  const int vorder = basis->order(1);
  const int worder = basis->order(2);

  const bool u_from_right = true;
  const bool v_from_right = true;
  const bool w_from_right = true;

  const double resolution = 1.0e-12;
  
  const int ncomp     = values.size()/(unum*vnum*wnum);
  const int dim       = basis->dimension();
  const bool rational = basis->rational();
  
  for (int i = 0; i < totpts; ++i) {
    if (pts[i].dimension() != ncomp) {
      pts[i].resize(ncomp);
    }
  }
  
  // Take care of the rational case
  const std::vector<double>::iterator co = rational ? basis->rcoefs_begin() : basis->coefs_begin();
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
    basis->basis(0).computeBasisValues(fe.u, &b0[0], derivs, resolution);
  else 
    basis->basis(1).computeBasisValuesLeft(fe.u, &b0[0], derivs, resolution);
  
  int uleft = basis->basis(0).lastKnotInterval();
  if (v_from_right) 
    basis->basis(1).computeBasisValues(fe.v, &b1[0], derivs, resolution);
  else 
    basis->basis(1).computeBasisValuesLeft(fe.v, &b1[0], derivs, resolution);
  
  int vleft = basis->basis(1).lastKnotInterval();
  if (w_from_right) 
    basis->basis(2).computeBasisValues(fe.w, &b2[0], derivs, resolution);
  else 
    basis->basis(2).computeBasisValuesLeft(fe.w, &b2[0], derivs, resolution);
  
  int wleft = basis->basis(2).lastKnotInterval();
  
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
	  int temp_ind = ncomp;
	  weight = co[coefind*kdim+dim];

	  for (int wder = 0; wder <= derivs; ++wder) {
	    for (int vder = 0; vder <= wder; ++vder) {
	      for (int uder = 0; uder <= vder; ++uder) {
		temp2[temp_ind] += b0[id + wder - vder]*weight;
		temp_ind+=ndim;
	      }
	    }
	  }
	}

	const double *val_p = &values[coefind*ncomp]; 
	for (int d = 0; d < ncomp; ++d,++val_p) {
	  int temp_ind = d;
	  
	  for (int wder = 0; wder <= derivs; ++wder) {
	    for (int vder = 0; vder <= wder; ++vder) {
	      for (int uder = 0; uder <= vder; ++uder) {
		temp2[temp_ind] += b0[id + wder - vder]*(*val_p)*weight;
		temp_ind+=kdim;
	      }
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
    std::vector<double> restemp2(totpts*ncomp);
    Go::volume_ratder(&restemp[0], dim, derivs, &restemp2[0]);
    for (int i = 0; i < totpts; ++i) {
      for (int d = 0; d < ncomp; ++d) {
	pts[i][d] = restemp2[i*ndim + d];
      }
    }
  } else {
    double* restemp_it=restemp.begin();
    for (int i = 0; i < totpts; ++i) {
      for (int d = 0; d < ncomp; ++d) {
	pts[i][d] = *restemp_it;
	++restemp_it;
      }
    }
  }  

 // Gradient of field wrt parametric coordinates 
  Matrix gradXi(nf,3);
  for (size_t i = 1;i <= 3;i++) 
    for (size_t n = 0;n < nf;n++) 
      gradXi(n+1,i) = pts[i][n];

  // Derivatives of coordinate mapping 
  std::vector<Go::Point> xpts(3);
  vol->point(xpts,fe.u,fe.v,fe.w,1);

  // Jacobian matrix
  Matrix Jac(3,3);
  for (size_t i = 1;i <= 3;i++)
    for (size_t n = 0;n < 3;n++)
      Jac(n+1,i) = xpts[i][n];
  Jac.inverse();
  
  // Gradient of solution in given parametric point
  // df/dX = df/dXi * Jac^-1 = df/dXi * [dX/dXi]^-1
  grad.multiply(gradXi,Jac);

  return true;
}


bool SplineFields3D::gradCoor(const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
