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
#include "ASMs2D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"

#include "GoTools/geometry/SplineSurface.h"


SplineFields2D::SplineFields2D (const ASMs2D* patch,
                                const RealArray& v, const char* name)
  : Fields(2,name), basis(patch->getBasis()), surf(patch->getSurface())
{
  const int n1 = basis->numCoefs_u();
  const int n2 = basis->numCoefs_v();
  nno = n1*n2;

  const int p1 = basis->order_u();
  const int p2 = basis->order_v();
  nelm = (n1-p1+1)*(n2-p2+1);

  // Ensure the values array has compatible length, pad with zeros if necessary
  nf = v.size()/nno;
  values.resize(nf*nno);
  RealArray::const_iterator end = v.size() > nf*nno ? v.begin()+nf*nno:v.end();
  std::copy(v.begin(),end,values.begin());
}


bool SplineFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool SplineFields2D::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  Go::BasisPtsSf spline;
  basis->computeBasis(fe.u,fe.v,spline);

  // Evaluate the solution field at the given point
  std::vector<int> ip;
  ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		     basis->order_u(),basis->order_v(),
		     spline.left_idx,ip);

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}

/* Old procedure, way too complex...
   The above code is more in line with how we do it in ASMs2D.
  const int uorder = basis->order_u();
  const int vorder = basis->order_v();
  const int unum = basis->numCoefs_u();
  const int vnum = basis->numCoefs_v();

  const int dim = basis->dimension();
  const bool rational = basis->rational();
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

  basis->basis_u().computeBasisValues(fe.u,Bu.begin());
  basis->basis_v().computeBasisValues(fe.v,Bv.begin());
  
  const int uleft = basis->basis_u().lastKnotInterval();
  const int vleft = basis->basis_v().lastKnotInterval();

  //compute the tensor product value
  const int val_start_ix = (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * ncomp;
  
  register double* vtemp;
  register const double* val_ptr = &values[val_start_ix]; 
  std::fill(result.begin(), result.end(), double(0));

  if (rational) {
      double wtemp;
      double w = 0.0;

      const int w_start_ix  = (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * kdim + dim;
      register const double* w_ptr  = &(basis->rcoefs_begin()[w_start_ix]);
      
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
*/


bool SplineFields2D::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool SplineFields2D::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
/*
  Go::BasisDerivsSf spline;
  surf->computeBasis(fe.u,fe.v,spline);
  TODO: The above is not available yet, the below is temporary workaround.
*/
  std::vector<Go::BasisDerivsSf> tmpSpline(1);
  surf->computeBasisGrid(RealArray(1,fe.u),RealArray(1,fe.v),tmpSpline);
  Go::BasisDerivsSf& spline = tmpSpline.front();

  const int uorder = surf->order_u();
  const int vorder = surf->order_v();
  const size_t nen = uorder*vorder;

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  std::vector<int> ip;
  ASMs2D::scatterInd(surf->numCoefs_u(),surf->numCoefs_v(),
		     uorder,vorder,spline.left_idx,ip);

  // Evaluate the Jacobian inverse
  Matrix Xnod, Jac;
  Vector Xctrl(&(*surf->coefs_begin()),surf->coefs_end()-surf->coefs_begin());
  utl::gather(ip,surf->dimension(),Xctrl,Xnod);
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
    /*
    basis->computeBasis(fe.u,fe.v,spline);
    TODO: The above is not available yet, the below is temporary workaround.
    */
    basis->computeBasisGrid(RealArray(1,fe.u),RealArray(1,fe.v),tmpSpline);
    Go::BasisDerivsSf& spline = tmpSpline.front();

    const size_t nbf = basis->order_u()*basis->order_v();
    dNdu.resize(nbf,2);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    ASMs2D::scatterInd(basis->numCoefs_u(),basis->numCoefs_v(),
		       basis->order_u(),basis->order_v(),
		       spline.left_idx,ip);
  }

  utl::gather(ip,nf,values,Xnod);
  grad.multiply(Xnod,dNdX); // grad = Xnod * dNdX

  return true;
}

//   // Derivatives of solution in parametric space
//   std::vector<Go::Point> pts(3);
//   basis->pointValue(pts,values,fe.u,fe.v,1);

//   // Gradient of field wrt parametric coordinates
//   Matrix gradXi(nf,2);
//   for (size_t i = 1;i <= 2;i++) 
//     for (size_t n = 0;n < nf;n++) 
//       gradXi(n+1,i) = pts[i][n];
  
//   // Derivatives of coordinate mapping 
//   std::vector<Go::Point> xpts(3);
//   basis->point(xpts,fe.u,fe.v,1);

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
/*
  // Gradient of field 
  Matrix gradXi(nf,2);

  const int uorder = basis->order_u();
  const int vorder = basis->order_v();
  const int unum = basis->numCoefs_u();
  const int vnum = basis->numCoefs_v();

  const int derivs = 1;
  const bool u_from_right = true;
  const bool v_from_right = true;
  const double resolution = 1.0e-12;

  const int totpts = (derivs + 1)*(derivs + 2)/2;
  const int ncomp = values.size()/(unum*vnum);
  
  // Dimension
  const int dim       = basis->dimension();
  const bool rational = basis->rational();
  
  // Take care of the rational case
  const std::vector<double>::iterator co = rational ? basis->rcoefs_begin() : basis->coefs_begin();
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
    basis->basis_u().computeBasisValues(fe.u, &b0[0], derivs, resolution);
  } else {
    basis->basis_u().computeBasisValuesLeft(fe.u, &b0[0], derivs, resolution);
  }
  int uleft  = basis->basis_u().lastKnotInterval();

  if (v_from_right) {
    basis->basis_v().computeBasisValues(fe.v, &b1[0], derivs, resolution);
  } else {
    basis->basis_v().computeBasisValuesLeft(fe.v, &b1[0], derivs, resolution);
  }
  int vleft = basis->basis_v().lastKnotInterval();

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
*/


bool SplineFields2D::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
