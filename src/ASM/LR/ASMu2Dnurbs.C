// $Id$
//==============================================================================
//!
//! \file ASMu2Dnurbs.C
//!
//! \date Nov 28 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D NURBS FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu2Dnurbs.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Profiler.h"


ASMu2Dnurbs::ASMu2Dnurbs (unsigned char n_s, unsigned char n_f)
  : ASMu2D(n_s, n_f)
{
  noNurbs = false;
}


ASMu2Dnurbs::ASMu2Dnurbs (const ASMu2Dnurbs& patch, unsigned char n_f)
  : ASMu2D(patch, n_f)
{
  noNurbs = patch.noNurbs;
}


bool ASMu2Dnurbs::read (std::istream& is)
{
  bool ok = this->ASMu2D::read(is);
  if (ok && !(tensorspline && tensorspline->rational()))
    std::cout <<"  ** LR-nurbs requested but input is a spline."<< std::endl;

  return ok;
}


bool ASMu2Dnurbs::evaluateBasis (int iel, FiniteElement& fe,
                                 int derivs) const
{
  if (noNurbs)
    return this->ASMu2D::evaluateBasis(iel,fe,derivs);

  const LR::Element* el = lrspline->getElement(iel);
  if (!el) return false;

  PROFILE3("ASMu2Dn::evalBasis");

  fe.xi  = 2.0*(fe.u - el->umin()) / (el->umax() - el->umin()) - 1.0;
  fe.eta = 2.0*(fe.v - el->vmin()) / (el->vmax() - el->vmin()) - 1.0;
  RealArray Nu, Nv;
#pragma omp critical
  {
    Nu = bezier_u.computeBasisValues(fe.xi, derivs);
    Nv = bezier_v.computeBasisValues(fe.eta,derivs);
  }
  const Matrix& C = bezierExtract[iel];

  RealArray w; w.reserve(el->nBasisFunctions());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  if (derivs < 1) {
    Matrix B;
    B.outer_product(Nu,Nv);
    fe.N = C*static_cast<const Vector&>(B);
    double W = fe.N.dot(w);
    for (size_t i = 0; i < fe.N.size(); i++)
      fe.N[i] *= w[i]/W;

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(static_cast<const Vector&>(B).sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  else {
    int p = lrspline->order(0)*lrspline->order(1);

    Vector B(p);
    Vector Bu(p); // Bezier basis functions differentiated wrt u
    Vector Bv(p); // Bezier basis functions differentiated wrt v

    size_t i, j, k = 0;
    for (j = 0; j < Nv.size(); j += derivs+1)
      for (i = 0; i < Nu.size(); i += derivs+1, k++) {
        B[k]  = Nu[i  ]*Nv[j  ];
        Bu[k] = Nu[i+1]*Nv[j  ];
        Bv[k] = Nu[i  ]*Nv[j+1];
      }

    fe.N = C*B;
    Vector dNxi  = C*Bu;
    Vector dNeta = C*Bv;

    double W       = fe.N.dot(w);
    double Wderxi  = dNxi.dot(w);
    double Wdereta = dNeta.dot(w);

    Matrix dNdu(w.size(),2);
    for (i = 1; i <= fe.N.size(); i++) {
      fe.N(i)  *= w[i-1]/W;
      dNdu(i,1) = (dNxi(i)*W  - fe.N(i)*Wderxi )*w[i-1]/(W*W);
      dNdu(i,2) = (dNeta(i)*W - fe.N(i)*Wdereta)*w[i-1]/(W*W);
    }

    Matrix Xnod, Jac;
    this->getElementCoordinates(Xnod,iel+1);
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(B.sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(1).sum()) > 1.0e-10)
      std::cerr <<"dNdu not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(2).sum()) > 1.0e-10)
      std::cerr <<"dNdv not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(Bu.sum()) > 1.0e-10 || fabs(Bv.sum()) > 1.0e-10)
      std::cerr <<"Bezier derivatives do not sum to zero at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  return true;
}


void ASMu2Dnurbs::computeBasis (double u, double v,
                                Go::BasisPtsSf& bas, int iel,
                                const LR::LRSplineSurface* spline) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel,spline);

  PROFILE3("ASMu2Dn::compBasis(0)");

  if (!spline)
    spline = lrspline.get();

  const LR::Element* el = spline->getElement(iel);

  Go::BasisPtsSf tmp;
  spline->computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W = w.dot(tmp.basisValues);

  bas.preparePts(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;
}


void ASMu2Dnurbs::computeBasis (double u, double v,
                                Go::BasisDerivsSf& bas, int iel,
                                const LR::LRSplineSurface* spline) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel,spline);

  PROFILE3("ASMu2Dn::compBasis(1)");

  if (!spline)
    spline = lrspline.get();

  const LR::Element* el = spline->getElement(iel);

  Go::BasisDerivsSf tmp;
  spline->computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W  = w.dot(tmp.basisValues);
  double Wx = w.dot(tmp.basisDerivs_u);
  double Wy = w.dot(tmp.basisDerivs_v);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
  {
    bas.basisValues[i]   =  tmp.basisValues[i]*w[i]/W;
    bas.basisDerivs_u[i] = (tmp.basisDerivs_u[i]*w[i] - bas.basisValues[i]*Wx)/W;
    bas.basisDerivs_v[i] = (tmp.basisDerivs_v[i]*w[i] - bas.basisValues[i]*Wy)/W;
  }
}


void ASMu2Dnurbs::computeBasis (double u, double v,
                                Go::BasisDerivsSf2& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  PROFILE3("ASMu2Dn::compBasis(2)");

  const LR::Element* el = lrspline->getElement(iel);

  Go::BasisDerivsSf2 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W   = w.dot(tmp.basisValues);
  double Wx  = w.dot(tmp.basisDerivs_u);
  double Wy  = w.dot(tmp.basisDerivs_v);
  double Wxx = w.dot(tmp.basisDerivs_uu);
  double Wyy = w.dot(tmp.basisDerivs_vv);
  double Wxy = w.dot(tmp.basisDerivs_uv);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i)
  {
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;

    double H1 = (tmp.basisDerivs_u[i]*W - tmp.basisValues[i]*Wx);
    double H2 = (tmp.basisDerivs_v[i]*W - tmp.basisValues[i]*Wy);
    bas.basisDerivs_u[i] = H1*w[i]/(W*W);
    bas.basisDerivs_v[i] = H2*w[i]/(W*W);

    double H1x = tmp.basisDerivs_uu[i]*W - tmp.basisValues[i]*Wxx;
    double H2y = tmp.basisDerivs_vv[i]*W - tmp.basisValues[i]*Wyy;
    double H1y = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_u[i]*Wy - tmp.basisDerivs_v[i]*Wx;
    bas.basisDerivs_uu[i] = (H1x*W - 2.0*H1*Wx)*w[i]/(W*W*W);
    bas.basisDerivs_vv[i] = (H2y*W - 2.0*H2*Wy)*w[i]/(W*W*W);
    bas.basisDerivs_uv[i] = (H1y*W - 2.0*H1*Wy)*w[i]/(W*W*W);
  }
}


void ASMu2Dnurbs::computeBasis (double u, double v,
                                Go::BasisDerivsSf3& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  PROFILE3("ASMu2Dn::compBasis(3)");

  const LR::Element* el = lrspline->getElement(iel);

  Go::BasisDerivsSf3 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W    = w.dot(tmp.basisValues);
  double Wx   = w.dot(tmp.basisDerivs_u);
  double Wy   = w.dot(tmp.basisDerivs_v);
  double Wxx  = w.dot(tmp.basisDerivs_uu);
  double Wyy  = w.dot(tmp.basisDerivs_vv);
  double Wxy  = w.dot(tmp.basisDerivs_uv);
  double Wxxx = w.dot(tmp.basisDerivs_uuu);
  double Wyyy = w.dot(tmp.basisDerivs_vvv);
  double Wxxy = w.dot(tmp.basisDerivs_uuv);
  double Wxyy = w.dot(tmp.basisDerivs_uvv);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
  {
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;

    double H1 = (tmp.basisDerivs_u[i]*W - tmp.basisValues[i]*Wx);
    double H2 = (tmp.basisDerivs_v[i]*W - tmp.basisValues[i]*Wy);
    bas.basisDerivs_u[i] = H1*w[i]/(W*W);
    bas.basisDerivs_v[i] = H2*w[i]/(W*W);

    double H1x = tmp.basisDerivs_uu[i]*W - tmp.basisValues[i]*Wxx;
    double H2y = tmp.basisDerivs_vv[i]*W - tmp.basisValues[i]*Wyy;
    double H1y = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_u[i]*Wy - tmp.basisDerivs_v[i]*Wx;
    double H2x = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_v[i]*Wx - tmp.basisDerivs_u[i]*Wy;
    double G1  = H1x*W - 2.0*H1*Wx;
    double G2  = H2y*W - 2.0*H2*Wy;
    bas.basisDerivs_uu[i] = G1*w[i]/(W*W*W);
    bas.basisDerivs_vv[i] = G2*w[i]/(W*W*W);
    bas.basisDerivs_uv[i] = (H1y*W - 2.0*H1*Wy)*w[i]/(W*W*W);

    double H1xx = tmp.basisDerivs_uuu[i]*W + tmp.basisDerivs_uu[i]*Wx
                - tmp.basisDerivs_u[i]*Wxx - tmp.basisValues[i]*Wxxx;
    double H2yy = tmp.basisDerivs_vvv[i]*W + tmp.basisDerivs_vv[i]*Wy
                - tmp.basisDerivs_v[i]*Wyy - tmp.basisValues[i]*Wyyy;
    double H1xy = tmp.basisDerivs_uuv[i]*W + tmp.basisDerivs_uu[i]*Wy
                - tmp.basisDerivs_v[i]*Wxx - tmp.basisValues[i]*Wxxy;
    double H2xy = tmp.basisDerivs_uvv[i]*W + tmp.basisDerivs_vv[i]*Wx
                - tmp.basisDerivs_u[i]*Wyy - tmp.basisValues[i]*Wxyy;

    double G1x = H1xx*W + H1x*Wx - 2.0*H1x*Wx - 2.0*H1*Wxx;
    double G1y = H1xy*W + H1x*Wy - 2.0*H1y*Wx - 2.0*H1*Wxy;
    double G2x = H2xy*W + H2y*Wx - 2.0*H2x*Wy - 2.0*H2*Wxy;
    double G2y = H2yy*W + H2y*Wy - 2.0*H2y*Wy - 2.0*H2*Wyy;

    bas.basisDerivs_uuu[i] = (G1x*W - 3.0*G1*Wx)*w[i]/(W*W*W*W);
    bas.basisDerivs_vvv[i] = (G2y*W - 3.0*G2*Wy)*w[i]/(W*W*W*W);
    bas.basisDerivs_uuv[i] = (G1y*W - 3.0*G1*Wy)*w[i]/(W*W*W*W);
    bas.basisDerivs_uvv[i] = (G2x*W - 3.0*G2*Wx)*w[i]/(W*W*W*W);
  }
}


LR::LRSplineSurface* ASMu2Dnurbs::createLRfromTensor ()
{
  // Creates a dim+1 dimensional LRSplineSurface from a tensor NURBS surface.
  auto&& createLRnurbs = [](const Go::SplineSurface* srf)
  {
    return new LR::LRSplineSurface(srf->numCoefs_u(), srf->numCoefs_v(),
                                   srf->order_u(), srf->order_v(),
                                   srf->basis_u().begin(),
                                   srf->basis_v().begin(),
                                   srf->rcoefs_begin(),
                                   srf->dimension()+1);
  };

  if (tensorspline)
  {
    if ((noNurbs = !tensorspline->rational()))
      lrspline.reset(new LR::LRSplineSurface(tensorspline));
    else
      lrspline.reset(createLRnurbs(tensorspline));
    delete tensorspline;
    tensorspline = nullptr;
  }

  return lrspline.get();
}
