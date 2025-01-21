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

#include "ASMu2D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Profiler.h"


bool ASMu2D::evaluateBasisNurbs (int iel, FiniteElement& fe,
                                 int derivs) const
{
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
    C.multiply(B,fe.N); // fe.N = C*B;
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


void ASMu2D::computeBasisNurbs (double u, double v,
                                Go::BasisPtsSf& bas, int iel,
                                const LR::LRSplineSurface& spline)
{
  PROFILE3("ASMu2D::compBasisN(0)");

  const LR::Element* el = spline.getElement(iel);

  Go::BasisPtsSf tmp;
  spline.computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W = w.dot(tmp.basisValues);

  bas.preparePts(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;
}


void ASMu2D::computeBasisNurbs (double u, double v,
                                Go::BasisDerivsSf& bas, int iel,
                                const LR::LRSplineSurface& spline)
{
  PROFILE3("ASMu2D::compBasisN(1)");

  const LR::Element* el = spline.getElement(iel);

  Go::BasisDerivsSf tmp;
  spline.computeBasis(u,v,tmp,iel);
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


void ASMu2D::computeBasisNurbs (double u, double v,
                                Go::BasisDerivsSf2& bas, int iel,
                                const LR::LRSplineSurface& spline)
{
  PROFILE3("ASMu2D::compBasisN(2)");

  const LR::Element* el = spline.getElement(iel);

  Go::BasisDerivsSf2 tmp;
  spline.computeBasis(u,v,tmp,iel);
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


void ASMu2D::computeBasisNurbs (double u, double v,
                                Go::BasisDerivsSf3& bas, int iel,
                                const LR::LRSplineSurface& spline)
{
  PROFILE3("ASMu2D::compBasisN(3)");

  const LR::Element* el = spline.getElement(iel);

  Go::BasisDerivsSf3 tmp;
  spline.computeBasis(u,v,tmp,iel);
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

void ASMu2D::writePostscriptElementsNurbs (std::shared_ptr<LR::LRSplineSurface> mesh,
                                           std::ostream& out, bool close,
                                           int nu, int nv)
{
  // get date
  time_t t = time(0);
  tm* lt = localtime(&t);
  char date[128];
  sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

  if (mesh != lrspline) {
    mesh.swap(lrspline);
    geomB = lrspline;
  }

  // This is done unconditionally as it will not have been performed
  // for the refined mesh. It will be on the next adaptive cycle though,
  // so we do it here instead of in general to avoid paying the price if not
  // outputting eps files.
  this->generateBezierBasis();
  this->generateBezierExtraction();

  // get bounding box (max/min of the control points)
  std::array<double,2> x{1e7, -1e7};
  std::array<double,2> y{1e7, -1e7};
  for (const LR::Basisfunction* b : lrspline->getAllBasisfunctions()) {
    std::vector<double>::const_iterator cp = b->cp();
    x[0] = (cp[0] < x[0]) ? cp[0] : x[0];
    x[1] = (cp[0] > x[1]) ? cp[0] : x[1];
    y[0] = (cp[1] < y[0]) ? cp[1] : y[0];
    y[1] = (cp[1] > y[1]) ? cp[1] : y[1];
  }

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double scale = dx > dy ? 1000.0/dx : 1000.0/dy;

  int xmin = (x[0] - dx/20.0)*scale;
  int ymin = (y[0] - dy/20.0)*scale;
  int xmax = (x[1] + dx/20.0)*scale;
  int ymax = (y[1] + dy/20.0)*scale;

  // print eps header
  out << "%!PS-Adobe-3.0 EPSF-3.0\n";
  out << "%%Creator: LRSplineHelpers.cpp object\n";
  out << "%%Title: LR-spline physical domain\n";
  out << "%%CreationDate: " << date << std::endl;
  out << "%%Origin: 0 0\n";
  out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
  out << "0 setgray\n";
  out << "1 setlinewidth\n";

  for (const LR::Element* el : lrspline->getAllElements()) {
    double umin = el->umin();
    double umax = el->umax();
    double vmin = el->vmin();
    double vmax = el->vmax();

    Vec3 pt;
    double u[2] = {umin, vmin};
    this->evalPoint(el->getId(), u, pt);

    out << "newpath\n";
    out <<  pt[0]*scale << " " << pt[1]*scale << " moveto\n";
    for (int i = 1; i < nu; i++) { // SOUTH side
      u[0] = umin + (umax-umin)*i/(nu-1);
      u[1] = vmin;
      this->evalPoint(el->getId(), u, pt);
      out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
    }
    for (int i = 1; i < nv; i++) { // EAST side
      u[0] = umax;
      u[1] = vmin + (vmax-vmin)*i/(nv-1);
      this->evalPoint(el->getId(), u, pt);
      out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
    }
    for (int i = nu-1 ; i-- > 0; ) { // NORTH side
      u[0] = umin + (umax-umin)*i/(nu-1);
      u[1] = vmax;
      this->evalPoint(el->getId(), u, pt);
      out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
    }
    for (int i = nv-1; i-- > 1; ) { // WEST side
      u[0] = umin;
      u[1] = vmin + (vmax-vmin)*i/(nv-1);
      this->evalPoint(el->getId(), u, pt);
      out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
    }
    out << "closepath\n";
    out << "stroke\n";
  }

  if (close)
    out << "%%EOF\n";

  if (mesh != lrspline) {
    mesh.swap(lrspline);
    geomB = lrspline;
    this->generateBezierBasis();
    this->generateBezierExtraction();
  }
}


void ASMu2D::writePostscriptMeshWithControlPointsNurbs (std::shared_ptr<LR::LRSplineSurface> mesh,
                                                        std::ostream& out, int nu, int nv)
{
  this->writePostscriptElementsNurbs(mesh, out, false, nu, nv);

  // get bounding box (max/min of the control points)
  std::array<double,2> x{1e7, -1e7};
  std::array<double,2> y{1e7, -1e7};
  for (const LR::Basisfunction* b : mesh->getAllBasisfunctions()) {
    std::vector<double>::const_iterator cp = b->cp();
    x[0] = (cp[0] < x[0]) ? cp[0] : x[0];
    x[1] = (cp[0] > x[1]) ? cp[0] : x[1];
    y[0] = (cp[1] < y[0]) ? cp[1] : y[0];
    y[1] = (cp[1] > y[1]) ? cp[1] : y[1];
  }

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double scale = dx > dy ? 1000.0/dx : 1000.0/dy;

  double circleSize = 15.0;

  // create the ellipse function
  out << "/ellipse {\n";
  out << "/endangle exch def\n";
  out << "/startangle exch def\n";
  out << "/yrad exch def\n";
  out << "/xrad exch def\n";
  out << "/y exch def\n";
  out << "/x exch def\n";
  out << "/savematrix matrix currentmatrix def\n";
  out << "x y translate\n";
  out << "xrad yrad scale\n";
  out << "0 0 1 startangle endangle arc\n";
  out << "savematrix setmatrix\n";
  out << "} def\n";

  // load the font to use
  out << "/Times-Roman findfont\n";
  out << "24 scalefont\n";
  out << "setfont\n";

  int i = -1;
  for (const LR::Basisfunction* b : mesh->getAllBasisfunctions()) {
    i++;
    double cp_x = b->cp(0);
    double cp_y = b->cp(1);
    // move C^{-1} text on internal functions
    int textX   = ((*b)[0][1] == (*b)[0][mesh->order(0)]) ? -2 : 1;
    int textY   = ((*b)[1][1] == (*b)[1][mesh->order(1)]) ? -2 : 1;
    // move text on edge functions
    if ((*b)[0][1] == mesh->endparam(0))
      textX = 1;
    else if ((*b)[0][mesh->order(0)-1] == mesh->startparam(0))
      textX = -2;
    if ((*b)[1][1] == mesh->endparam(1))
      textY = 1;
    else if ((*b)[1][mesh->order(1)-1] == mesh->startparam(1))
      textY = -2;

    out << "newpath\n";
    out << "0.45 0.45 0.45 setrgbcolor \n";
    out << cp_x*scale << " " << cp_y*scale << " " << circleSize << " " << circleSize << " 0 360 ellipse\n";
    out << "closepath fill\n";
    out << "0 setgray\n";
    out << cp_x*scale << " " << cp_y*scale << " " << circleSize << " " << circleSize << " 0 360 ellipse\n";
    out << "closepath stroke\n";
    out << "\n";
    out << "newpath\n";
    out << cp_x*scale + textX*circleSize << " " << cp_y*scale + textY*circleSize << " moveto\n";
    out << "(" << i << ") show\n";
    out << "\n";
  }
  out << "%%EOF\n";
}
