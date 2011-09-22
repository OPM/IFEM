// $Id: ASMu2D.C 1108 2011-08-23 14:51:03Z kjetijo $
//==============================================================================
//!
//! \file ASMu2D.C
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"
#include <fstream>

#include "ASMu2D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include <ctype.h>

ASMu2D::ASMu2D (const char* fName, unsigned char n_s, unsigned char n_f)
	: ASMunstruct(2,n_s,n_f), lrspline(0)
{
	if (fName)
	{
		std::cout <<"\nReading patch file "<< fName << std::endl;
		std::ifstream is(fName);
		if (!is.good())
			std::cerr <<" *** ASMu2D: Failure opening patch file"<< std::endl;
		else
			this->read(is);
	}
}


ASMu2D::ASMu2D (std::istream& is, unsigned char n_s, unsigned char n_f)
	: ASMunstruct(2,n_s,n_f), lrspline(0)
{
	this->read(is);
}


bool ASMu2D::read (std::istream& is)
{
	if (lrspline) delete lrspline;

	// read inputfile as either an LRSpline file directly or a tensor product B-spline and convert
	char firstline[256];
	is.getline(firstline, 256);
	if(strncmp(firstline, "# LRSPLINE", 10) == 0) {
		lrspline = new LR::LRSplineSurface();
		is >> *lrspline;
	} else { // probably a SplineSurface, so we'll read that and convert
		Go::SplineSurface *surf = new Go::SplineSurface();
		is >> *surf;
		lrspline = new LR::LRSplineSurface(surf);
	}

	// Eat white-space characters to see if there is more data to read
	char c;
	while (is.get(c))
		if (!isspace(c))
		{
			is.putback(c);
			break;
		}

	if (!is.good() && !is.eof())
	{
		std::cerr <<" *** ASMu2D::read: Failure reading spline data"<< std::endl;
		delete lrspline;
		lrspline = 0;
		return false;
	}
	else if (lrspline->dimension() < 2)
	{
		std::cerr <<" *** ASMu2D::read: Invalid spline lrsplineace patch, dim="
		          << lrspline->dimension() << std::endl;
		delete lrspline;
		lrspline = 0;
		return false;
	}
	else if (lrspline->dimension() < nsd)
	{
		std::cout <<"  ** ASMu2D::read: The dimension of this lrsplineace patch "
		          << lrspline->dimension() <<" is less than nsd="<< nsd
		          <<".\n                   Resetting nsd to "<< lrspline->dimension()
		          <<" for this patch."<< std::endl;
		nsd = lrspline->dimension();
	}

	geo = lrspline;
	return true;
}


bool ASMu2D::write (std::ostream& os, int) const
{
	if (!lrspline) return false;

	os << *lrspline;

	return os.good();
}


void ASMu2D::clear ()
{
	// Erase spline data
	if (lrspline) delete lrspline;
	lrspline = 0;
	geo = 0;

	// Erase the FE data
	ASMbase::clear();
}

bool ASMu2D::cornerRefine (int minBasisfunctions) 
{
	if (!lrspline ) return false;

	double h = 1.0;
	int nBasis = lrspline->nBasisFunctions();
	double unif_step_h = 1.0 / ((minBasisfunctions - nBasis) / 3.0 + 1.0);
	while(lrspline->nBasisFunctions() < minBasisfunctions) {
		lrspline->insert_const_u_edge(h-unif_step_h, 0, h);
		lrspline->insert_const_v_edge(h-unif_step_h, 0, h);
		h -= unif_step_h;
	}
	std::ofstream meshFile;
	meshFile.open("mesh.eps");
	lrspline->writePostscriptMesh(meshFile);
	meshFile.close();
	return true;
}

bool ASMu2D::diagonalRefine (int minBasisfunctions) 
{
	if (!lrspline ) return false;

	double end1 = lrspline->endparam_u();
	double end2 = lrspline->endparam_v();
	double h = 1.0;
	int iter = 0;
	double u = h/2.0;
	double v = h/2.0;
	while(lrspline->nBasisFunctions() < minBasisfunctions) {
		lrspline->insert_const_u_edge(u, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end2) ? end2 : (iter+2)*h);
		lrspline->insert_const_v_edge(v, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end1) ? end1 : (iter+2)*h);
		u += h;
		v += h;
		iter++;
		if( u>end1 ) {
			h /= 2.0;
			iter = 0;
			u = h/2.0;
			v = h/2.0;
		}
	}
	std::ofstream meshFile;
	meshFile.open("mesh.eps");
	lrspline->writePostscriptMesh(meshFile);
	meshFile.close();

	return true;
}

bool ASMu2D::uniformRefine (int minBasisfunctions) 
{
	if (!lrspline ) return false;

	double end1 = lrspline->endparam_u();
	double end2 = lrspline->endparam_v();
	double h = 1.0;
	bool step_u = true;
	double u = h/2.0;
	double v = h/2.0;
	while(lrspline->nBasisFunctions() < minBasisfunctions) {
		if(step_u) {
			lrspline->insert_const_u_edge(u, 0, end2);
			u += h;
			if(u > end1) {
				step_u = !step_u;
				u = h/4.0;
			}
		} else {
			lrspline->insert_const_v_edge(v, 0, end1);
			v += h;
			if(v > end2) {
				step_u = !step_u;
				v = h/4.0;
				h /= 2.0;
			}
		}
	}
	std::ofstream meshFile;
	meshFile.open("mesh.eps");
	lrspline->writePostscriptMesh(meshFile);
	meshFile.close();
	return true;
}


bool ASMu2D::generateFEMTopology ()
{
	if (!lrspline) return false;

	const int nBasis    = lrspline->nBasisFunctions();
	const int nElements = lrspline->nElements();

	const int p1 = lrspline->order_u();
	const int p2 = lrspline->order_v();

	// Consistency checks, just to be fool-proof
	if (nBasis < 4)         return false;
	if (p1 <  1 || p1 <  1) return false;
	if (p1*p2 > nBasis)     return false;

	MLGE.resize(nElements,0);
	MLGN.resize(nBasis);
	MNPC.resize(MLGE.size());
	lrspline->generateIDs();

	int iel = 0;
	std::vector<LR::Element*>::iterator el_it = lrspline->elementBegin();
	std::vector<LR::Basisfunction*>::iterator b_it;
	for(iel=0; iel<nElements; iel++, el_it++)
	{
		int nSupportFunctions =  (*el_it)->nBasisFunctions();
		MLGE[iel] = ++gEl; // global element number over all patches
		MNPC[iel].resize(nSupportFunctions, 0);

		b_it = (*el_it)->supportBegin();
		for(int lnod=0; lnod<nSupportFunctions; lnod++, b_it++)
			MNPC[iel][lnod] = (*b_it)->getId();
	}

	for(size_t i=0; i<MLGN.size(); i++)
		MLGN[i] = i+1;

	gNod += nBasis;

	return true;
}

// this is connecting multiple patches and handling deformed geometries.
// We'll deal with at a later time, for now we only allow single patch models

#if 0

int ASMu2D::Edge::next ()
{
	int ret = icnod;
	icnod += incr;

	return ret;
}


int ASMu2D::BlockNodes::next ()
{
	int ret = iinod;
	iinod += inc[0];

	if (++indxI >= nnodI-1)
	{
		indxI = 1;
		iinod += inc[1] - inc[0]*(nnodI-2);
	}

	return ret;
}


bool ASMu2D::assignNodeNumbers ()
{
	for(size_t i=0; i<MLGN.size(); i++)
		MLGN[i] = i+1;
	return true;
} 

bool ASMu2D::connectPatch (int edge, ASMu2D& neighbor, int nedge, bool revers)
{
	return this->connectBasis(edge,neighbor,nedge,revers);
}


bool ASMu2D::connectBasis (int edge, ASMu2D& neighbor, int nedge, bool revers,
                           int basis, int slave, int master)
{
	// Set up the slave node numbers for this surface patch

	int n1, n2;
	if (!this->getSize(n1,n2,basis)) return false;
	int node = slave+1, i1 = 0;

	switch (edge)
		{
		case 2: // Positive I-direction
			node += n1-1;
		case 1: // Negative I-direction
			i1 = n1;
			n1 = n2;
			break;

		case 4: // Positive J-direction
			node += n1*(n2-1);
		case 3: // Negative J-direction
			i1 = 1;
			break;

		default:
			std::cerr <<" *** ASMu2D::connectPatch: Invalid slave edge "
			          << edge << std::endl;
			return false;
		}

	int i;
	IntVec slaveNodes(n1,0);
	for (i = 0; i < n1; i++, node += i1)
		slaveNodes[i] = node;

	// Set up the master node numbers for the neighboring surface patch

	if (!neighbor.getSize(n1,n2,basis)) return false;
	node = master+1; i1 = 0;

	switch (nedge)
		{
		case 2: // Positive I-direction
			node += n1-1;
		case 1: // Negative I-direction
			i1 = n1;
			n1 = n2;
			break;

		case 4: // Positive J-direction
			node += n1*(n2-1);
		case 3: // Negative J-direction
			i1 = 1;
			break;

		default:
			std::cerr <<" *** ASMu2D::connectPatch: Invalid master edge "
			          << nedge << std::endl;
			return false;
		}

	if (n1 != (int)slaveNodes.size())
	{
		std::cerr <<" *** ASMu2D::connectPatch: Non-matching edges, sizes "
		          << n1 <<" and "<< slaveNodes.size() << std::endl;
		return false;
	}

	const double xtol = 1.0e-4;
	for (i = 0; i < n1; i++, node += i1)
	{
		int k = revers ? n1-i-1 : i;
		if (!neighbor.getCoord(node).equal(this->getCoord(slaveNodes[k]),xtol))
		{
			std::cerr <<" *** ASMu2D::connectPatch: Non-matching nodes "
			          << node <<": "<< neighbor.getCoord(node)
			          <<"\n                                          and "
			          << slaveNodes[k] <<": "<< this->getCoord(slaveNodes[k])
			          << std::endl;
			return false;
		}
		else
			ASMbase::collapseNodes(neighbor.MLGN[node-1],MLGN[slaveNodes[k]-1]);
	}

	return true;
}


void ASMu2D::closeEdges (int dir, int basis, int master)
{
	int n1, n2;
	if (basis < 1) basis = 1;
	if (!this->getSize(n1,n2,basis)) return;

	switch (dir)
		{
		case 1: // Edges are closed in I-direction
			for (int i2 = 1; i2 <= n2; i2++, master += n1)
				this->makePeriodic(master,master+n1-1);
			break;

		case 2: // Edges are closed in J-direction
			for (int i1 = 1; i1 <= n1; i1++, master++)
				this->makePeriodic(master,master+n1*(n2-1));
			break;
		}
}

#endif


void ASMu2D::constrainEdge (int dir, int dof, int code)
{
	std::vector<LR::Basisfunction*> edgeFunctions;
	switch (dir)
		{
		case  1: // Right edge (positive I-direction)
			lrspline->getEdgeFunctions(edgeFunctions, LR::EAST);
			break;
		case -1: // Left edge (negative I-direction)
			lrspline->getEdgeFunctions(edgeFunctions, LR::WEST);
			break;
		case  2: // Back edge (positive J-direction)
			lrspline->getEdgeFunctions(edgeFunctions, LR::NORTH);
			break;
		case -2: // Front edge (negative J-direction)
			lrspline->getEdgeFunctions(edgeFunctions, LR::SOUTH);
			break;
		}
	for(size_t i=0; i<edgeFunctions.size(); i++) 
		this->prescribe(edgeFunctions[i]->getId()+1,dof,code);
}


void ASMu2D::constrainCorner (int I, int J, int dof, int code)
{
	std::vector<LR::Basisfunction*> edgeFunctions;
	if     (I == 0 && J == 0)
		lrspline->getEdgeFunctions(edgeFunctions, LR::SOUTH_WEST);
	else if(I >  0 && J == 0)
		lrspline->getEdgeFunctions(edgeFunctions, LR::SOUTH_EAST);
	else if(I == 0 && J >  0)
		lrspline->getEdgeFunctions(edgeFunctions, LR::NORTH_WEST);
	else 
		lrspline->getEdgeFunctions(edgeFunctions, LR::NORTH_EAST);

#ifdef SP_DEBUG
	if(edgeFunctions.size() != 1) {
		std::cerr <<" *** ASMu2D::constrainCorner: more than one corner"
		          <<" returned from lrspline->getEdgeFunctions()" << std::endl;
		return;
	}
#endif
	
	this->prescribe(edgeFunctions[0]->getId()+1,dof,code);
}


// Hopefully we don't have to constrain non-corner singlenodes inside patches
#if 0
void ASMu2D::constrainNode (double xi, double eta, int dof, int code)
{
	if (xi  < 0.0 || xi  > 1.0) return;
	if (eta < 0.0 || eta > 1.0) return;

	int n1, n2;
	if (!this->getSize(n1,n2,1)) return;

	int node = 1;
	if (xi  > 0.0) node += int(0.5+(n1-1)*xi);
	if (eta > 0.0) node += n1*int(0.5+(n2-1)*eta);

	this->prescribe(node,dof,code);
}
#endif


#define DERR -999.99

double ASMu2D::getParametricArea (int iel) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu2D::getParametricArea: Element index "<< iel
		          <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return DERR;
	}
#endif
	if (MNPC[iel-1].empty())
		return 0.0;

	return lrspline->getElement(iel-1)->area();
}


double ASMu2D::getParametricLength (int iel, int dir) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu2D::getParametricLength: Element index "<< iel
		          <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return DERR;
	}
#endif
	if (MNPC[iel-1].empty())
		return 0.0;

	LR::Element *el = lrspline->getElement(iel-1);
	switch (dir)
		{
		case 1: return el->vmax() - el->vmin();
		case 2: return el->umax() - el->umin();
		}

	std::cerr <<" *** ASMu2D::getParametricLength: Invalid edge direction "
	          << dir << std::endl;
	return DERR;
}


bool ASMu2D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu2D::getElementCoordinates: Element index "<< iel
		          <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return false;
	}
#endif

	const IntVec& mnpc = MNPC[iel-1];
	X.resize(nsd,mnpc.size());

	std::vector<LR::Basisfunction*>::iterator bit = lrspline->getElement(iel-1)->supportBegin();
	for (size_t n = 0; n < mnpc.size(); n++, bit++)
		for (size_t i = 0; i < nsd; i++)
			X(i+1,n+1) = (*bit)->controlpoint_[i];

#if SP_DEBUG > 2
	std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
	return true;
}


void ASMu2D::getNodalCoordinates (Matrix& X) const
{
	const int nBasis = lrspline->nBasisFunctions();
	X.resize(nsd,nBasis);

	std::vector<LR::Basisfunction*>::iterator bit = lrspline->basisBegin();

	for(int inod=0; inod<nBasis; inod++, bit++)
		for(size_t i=0; i<nsd; i++)
			X(i+1,inod+1) = (*bit)->controlpoint_[i];
}


Vec3 ASMu2D::getCoord (size_t inod) const
{
	Vec3 X;
	LR::Basisfunction* basis = lrspline->getBasisfunction(inod-1);
	for (unsigned char i = 0; i < nsd; i++)
		X[i] = basis->controlpoint_[i];

	return X;
}


void ASMu2D::getGaussPointParameters (Vector& uGP, int dir, int nGauss,
                                      int iel, const double* xi) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu2D::getGaussPointParameters: Element index "<< iel
		          <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return;
	}
#endif
	LR::Element* el = lrspline->getElement(iel-1);

	uGP.resize(nGauss);
	double ustart = (dir==0) ? el->umin() : el->vmin();
	double ustop  = (dir==0) ? el->umax() : el->vmax();

	for (int i = 1; i <= nGauss; i++)
		uGP(i) = 0.5*((ustop-ustart)*xi[i-1] + ustop+ustart);
}


/*!
	\brief Computes the characteristic element length from nodal coordinates.
*/

#if 0
static double getElmSize (int p1, int p2, const Matrix& X)
{
	int n = X.rows();
	int i, j, id1, id2;
	double value, v1, h = 1.0e12;

	// Y-direction
	for (i = 1; i <= p1; i++)
	{
		id1 = i;
		id2 = id1 + (p2-1)*p1;
		value = 0.0;
		for (j = 1; j <= n; j++)
		{
			v1 = X(j,id2) - X(j,id1);
			value += v1*v1;
		}
		if (value < h) h = value;
	}

	// X-direction
	for (j = 0; j < p2; j++)
	{
		id1 = j*p1 + 1;
		id2 = id1 + p1 - 1;
		value = 0.0;
		for (i = 1; i <= n; i++)
		{
			v1 = X(i,id2) - X(i,id1);
			value += v1*v1;
		}
		if (value < h) h = value;
	}

	return sqrt(h);
}
#endif


void ASMu2D::extractBasis (const Go::BasisDerivsSf& spline,
                           Vector& N, Matrix& dNdu)
{
	dNdu.resize(N.size(),2);

	size_t jp, n = 1;
	for (jp = 0; jp < N.size(); jp++, n++)
	{
		 N  (n)   = spline.basisValues[jp];
		dNdu(n,1) = spline.basisDerivs_u[jp];
		dNdu(n,2) = spline.basisDerivs_v[jp];
	}
}


void ASMu2D::extractBasis (const Go::BasisDerivsSf2& spline,
                           Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2)
{
	 dNdu .resize(N.size(),2);
	d2Ndu2.resize(N.size(),2,2);

	size_t jp, n = 1;
	for (jp = 0; jp < N.size(); jp++, n++)
	{
			N   (n)     = spline.basisValues[jp];
		 dNdu (n,1)   = spline.basisDerivs_u[jp];
		 dNdu (n,2)   = spline.basisDerivs_v[jp];
		d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
		d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
		d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
	}
}


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const LintegralVec& locInt)
{
	if (!lrspline) return true; // silently ignore empty patches

	PROFILE2("ASMu2D::integrate(I)");

	// Get Gaussian quadrature points and weights
	const double* xg = GaussQuadrature::getCoord(nGauss);
	const double* wg = GaussQuadrature::getWeight(nGauss);
	if (!xg || !wg) return false;

	// Get the reduced integration quadrature points, if needed
	const double* xr = 0;
	int nRed = integrand.getIntegrandType() - 10;
	if (nRed < 1)
		nRed = nRed < 0 ? nGauss : 0;
	else if (!(xr = GaussQuadrature::getCoord(nRed)))
		return false;

	Matrix   dNdu, Xnod, Jac;
	Matrix3D d2Ndu2, Hess;
	Vec4     X;


	// === Assembly loop over all elements in the patch ==========================

	std::vector<LR::Element*>::iterator el;
	int iel = 1;
	for(el = lrspline->elementBegin(); el!=lrspline->elementEnd(); el++, iel++)
	{
		FiniteElement fe(MNPC[iel-1].size());
		fe.iel = MLGE[iel-1];
		if (fe.iel < 1) continue; // zero-area element

		// Get element area in the parameter space
		double dA = this->getParametricArea(iel);
		if (dA < 0.0) return false; // topology error (probably logic error)

		// Set up control point (nodal) coordinates for current element
		if (!this->getElementCoordinates(Xnod,iel)) return false;

		// Compute parameter values of the Gauss points over this element
		Vector gpar[2], redpar[2];
		for (int d = 0; d < 2; d++)
		{
			this->getGaussPointParameters(gpar[d],d,nGauss,iel,xg);
			if (integrand.getIntegrandType() > 10)
				this->getGaussPointParameters(redpar[d],d,nRed,iel,xr);
		}


#if 0
		// Compute characteristic element length, if needed
		if (integrand.getIntegrandType() == 2)
			fe.h = getElmSize(p1,p2,Xnod);

		else if (integrand.getIntegrandType() == 3)
		{
		// --- Compute average value of basis functions over the element -------

			fe.Navg.resize(p1*p2,true);
			double area = 0.0;
			int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
			for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
				for (int i = 0; i < nGauss; i++, ip++)
				{
				// Fetch basis function derivatives at current integration point
					extractBasis(spline[ip],fe.N,dNdu);

					// Compute Jacobian determinant of coordinate mapping
					// and multiply by weight of current integration point
					double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
					double weight = 0.25*dA*wg[i]*wg[j];

					// Numerical quadrature
					fe.Navg.add(fe.N,detJac*weight);
					area += detJac*weight;
				}

			// Divide by element area
			fe.Navg /= area;
		}

		else if (integrand.getIntegrandType() == 4)
		{
		// Compute the element center
			Go::Point X0;
			double u0 = 0.5*(gpar[0](1,i1-p1+1) + gpar[0](nGauss,i1-p1+1));
			double v0 = 0.5*(gpar[1](1,i2-p2+1) + gpar[1](nGauss,i2-p2+1));
			lrspline->point(X0,u0,v0);
			for (unsigned char i = 0; i < nsd; i++)
				X[i] = X0[i];
		}
#endif

		// Initialize element quantities
		if (!integrand.initElement(MNPC[iel-1],X,nRed*nRed)) return false;

		// Caution: Unless locInt is empty, we assume it points to an array of
		// LocalIntegral pointers, of length at least the number of elements in
		// the model (as defined by the highest number in the MLGE array).
		// If the array is shorter than this, expect a segmentation fault.
		LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


#if 0
		if (integrand.getIntegrandType() > 10)
		{
		// --- Selective reduced integration loop ------------------------------

			int ip = ((i2-p2)*nRed*nel1 + i1-p1)*nRed;
			for (int j = 0; j < nRed; j++, ip += nRed*(nel1-1))
				for (int i = 0; i < nRed; i++, ip++)
				{
					// Local element coordinates of current integration point
					fe.xi  = xr[i];
					fe.eta = xr[j];

					// Parameter values of current integration point
					fe.u = redpar[0](i+1,i1-p1+1);
					fe.v = redpar[1](j+1,i2-p2+1);

					// Fetch basis function derivatives at current point
					extractBasis(splineRed[ip],fe.N,dNdu);

					// Compute Jacobian inverse and derivatives
					fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

					// Compute the reduced integration terms of the integrand
					if (!integrand.reducedInt(fe))
					  return false;
				}
		}
#endif


		// --- Integration loop over all Gauss points in each direction ----------

		for (int j = 0; j < nGauss; j++)
			for (int i = 0; i < nGauss; i++)
			{
				// Local element coordinates of current integration point
				fe.xi  = xg[i];
				fe.eta = xg[j];

				// Parameter values of current integration point
				fe.u = gpar[0](i+1);
				fe.v = gpar[1](j+1);

				// compute basis function values at integration point
				Go::BasisDerivsSf  spline;
				Go::BasisDerivsSf2 spline2;
				Go::BasisDerivsSf  splineRed;
// I'll impement second order differentiation in LRSplines soon - promise
#if 0
				if (integrand.getIntegrandType() == 2)
					lrspline->computeBasis(fe.u,fe.v,spline2, iel);
				else
#endif
					lrspline->computeBasis(fe.u,fe.v,spline, iel-1);
				if (integrand.getIntegrandType() > 10)
					lrspline->computeBasis(fe.u,fe.v,splineRed, iel-1);
		
				// Fetch basis function derivatives at current integration point
				if (integrand.getIntegrandType() == 2)
					extractBasis(spline2,fe.N,dNdu,d2Ndu2);
				else
					extractBasis(spline,fe.N,dNdu);

#if SP_DEBUG > 4
   				std::cout <<"\nBasis functions at a integration point ";
   				std::cout <<" : (u,v) = "<< spline.param[0] <<" "<< spline.param[1]
   				          <<"  left_idx = "<< spline.left_idx[0] <<" "<< spline.left_idx[1] << std::endl;
   				for (size_t ii = 0; ii < spline.basisValues.size(); ii++)
   					std::cout << 1+ii <<'\t'<< spline.basisValues[ii] <<'\t'
					          << spline.basisDerivs_u[ii] <<'\t'<< spline.basisDerivs_v[ii] << std::endl;
#endif


				// Compute Jacobian inverse of coordinate mapping and derivatives
				fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
				if (fe.detJxW == 0.0) continue; // skip singular points
		
				// Compute Hessian of coordinate mapping and 2nd order derivatives
				if (integrand.getIntegrandType() == 2)
					if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
					  return false;

#if SP_DEBUG > 4
				std::cout  <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX << std::endl;
#endif

				// Cartesian coordinates of current integration point
				X = Xnod * fe.N;
				X.t = time.t;

				// Evaluate the integrand and accumulate element contributions
				fe.detJxW *= 0.25*dA*wg[i]*wg[j];
				if (!integrand.evalInt(elmInt,fe,time,X))
				return false;
			}

		// Finalize the element quantities
		if (!integrand.finalizeElement(elmInt,time))
			return false;

		// Assembly of global system integral
		if (!glInt.assemble(elmInt,fe.iel))
			return false;
	}

	return true;
}


bool ASMu2D::integrate (Integrand& integrand, int lIndex,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const LintegralVec& locInt)
{
	if (!lrspline) return true; // silently ignore empty patches

#if 0
	PROFILE2("ASMu2D::integrate(B)");

	// Get Gaussian quadrature points and weights
	const double* xg = GaussQuadrature::getCoord(nGauss);
	const double* wg = GaussQuadrature::getWeight(nGauss);
	if (!xg || !wg) return false;

	// Find the parametric direction of the edge normal {-2,-1, 1, 2}
	const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

	const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
	const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

	// Compute parameter values of the Gauss points along the whole patch edge
	Matrix gpar[2];
	for (int d = 0; d < 2; d++)
		if (-1-d == edgeDir)
		{
			gpar[d].resize(1,1);
			gpar[d].fill(d == 0 ? lrspline->startparam_u() : lrspline->startparam_v());
		}
		else if (1+d == edgeDir)
		{
			gpar[d].resize(1,1);
			gpar[d].fill(d == 0 ? lrspline->endparam_u() : lrspline->endparam_v());
		}
		else
			this->getGaussPointParameters(gpar[d],d,nGauss,iel,xg);

	// Evaluate basis function derivatives at all integration points
	std::vector<Go::BasisDerivsSf> spline;
	lrspline->computeBasisGrid(gpar[0],gpar[1],spline);

	const int p1 = lrspline->order_u();
	const int p2 = lrspline->order_v();
	const int n1 = lrspline->numCoefs_u();
	const int n2 = lrspline->numCoefs_v();

	FiniteElement fe(p1*p2);
	fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
	fe.u = gpar[0](1,1);
	fe.v = gpar[1](1,1);

	Matrix dNdu, Xnod, Jac;
	Vec4   X;
	Vec3   normal;


	// === Assembly loop over all elements on the patch edge =====================

	std::vector<LR::Element*>::iterator el
	int iel = 1;
	for(el = lrspline->elementBegin(); el!=lrspline->elementEnd(); el++, iel++)
	{
		fe.iel = MLGE[iel-1];
		if (fe.iel < 1) continue; // zero-area element

		// Skip elements that are not on current boundary edge
		bool skipMe = false;
		switch (edgeDir)
		{
			case -1: if (i1 > p1) skipMe = true; break;
			case  1: if (i1 < n1) skipMe = true; break;
			case -2: if (i2 > p2) skipMe = true; break;
			case  2: if (i2 < n2) skipMe = true; break;
		}
			if (skipMe) continue;

			// Get element edge length in the parameter space
			double dS = this->getParametricLength(iel,t1);
			if (dS < 0.0) return false; // topology error (probably logic error)

			// Set up control point coordinates for current element
			if (!this->getElementCoordinates(Xnod,iel)) return false;

			// Initialize element quantities
			if (!integrand.initElementBou(MNPC[iel-1])) return false;

			// Caution: Unless locInt is empty, we assume it points to an array of
			// LocalIntegral pointers, of length at least the number of elements in
			// the model (as defined by the highest number in the MLGE array).
			// If the array is shorter than this, expect a segmentation fault.
			LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


			// --- Integration loop over all Gauss points along the edge -------------

			int ip = (t1 == 1 ? i2-p2 : i1-p1)*nGauss;
			for (int i = 0; i < nGauss; i++, ip++)
			{
	// Local element coordinates and parameter values
	// of current integration point
	if (gpar[0].size() > 1)
	{
		fe.xi = xg[i];
		fe.u = gpar[0](i+1,i1-p1+1);
	}
	if (gpar[1].size() > 1)
	{
		fe.eta = xg[i];
		fe.v = gpar[1](i+1,i2-p2+1);
	}

	// Fetch basis function derivatives at current integration point
	extractBasis(spline[ip],fe.N,dNdu);

	// Compute basis function derivatives and the edge normal
	fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	if (fe.detJxW == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.5*dS*wg[i];
	if (!integrand.evalBou(elmInt,fe,time,X,normal))
		return false;
			}

			// Assembly of global system integral
			if (!glInt.assemble(elmInt,fe.iel))
	return false;
		}

#endif
	return true;
}


int ASMu2D::evalPoint (const double* xi, double* param, Vec3& X) const
{
	if (!lrspline) return -2;

	param[0] = (1.0-xi[0])*lrspline->startparam_u() + xi[0]*lrspline->endparam_u();
	param[1] = (1.0-xi[1])*lrspline->startparam_v() + xi[1]*lrspline->endparam_v();

	Go::Point X0;
	lrspline->point(X0,param[0],param[1]);
	for (unsigned char d = 0; d < nsd; d++)
		X[d] = X0[d];

	return 0;
}


bool ASMu2D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
#ifdef SP_DEBUG
	std::cout << "ASMu2D::getGridParameters(  )\n";
#endif
	if(nSegPerSpan != 1) {
		std::cerr << "ASMu2D::getGridParameters called with nSegPerSpan != 2\n";
		return false;
	}

	// output is written once for each element resulting in a lot of unnecessary storage
	// this is preferable to figuring out all element topology information

	std::vector<LR::Element*>::iterator el;
	for(el=lrspline->elementBegin(); el<lrspline->elementEnd(); el++) {
		// evaluate element at element corner points
		double umin = (**el).umin();
		double umax = (**el).umax();
		double vmin = (**el).vmin();
		double vmax = (**el).vmax();
		for(int iv=0; iv<2; iv++) {
			for(int iu=0; iu<2; iu++) {
				double u = (iu==0) ? umin : umax;
				double v = (iv==0) ? vmin : vmax;
				if(dir==0)
					prm.push_back(u);
				else
					prm.push_back(v);
			}
		}
	}
	return true;
}

#if 0
	if (!lrspline) return false;

	if (nSegPerSpan < 1)
	{
		std::cerr <<" *** ASMu2D::getGridParameters: Too few knot-span points "
		          << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
		return false;
	}

	RealArray::const_iterator uit = lrspline->basis(dir).begin();
	double ucurr = 0.0, uprev = *(uit++);
	while (uit != lrspline->basis(dir).end())
	{
		ucurr = *(uit++);
		if (ucurr > uprev)
			if (nSegPerSpan == 1)
				prm.push_back(uprev);
			else for (int i = 0; i < nSegPerSpan; i++)
			{
				double xg = (double)(2*i-nSegPerSpan)/(double)nSegPerSpan;
				prm.push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
			}
		uprev = ucurr;
	}

	if (ucurr > prm.back())
		prm.push_back(ucurr);
	return true;
}


bool ASMu2D::getGrevilleParameters (RealArray& prm, int dir) const
{
	if (!lrspline) return false;

	const Go::BsplineBasis& basis = lrspline->basis(dir);

	prm.resize(basis.numCoefs());
	for (size_t i = 0; i < prm.size(); i++)
		prm[i] = basis.grevilleParameter(i);

	return true;
}
#endif


bool ASMu2D::tesselate (ElementBlock& grid, const int* npe) const
{
#ifdef SP_DEBUG
	std::cout << "ASMu2D::tesselate(  )\n";
#endif
	if(!lrspline) return false;

	if(npe[0] != 2 || npe[1] != 2) {
		std::cerr <<" *** ASMu2D::tesselate: unstructured tesselation only allows for 2 tesselation points\n";
		return false;
	}

	// output is written once for each element resulting in a lot of unnecessary storage
	// this is preferable to figuring out all element topology information
	grid.unStructResize((size_t) lrspline->nElements(), (size_t) lrspline->nElements()*4);

	std::vector<LR::Element*>::iterator el;
	int inod = 0;
	int iel = 0;
	for(el=lrspline->elementBegin(); el<lrspline->elementEnd(); el++, iel++) {
		// evaluate element at element corner points
		double umin = (**el).umin();
		double umax = (**el).umax();
		double vmin = (**el).vmin();
		double vmax = (**el).vmax();
		for(int iv=0; iv<2; iv++) {
			for(int iu=0; iu<2; iu++) {
				double u = (iu==0) ? umin : umax;
				double v = (iv==0) ? vmin : vmax;
				Go::Point pt;
				lrspline->point(pt, u,v, iel);
				for(int dim=0; dim<nsd; dim++)
					grid.setCoor(inod, dim, pt[dim]);
				inod++;
			}
		}
	}

	// due to duplicate point storage, the element-to-node array is sort of trivial
	for(int i=0; i<lrspline->nElements()*4; i+=4) {
		// enumerate nodes counter clockwise around the quad
		grid.setNode(i  , i  );
		grid.setNode(i+1, i+1);
		grid.setNode(i+2, i+3);
		grid.setNode(i+3, i+2);
	}

	return true;
}

#if 0


void ASMu2D::scatterInd (int n1, int n2, int p1, int p2,
                         const int* start, IntVec& index)
{
	index.reserve(p1*p2);
	int ip = ((start[1]-p2+1))*n1 + (start[0]-p1+1);
	for (int i2 = 0; i2 < p2; i2++, ip += n1-p1)
		for (int i1 = 0; i1 < p1; i1++, ip++)
			index.push_back(ip);
}

#endif


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe) const
{
#ifdef SP_DEBUG
	std::cout << "ASMu2D::evalSolution(Matrix, Vector, int* )\n";
#endif
	// Compute parameter values of the result sampling points
	RealArray gpar[2];
	for (int dir = 0; dir < 2; dir++)
		if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
			return false;

	// Evaluate the primary solution at all sampling points
	return this->evalSolution(sField,locSol,gpar);
}


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool regular) const
{
#ifdef SP_DEBUG
	std::cout << "ASMu2D::evalSolution(Matrix, Vector, RealArray )\n";
#endif
	size_t nComp = locSol.size() / this->getNoNodes();
	if (nComp*this->getNoNodes() != locSol.size())
		return false;

	if(gpar[0].size() != gpar[1].size())
		return false;

	Matrix Xtmp;
	Vector Ytmp;

	// Evaluate the primary solution field at each point
	size_t nPoints = gpar[0].size();
	sField.resize(nComp,nPoints);
	for (size_t i = 0; i < nPoints; i++)
	{
		// fetch element containing evaluation point
		int iel = i/4; // points are always listed in the same order as the elemnts, 4 pts per element

		// fetch index of non-zero basis functions on this element
		IntVec ip = MNPC[iel];

		// evaluate the basis functions at current parametric point
		Go::BasisPtsSf spline;
		lrspline->computeBasis(gpar[0][i], gpar[1][i], spline, iel);

		// Now evaluate the solution field
		utl::gather(ip,nComp,locSol,Xtmp);
		Xtmp.multiply(spline.basisValues, Ytmp);
		sField.fillColumn(1+i,Ytmp);
		// std::cout << "u(" << gpar[0][i] << ", " << gpar[1][i] << ") = " << Ytmp << "\n";
		// std::cout << "Xtmp = " << Xtmp << std::endl;
		// std::cout << "ip   = [" ;
		// for(uint i=0; i<ip.size(); i++) std::cout << ip[i] << ", ";
		// std::cout << "]\n";
		// std::vector<LR::Basisfunction*>::iterator it;
		// for(it=lrspline->getElement(iel)->supportBegin(); it<lrspline->getElement(iel)->supportEnd(); it++)
			// std::cout << **it << std::endl;
	}

	return true;
}

bool ASMu2D::evalSolution (Matrix& sField, const Integrand& integrand,
                           const int* npe, bool project) const
{
	if(project) {
		std::cerr << "ASMu2D::evalSolution unstructured LR-spline projection\n";
		std::cerr << "techniques not implemented yet\n";
		return false;
	}
	if (npe)
	{
		// Compute parameter values of the result sampling points
		RealArray gpar[2];
		if (this->getGridParameters(gpar[0],0,npe[0]-1) &&
		    this->getGridParameters(gpar[1],1,npe[1]-1))
			// Evaluate the secondary solution directly at all sampling points
			return this->evalSolution(sField,integrand,gpar);
	}
	else
	{
		std::cerr << "ASMu2D::evalSolution unstructured LR-spline projection\n";
		std::cerr << "techniques not implemented yet\n";
		return false;
	}

	std::cerr <<" *** ASMu2D::evalSolution: Failure!"<< std::endl;
	return false;
}

#if 0

Go::GeomObject* ASMu2D::evalSolution (const Integrand& integrand) const
{
	return this->projectSolution(integrand);
}


Go::SplineSurface* ASMu2D::projectSolution (const Integrand& integrand) const
{
	// Compute parameter values of the result sampling points (Greville points)
	RealArray gpar[2];
	for (int dir = 0; dir < 2; dir++)
		if (!this->getGrevilleParameters(gpar[dir],dir))
			return 0;

	// Evaluate the secondary solution at all sampling points
	Matrix sValues;
	if (!this->evalSolution(sValues,integrand,gpar))
		return 0;

	// Project the results onto the spline basis to find control point
	// values based on the result values evaluated at the Greville points.
	// Note that we here implicitly assume the the number of Greville points
	// equals the number of control points such that we don't have to resize
	// the result array. Think that is always the case, but beware if trying
	// other projection schemes later.

	RealArray weights;
	if (lrspline->rational())
		lrspline->getWeights(weights);

	const Vector& vec = sValues;
	return Go::SurfaceInterpolator::regularInterpolation(lrspline->basis(0),
	                                                     lrspline->basis(1),
	                                                     gpar[0], gpar[1],
	                                                     const_cast<Vector&>(vec),
	                                                     sValues.rows(),
	                                                     lrspline->rational(),
	                                                     weights);
}

#endif

bool ASMu2D::evalSolution (Matrix& sField, const Integrand& integrand,
                           const RealArray* gpar, bool regular) const
{
#ifdef SP_DEBUG
	std::cout << "ASMu2D::evalSolution(  )\n";
#endif
	// TODO: investigate the possibility of doing "regular" refinement by unfiorm
	//       tesselation grid and ignoring LR meshlines

	if(gpar[0].size() != gpar[1].size())
		return false;

	sField.resize(0,0);

	// Fetch nodal (control point) coordinates
	Matrix Xnod;

	IntVec ip;
	Vector solPt;
	Matrix dNdu, dNdX, Jac;

	// Evaluate the secondary solution field at each point
	size_t nPoints = gpar[0].size();

	for (size_t i = 0; i < nPoints; i++)
	{
		// fetch element containing evaluation point
		int iel = i/4; // points are always listed in the same order as the elemnts, 4 pts per element

		// fetch number of nonzero basis functions on this element as well as their indices
		int nBasis = lrspline->getElement(iel)->nBasisFunctions();
		Vector N(nBasis);
		ip = MNPC[iel];

		// evaluate the basis functions at current parametric point
		Go::BasisDerivsSf spline;
		lrspline->computeBasis(gpar[0][i], gpar[1][i], spline, iel);
		extractBasis(spline, N, dNdu);

		// Set up control point (nodal) coordinates for current element
		if (!this->getElementCoordinates(Xnod,iel+1)) return false;

		// Compute the Jacobian inverse
		if (utl::Jacobian(Jac,dNdX,Xnod,dNdu) == 0.0) // Jac = (Xnod * dNdu)^-1
			continue; // skip singular points

		// Now evaluate the solution field
		if (!integrand.evalSol(solPt,N,dNdX,Xnod*N,ip))
			return false;
		else if (sField.empty())
			sField.resize(solPt.size(),nPoints,true);

		sField.fillColumn(1+i,solPt);
	}

	return true;
}

