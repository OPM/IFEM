// $Id$
//==============================================================================
//!
//! \file ASMu3D.C
//!
//! \date January 2013
//!
//! \author Kjetil A. Johannessen / NTNU
//!
//! \brief Driver for assembly of unstructured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu3D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "MPC.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


ASMu3D::ASMu3D (unsigned char n_f)
	: ASMunstruct(3,3,n_f), lrspline(0), tensorspline(0), workingEl(-1)
{
	ASMunstruct::resetNumbering(); // Replace this when going multi-patch...
	swapW = false;
}


ASMu3D::ASMu3D (const ASMu3D& patch, unsigned char n_f)
	: ASMunstruct(patch,n_f), lrspline(patch.lrspline), tensorspline(0), workingEl(-1)
{
	swapW = patch.swapW;
}

size_t ASMu3D::getNodeIndex (int globalNum, bool noAddedNodes) const
{
	std::cerr << "ASM3D::getNodeIndex not implemented properly\n";
	exit(123213);

#if 0
	IntVec::const_iterator it = std::find(MLGN.begin(),MLGN.end(),globalNum);
	if (it == MLGN.end()) return 0;
	
	size_t inod = 1 + (it-MLGN.begin());
	if (noAddedNodes && !xnMap.empty() && inod > nodeInd.size())
	{
		std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
		if (it != xnMap.end()) return it->second;
	}
	
	return inod;
#endif
	return 0;
}

/*
char ASMu3D::getNodeType (size_t inod) const
{
	std::cerr << "ASM3D::getNodeType not implemented properly\n";
	exit(123213);
	if (this->isLMn(inod)) return 'L';
	// return inod > nodeInd.size() ? 'X' : 'D';
	return 0;
}
*/


LR::LRSplineSurface* ASMu3D::getBoundary (int dir)
{
	std::cerr << "ASMu3D::getBoundary not implemented properly yet" << std::endl;
	exit(776654);
	if (dir < -3 || dir == 0 || dir > 3)
		return NULL;

	// The boundary surfaces are stored internally in the SplineVolume object
	// int iface = dir > 0 ? 2*dir-1 : -2*dir-2;
	return NULL;
	// return lrspline->getBoundarySurface(iface).get();
}


bool ASMu3D::read (std::istream& is)
{
	if (shareFE) return true;
	if (lrspline) delete lrspline;

	// read inputfile as either an LRSpline file directly or a tensor product B-spline and convert
	char firstline[256];
	is.getline(firstline, 256);
	if(strncmp(firstline, "# LRSPLINE", 10) == 0) {
		lrspline = new LR::LRSplineVolume();
		is >> *lrspline;
	} else { // probably a SplineVolume, so we'll read that and convert
		tensorspline = new Go::SplineVolume();
		is >> *tensorspline;
		lrspline = new LR::LRSplineVolume(tensorspline);
	}

	// Eat white-space characters to see if there is more data to read
	char c;
	while (is.get(c))
		if (!isspace(c)) {
			is.putback(c);
			break;
	}

	if (!is.good() && !is.eof())
	{
		std::cerr <<" *** ASMu3D::read: Failure reading spline data"<< std::endl;
		delete lrspline;
		lrspline = 0;
		return false;
	}
	else if (lrspline->dimension() < 3)
	{
		std::cerr <<" *** ASMu3D::read: Invalid spline volume patch, dim="
	      << lrspline->dimension() << std::endl;
		delete lrspline;
		lrspline = 0;
		return false;
	}

	geo = lrspline;
	return true;
}


bool ASMu3D::write (std::ostream& os, int) const
{
	if (!lrspline) return false;

	os << *lrspline;

	return os.good();
}


void ASMu3D::clear (bool retainGeometry)
{
	if (!retainGeometry)
	{
		// Erase spline data
		if (lrspline && !shareFE) delete lrspline;
		lrspline = 0;
		geo = 0;
		tensorspline = 0;
		tensorspline = 0;
	}

	// Erase the FE data
	this->ASMbase::clear(retainGeometry);
	xnMap.clear();
	nxMap.clear();
}

size_t ASMu3D::getNoNodes (int basis) const
{
	return lrspline->nBasisFunctions();
}

bool ASMu3D::checkRightHandSystem ()
{
	if (!lrspline || shareFE) return false;
	std::cerr << "ASMu3D::checkRightHandSystem() not properly implemented yet" << std::endl;
	exit(776654);

	// Evaluate the spline volume at its center
	double u = 0.5*(lrspline->startparam(0) + lrspline->endparam(0));
	double v = 0.5*(lrspline->startparam(1) + lrspline->endparam(1));
	double w = 0.5*(lrspline->startparam(2) + lrspline->endparam(2));
	std::vector<Go::Point> pts;
	lrspline->point(pts,u,v,w,1);

	// Check that |J| = (dXdu x dXdv) * dXdw > 0.0
	if ((pts[1] % pts[2]) * pts[3] > 0.0) return false;

	// This patch has a negative Jacobian determinant. Probably it is modelled
	// in a left-hand-system. Swap the w-parameter direction to correct for this.

	// lrspline->reverseParameterDirection(2);

	std::cout <<"\tSwapped."<< std::endl;
	return swapW = true;
}


bool ASMu3D::refine (int dir, const RealArray& xi)
{
	if (!tensorspline || dir < 0 || dir > 2 || xi.empty()) return false;
	if (xi.front() < 0.0 || xi.back() > 1.0) return false;
	if (shareFE) return true;

	RealArray extraKnots;
	RealArray::const_iterator uit = tensorspline->basis(dir).begin();
	double ucurr, uprev = *(uit++);
	while (uit != tensorspline->basis(dir).end())
	{
		ucurr = *(uit++);
		if (ucurr > uprev)
			for (size_t i = 0; i < xi.size(); i++)
	if (i > 0 && xi[i] < xi[i-1])
	  return false;
	else
	  extraKnots.push_back(ucurr*xi[i] + uprev*(1.0-xi[i]));

		uprev = ucurr;
	}

	tensorspline->insertKnot(dir,extraKnots);
	if(lrspline) delete lrspline;
	geo = lrspline = new LR::LRSplineVolume(tensorspline);
	return true;
}

bool ASMu3D::uniformRefine (int dir, int nInsert)
{
	if (!tensorspline || dir < 0 || dir > 2 || nInsert < 1) return false;
	if (shareFE) return true;

	RealArray extraKnots;
	RealArray::const_iterator uit = tensorspline->basis(dir).begin();
	double ucurr, uprev = *(uit++);
	while (uit != tensorspline->basis(dir).end())
	{
		ucurr = *(uit++);
		if (ucurr > uprev)
			for (int i = 0; i < nInsert; i++)
			{
				double xi = (double)(i+1)/(double)(nInsert+1);
				extraKnots.push_back(ucurr*xi + uprev*(1.0-xi));
			}
		uprev = ucurr;
	}

	tensorspline->insertKnot(dir,extraKnots);
	if(lrspline) delete lrspline;
	geo = lrspline = new LR::LRSplineVolume(tensorspline);
	return true;
}

bool ASMu3D::raiseOrder (int ru, int rv, int rw)
{
	if (!tensorspline) return false;
	if (shareFE) return true;

	tensorspline->raiseOrder(ru,rv,rw);
	delete lrspline;
	geo = lrspline = new LR::LRSplineVolume(tensorspline);
	return true;
}

bool ASMu3D::generateFEMTopology ()
{
	// At this point we are through with the tensor spline object.
	// So release it to avoid memory leakage.
	delete tensorspline;
	tensorspline = 0;

	if (!lrspline) return false;

	const int nBasis    = lrspline->nBasisFunctions();
	const int nElements = lrspline->nElements();

	if ((size_t)nBasis == MLGN.size())
		return true;
	else if (!MLGN.empty())
	{
		std::cerr <<" *** ASMu3D::generateFEMTopology: Inconsistency"
		          <<" between the number of FE nodes "<< MLGN.size()
		          <<"\n     and the number of basis functions "<< nBasis
		          <<" in the patch."<< std::endl;
		return false;
	}
	else if (shareFE)
		return true;

	const int p1 = lrspline->order(0);
	const int p2 = lrspline->order(1);
	const int p3 = lrspline->order(2);

	// Consistency checks, just to be fool-proof
	if (nBasis < 8)                 return false;
	if (p1 < 1 || p2 < 1 || p3 < 1) return false;

	myMLGN.resize(nBasis);
	myMLGE.resize(nElements);
	myMNPC.resize(nElements);
	lrspline->generateIDs();

	std::vector<LR::Element*>::iterator el_it = lrspline->elementBegin();
	for (int iel=0; iel<nElements; iel++, el_it++)
	{
		LR::Element *el = *el_it;
		int nSupportFunctions = el->nBasisFunctions();
		myMLGE[iel] = ++gEl; // global element number over all patches
		myMNPC[iel].resize(nSupportFunctions);

		int lnod = 0;
		for(LR::Basisfunction *b : el->support())
			myMNPC[iel][lnod++] = b->getId();
	}

	for (int inod = 0; inod < nBasis; inod++)
		myMLGN[inod] = ++gNod;

	return true;
}



bool ASMu3D::connectPatch (int face, ASMu3D& neighbor, int nface, int norient)
{
	std::cerr << "ASMu3D::connectPatch(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	if (swapW && face > 4) // Account for swapped parameter direction
		face = 11-face;

	if (neighbor.swapW && face > 4) // Account for swapped parameter direction
		nface = 11-nface;

	return this->connectBasis(face,neighbor,nface,norient);
#endif
	return false;
}


bool ASMu3D::connectBasis (int face, ASMu3D& neighbor, int nface, int norient,
				                   int basis, int slave, int master)
{
	std::cerr << "ASMu3D::connectBasis(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	if (shareFE && neighbor.shareFE)
		return true;
	else if (shareFE || neighbor.shareFE)
	{
		std::cerr <<" *** ASMu3D::connectPatch: Logic error, cannot"
	      <<" connect a sharedFE patch with an unshared one"<< std::endl;
		return false;
	}

	// Set up the slave node numbers for this volume patch

	int n1, n2, n3;
	if (!this->getSize(n1,n2,n3,basis)) return false;
	int node = slave+1, i1 = 0, i2 = 0;

	switch (face)
		{
		case 2: // Positive I-direction
			node += n1-1;
		case 1: // Negative I-direction
			i1 = n1;
			n1 = n2;
			n2 = n3;
			break;

		case 4: // Positive J-direction
			node += n1*(n2-1);
		case 3: // Negative J-direction
			i2 = n1*(n2-1);
			i1 = 1;
			n2 = n3;
			break;

		case 6: // Positive K-direction
			node += n1*n2*(n3-1);
		case 5: // Negative K-direction
			i1 = 1;
			break;

		default:
			std::cerr <<" *** ASMu3D::connectPatch: Invalid slave face "
		<< face << std::endl;
			return false;
		}

	int i, j;
	IntMat slaveNodes(n1,IntVec(n2,0));
	for (j = 0; j < n2; j++, node += i2)
		for (i = 0; i < n1; i++, node += i1)
			slaveNodes[i][j] = node;

	// Set up the master node numbers for the neighboring volume patch

	if (!neighbor.getSize(n1,n2,n3,basis)) return false;
	node = master+1; i1 = i2 = 0;

	switch (nface)
		{
		case 2: // Positive I-direction
			node += n1-1;
		case 1: // Negative I-direction
			i1 = n1;
			n1 = n2;
			n2 = n3;
			break;

		case 4: // Positive J-direction
			node += n1*(n2-1);
		case 3: // Negative J-direction
			i2 = n1*(n2-1);
			i1 = 1;
			n2 = n3;
			break;

		case 6: // Positive K-direction
			node += n1*n2*(n3-1);
		case 5: // Negative K-direction
			i1 = 1;
			break;

		default:
			std::cerr <<" *** ASMu3D::connectPatch: Invalid master face "
		<< nface << std::endl;
			return false;
		}

	if (norient < 0 || norient > 7)
	{
		std::cerr <<" *** ASMu3D::connectPatch: Orientation flag "
	      << norient <<" is out of range [0,7]"<< std::endl;
		return false;
	}

	int m1 = slaveNodes.size();
	int m2 = slaveNodes.front().size();
	if (norient < 4 ? (n1 != m1 || n2 != m2) : (n2 != m1 || n1 != m2))
	{
		std::cerr <<" *** ASMu3D::connectPatch: Non-matching faces, sizes "
	      << n1 <<","<< n2 <<" and "<< m1 <<","<< m2 << std::endl;
		return false;
	}

	const double xtol = 1.0e-4;
	for (j = 0; j < n2; j++, node += i2)
		for (i = 0; i < n1; i++, node += i1)
		{
			int k = i, l = j;
			switch (norient)
	{
	case  1: k =    i  ; l = n2-j-1; break;
	case  2: k = n1-i-1; l =    j  ; break;
	case  3: k = n1-i-1; l = n2-j-1; break;
	case  4: k =    j  ; l =    i  ; break;
	case  5: k =    j  ; l = n1-i-1; break;
	case  6: k = n2-j-1; l =    i  ; break;
	case  7: k = n2-j-1; l = n1-i-1; break;
	default: k =    i  ; l = j     ;
	}

			int slave = slaveNodes[k][l];
			if (!neighbor.getCoord(node).equal(this->getCoord(slave),xtol))
			{
	std::cerr <<" *** ASMu3D::connectPatch: Non-matching nodes "
		  << node <<": "<< neighbor.getCoord(node)
		  <<"\n                                          and "
		  << slave <<": "<< this->getCoord(slave) << std::endl;
	return false;
			}
			else
	ASMbase::collapseNodes(neighbor,node,*this,slave);
		}

	return true;
#endif
	return false;
}


void ASMu3D::closeFaces (int dir, int basis, int master)
{
	std::cerr << "ASMu3D::closeFaces(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	int n1, n2, n3;
	if (basis < 1) basis = 1;
	if (!this->getSize(n1,n2,n3,basis)) return;

	switch (dir)
		{
		case 1: // Faces are closed in I-direction
			for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, master += n1)
	  this->makePeriodic(master,master+n1-1);
			break;

		case 2: // Faces are closed in J-direction
			for (int i3 = 1; i3 <= n3; i3++, master += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, master++)
	  this->makePeriodic(master,master+n1*(n2-1));
			break;

		case 3: // Faces are closed in K-direction
			for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, master++)
	  this->makePeriodic(master,master+n1*n2*(n3-1));
			break;
		}
	#endif
}

/*!
	A negative \a code value implies direct evaluation of the Dirichlet condition
	function at the control point. Positive \a code implies projection onto the
	spline basis representing the boundary surface (needed for curved faces and/or
	non-constant functions).
*/

void ASMu3D::constrainFace (int dir, bool open, int dof, int code)
{
	if (swapW) // Account for swapped parameter direction
		if (dir == 3 || dir == -3) dir = -dir;
	
	if(open)
		std::cerr << "\nWARNING: ASMu3D::constrainFace, open boundary conditions not supported yet. Treating it as closed" << std::endl;

	int bcode = code;
	if (code > 0) {// Dirichlet projection will be performed
		std::cerr << "\nWARNING: Projective (nonhomogenuous) dirichlet boundary conditions not implemented. ";
		std::cerr << "\n         Performing variational diminishing approximation instead" << std::endl;
		// dirich.push_back(DirichletFace(this->getBoundary(dir),dof,code));
	}
	else if (code < 0)
		bcode = -code;

	// get all the boundary functions from the LRspline object
	std::vector<LR::Basisfunction*> thisEdge;
	if(dir == -1) 
		lrspline->getEdgeFunctions(thisEdge, LR::WEST, 1);
	else if(dir == 1) 
		lrspline->getEdgeFunctions(thisEdge, LR::EAST, 1);
	else if(dir == -2) 
		lrspline->getEdgeFunctions(thisEdge, LR::SOUTH, 1);
	else if(dir == 2) 
		lrspline->getEdgeFunctions(thisEdge, LR::NORTH, 1);
	else if(dir == -3) 
		lrspline->getEdgeFunctions(thisEdge, LR::BOTTOM, 1);
	else if(dir == 3) 
		lrspline->getEdgeFunctions(thisEdge, LR::TOP, 1);

	std::cout << "\nNumber of constraints: " << thisEdge.size() << std::endl;

	// enforce the boundary conditions
	for(LR::Basisfunction* b  : thisEdge)
		this->prescribe(b->getId()+1,dof,bcode);
}

size_t constrainFaceLocal(int dir, bool open, int dof, int code, bool project, char T1)
{
	std::cerr << "ASMu3D::constrainFaceLocal not implemented properly yet" << std::endl;
	exit(776654);
	return 0;
}

void ASMu3D::constrainEdge (int lEdge, bool open, int dof, int code)
{
	if(open)
		std::cerr << "\nWARNING: ASMu3D::constrainEdge, open boundary conditions not supported yet. Treating it as closed" << std::endl;

	if (swapW && lEdge <= 8) // Account for swapped parameter direction
		lEdge += (lEdge-1)%4 < 2 ? 2 : -2;

	// lEdge = 1-4, running index is u (vmin,wmin), (vmax,wmin), (vmin,wmax), (vmax,wmax)
	// lEdge = 5-8, running index is v (umin,wmin), (umax,wmin), (umin,wmax), (umax,wmax)
	// lEdge = 9-12, running index is w 

	int edge = LR::NONE;
	if(lEdge == 1) 
		edge = LR::BOTTOM | LR::SOUTH;
	else if(lEdge == 2) 
		edge = LR::BOTTOM | LR::NORTH;
	else if(lEdge == 3) 
		edge = LR::TOP    | LR::SOUTH;
	else if(lEdge == 4) 
		edge = LR::TOP    | LR::NORTH;
	else if(lEdge == 5) 
		edge = LR::BOTTOM | LR::WEST;
	else if(lEdge == 6) 
		edge = LR::BOTTOM | LR::EAST;
	else if(lEdge == 7) 
		edge = LR::TOP    | LR::WEST;
	else if(lEdge == 8) 
		edge = LR::TOP    | LR::EAST;
	else if(lEdge == 9) 
		edge = LR::SOUTH  | LR::WEST;
	else if(lEdge == 10) 
		edge = LR::SOUTH  | LR::EAST;
	else if(lEdge == 11) 
		edge = LR::NORTH  | LR::WEST;
	else if(lEdge == 12) 
		edge = LR::NORTH  | LR::EAST;

	// get all the boundary functions from the LRspline object
	std::vector<LR::Basisfunction*> thisEdge;
	lrspline->getEdgeFunctions(thisEdge, (LR::parameterEdge) edge, 1);

	// enforce the boundary conditions
	for(LR::Basisfunction* b  : thisEdge)
		this->prescribe(b->getId(),dof,code);
}


void ASMu3D::constrainLine (int fdir, int ldir, double xi, int dof, int code)
{
	std::cerr << "ASMu3D::constrainLine not implemented properly yet" << std::endl;
	exit(776654);
	#if 0
	if (xi < 0.0 || xi > 1.0) return;

	int n1, n2, n3, node = 1;
	if (!this->getSize(n1,n2,n3,1)) return;

	if (swapW) // Account for swapped parameter direction
		if (fdir == 3 || fdir == -3) fdir = -fdir;

	switch (fdir)
		{
		case  1: // Right face (positive I-direction)
			node += n1-1;
		case -1: // Left face (negative I-direction)
			if (ldir == 2)
			{
	// Line goes in J-direction
	node += n1*n2*int(0.5+(n3-1)*(swapW ? 1.0-xi : xi));
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  this->prescribe(node,dof,code);
			}
			else if (ldir == 3)
			{
	// Line goes in K-direction
	node += n1*int(0.5+(n2-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  this->prescribe(node,dof,code);
			}
			break;

		case  2: // Back face (positive J-direction)
			node += n1*(n2-1);
		case -2: // Front face (negative J-direction)
			if (ldir == 1)
			{
	// Line goes in I-direction
	node += n1*n2*int(0.5+(n3-1)*(swapW ? 1.0-xi : xi));
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  this->prescribe(node,dof,code);
			}
			else if (ldir == 3)
			{
	// Line goes in K-direction
	node += int(0.5+(n1-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  this->prescribe(node,dof,code);
			}
			break;

		case  3: // Top face (positive K-direction)
			node += n1*n2*(n3-1);
		case -3: // Bottom face (negative K-direction)
			if (ldir == 1)
			{
	// Line goes in I-direction
	node += n1*int(0.5+(n2-1)*xi);
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  this->prescribe(node,dof,code);
			}
			else if (ldir == 2)
			{
	// Line goes in J-direction
	node += int(0.5+(n1-1)*xi);
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  this->prescribe(node,dof,code);
			}
			break;
		}
	#endif
}


void ASMu3D::constrainCorner (int I, int J, int K, int dof, int code)
{
	std::cerr << "ASMu3D::constrainCorner not implemented properly yet" << std::endl;
	exit(776654);
	#if 0
	int n1, n2, n3;
	if (!this->getSize(n1,n2,n3,1)) return;

	if (swapW) // Account for swapped parameter direction
		K = -K;

	int node = 1;
	if (I > 0) node += n1-1;
	if (J > 0) node += n1*(n2-1);
	if (K > 0) node += n1*n2*(n3-1);

	this->prescribe(node,dof,code);
	#endif
}


void ASMu3D::constrainNode (double xi, double eta, double zeta,
			    int dof, int code)
{
	std::cerr << "ASMu3D::constrainNode not implemented properly yet" << std::endl;
	exit(776654);
	if (xi   < 0.0 || xi   > 1.0) return;
	if (eta  < 0.0 || eta  > 1.0) return;
	if (zeta < 0.0 || zeta > 1.0) return;

	if (swapW) // Account for swapped parameter direction
		zeta = 1.0-zeta;

#if 0
	int n1, n2, n3;
	if (!this->getSize(n1,n2,n3,1)) return;

	int node = 1;
	if (xi   > 0.0) node += int(0.5+(n1-1)*xi);
	if (eta  > 0.0) node += n1*int(0.5+(n2-1)*eta);
	if (zeta > 0.0) node += n1*n2*int(0.5+(n3-1)*zeta);

	this->prescribe(node,dof,code);
#endif
}


/*!
	This method projects the function describing the in-homogeneous Dirichlet
	boundary condition onto the spline basis defining the boundary surface,
	in order to find the control point values which are used as the prescribed
	values of the boundary DOFs.
*/

bool ASMu3D::updateDirichlet (const std::map<int,RealFunc*>& func,
			      const std::map<int,VecFunc*>& vfunc, double time)
{
	// std::cerr << "ASMu3D::updateDirichlet not implemented properly yet\n";
	std::cerr << "\nWARNING: ASMu3D::updateDirichlet ignored due to non-projecting boundary condition implementation" << std::endl;
	return true;
	#if 0
	std::map<int,RealFunc*>::const_iterator fit;
	std::map<int,VecFunc*>::const_iterator vfit;
	std::vector<DirichletFace>::const_iterator dit;
	std::vector<Ipair>::const_iterator nit;

	for (size_t i = 0; i < dirich.size(); i++)
	{
		// Project the function onto the spline surface basis
		Go::SplineSurface* dsurf = 0;
		if ((fit = func.find(dirich[i].code)) != func.end())
			dsurf = SplineUtils::project(dirich[i].surf,*fit->second,time);
		else if ((vfit = vfunc.find(dirich[i].code)) != vfunc.end())
			dsurf = SplineUtils::project(dirich[i].surf,*vfit->second,nf,time);
		else
		{
			std::cerr <<" *** ASMu3D::updateDirichlet: Code "<< dirich[i].code
		<<" is not associated with any function."<< std::endl;
			return false;
		}
		if (!dsurf)
		{
			std::cerr <<" *** ASMu3D::updateDirichlet: Projection failure."
		<< std::endl;
			return false;
		}

		// Loop over the (interior) nodes (control points) of this boundary surface
		for (nit = dirich[i].nodes.begin(); nit != dirich[i].nodes.end(); nit++)
			for (int dofs = dirich[i].dof; dofs > 0; dofs /= 10)
			{
				int dof = dofs%10;
				// Find the constraint equation for current (node,dof)
				MPC pDOF(MLGN[nit->second-1],dof);
				MPCIter mit = mpcs.find(&pDOF);
				if (mit == mpcs.end()) continue; // probably a deleted constraint

				// Find index to the control point value for this (node,dof) in dsurf
				RealArray::const_iterator cit = dsurf->coefs_begin();
				if (dsurf->dimension() > 1) // A vector field is specified
				  cit += (nit->first-1)*dsurf->dimension() + (dof-1);
				else // A scalar field is specified at this dof
				  cit += (nit->first-1);

				// Now update the prescribed value in the constraint equation
				(*mit)->setSlaveCoeff(*cit);
#if SP_DEBUG > 1
				std::cout <<"Updated constraint: "<< **mit;
#endif
			}
		delete dsurf;
	}

	// The parent class method takes care of the corner nodes with direct
	// evaluation of the Dirichlet functions (since they are interpolatory)
	return this->ASMbase::updateDirichlet(func,vfunc,time);
	#endif
}

#if 0

double ASMu3D::getParametricArea (int iel, int dir) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu3D::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return DERR;
	}
#endif
	if (MNPC[iel-1].empty())
		return 0.0;

	int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
	if (inod1 < 0 || (size_t)inod1 >= nodeInd.size())
	{
		std::cerr <<" *** ASMu3D::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
		return DERR;
	}
#endif

	const int ni = nodeInd[inod1].I;
	const int nj = nodeInd[inod1].J;
	const int nk = nodeInd[inod1].K;
	switch (dir)
		{
		case 1: return lrspline->knotSpan(1,nj)*lrspline->knotSpan(2,nk);
		case 2: return lrspline->knotSpan(0,ni)*lrspline->knotSpan(2,nk);
		case 3: return lrspline->knotSpan(0,ni)*lrspline->knotSpan(1,nj);
		}

	std::cerr <<" *** ASMu3D::getParametricArea: Invalid face direction "
	    << dir << std::endl;
	return DERR;
}


int ASMu3D::coeffInd (size_t inod) const
{
#ifdef INDEX_CHECK
	if (inod >= nodeInd.size())
	{
		std::cerr <<" *** ASMu3D::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
		return -1;
	}
#endif

	const int ni = nodeInd[inod].I;
	const int nj = nodeInd[inod].J;
	const int nk = nodeInd[inod].K;
	return (nk*lrspline->numCoefs(1) + nj)*lrspline->numCoefs(0) + ni;
}
#endif

Vec3 ASMu3D::getCoord (size_t inod) const
{
	if (inod <= MLGN.size())
	{
		// This is a node added due to constraints in local directions.
		// Find the corresponding original node (see constrainEdgeLocal)
		std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
		if (it != xnMap.end()) inod = it->second;
	}
	
	if (inod == 0) return Vec3();
	LR::Basisfunction *b = lrspline->getBasisfunction(inod-1);

	RealArray::const_iterator cit = b->cp();
	return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMu3D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
	if (iel < 1 || (size_t)iel > MNPC.size())
	{
		std::cerr <<" *** ASMu3D::getElementCoordinates: Element index "<< iel
		          <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
		return false;
	}
#endif

	const IntVec& mnpc = MNPC[iel-1];
	X.resize(3,mnpc.size());

	int n = 1;
	for (LR::Basisfunction *b : lrspline->getElement(iel-1)->support() ) {
		for (size_t i = 1; i <= 3; i++)
			X(i,n) = b->cp(i-1);
		n++;
	}

#if SP_DEBUG > 2
	std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
	return true;
}


void ASMu3D::getNodalCoordinates (Matrix& X) const
{
	const int n = lrspline->nBasisFunctions();
	X.resize(3,n);

	size_t inod = 1;
	for(LR::Basisfunction *b : lrspline->getAllBasisfunctions() )
		for (size_t i = 0; i < 3; i++)
			X(i+1,inod++) = b->cp(i);
}


bool ASMu3D::updateCoords (const Vector& displ)
{
	std::cerr <<"ASMu3D::updateCoords not implemented properly yet" << std::endl;
	exit(776654);
	if (!lrspline) return true; // silently ignore empty patches
	if (shareFE) return true;

	if (displ.size() != 3*MLGN.size())
	{
		std::cerr <<" *** ASMu3D::updateCoords: Invalid dimension "
		          << displ.size() <<" on displacement vector, should be "
		          << 3*MLGN.size() << std::endl;
		return false;
	}

	// lrspline->deform(displ,3);
	return true;
}


bool ASMu3D::getOrder (int& p1, int& p2, int& p3) const
{
	if (!lrspline) return false;

	p1 = lrspline->order(0);
	p2 = lrspline->order(1);
	p3 = lrspline->order(2);
	return true;
}


#if  0
bool ASMu3D::getSize (int& n1, int& n2, int& n3, int) const
{
	if (!lrspline) return false;

	n1 = lrspline->nBasisFunctions();
	n2 = lrspline->nBasisFunctions();
	n3 = lrspline->nBasisFunctions();
	return true;
}


size_t ASMu3D::getNoBoundaryElms (char lIndex, char ldim) const
{
	if (!lrspline) return 0;

	if (ldim < 1 && lIndex > 0)
		return 1;
	else if (ldim < 2 && lIndex > 0 && lIndex <= 12)
		return lrspline->numCoefs((lIndex-1)/4) - lrspline->order((lIndex-1)/4) + 1;

	int n1 = lrspline->numCoefs(0) - lrspline->order(0) + 1;
	int n2 = lrspline->numCoefs(1) - lrspline->order(1) + 1;
	int n3 = lrspline->numCoefs(2) - lrspline->order(2) + 1;

	switch (lIndex)
		{
		case 1:
		case 2:
			return n2*n3;
		case 3:
		case 4:
			return n1*n3;
		case 5:
		case 6:
			return n1*n2;
		}

	return 0;
}
#endif


void ASMu3D::getGaussPointParameters (RealArray& uGP, int dir, int nGauss,
				                              int iEl, const double* xi) const
{
	LR::Element *el = lrspline->getElement(iEl);
	double start = el->getParmin(dir);
	double stop  = el->getParmax(dir);

	uGP.resize(nGauss);
	
	for(int i=0; i< nGauss; i++)	
		uGP[i] = 0.5*((stop-start)*xi[i] + stop+start);
}


void ASMu3D::getElementCorners (int iEl, std::vector<Vec3>& XC) const
{
	LR::Element *el = lrspline->getElement(iEl);
	double u[] = { el->getParmin(0), el->getParmax(0) };
	double v[] = { el->getParmin(1), el->getParmax(1) };
	double w[] = { el->getParmin(2), el->getParmax(2) };

	XC.resize(8);
	int ip=0;
	for(int k=0; k<2; k++) 
		for(int j=0; j<2; j++) 
			for(int i=0; i<2; i++) {
				Go::Point pt;
				lrspline->point(pt, u[i], v[j], w[k], iEl);
				XC[ip++] = Vec3(pt);
	}
}


void ASMu3D::evaluateBasis (const FiniteElement &el, Vector &N, Matrix &dNdu)
{
	PROFILE2("Spline evaluation");
	size_t nBasis = lrspline->getElement(el.iel)->nBasisFunctions();

	std::vector<std::vector<double> > result;
	lrspline->computeBasis(el.u, el.v, el.w, result, 1, el.iel );

	N.resize(nBasis);
	dNdu.resize(nBasis,3);
	size_t jp, n = 1;
	for (jp = 0; jp < nBasis; jp++, n++) {
		N     (n)     = result[jp][0];
		dNdu (n,1)    = result[jp][1];
		dNdu (n,2)    = result[jp][2];
		dNdu (n,3)    = result[jp][3];
	}
}

void ASMu3D::evaluateBasis (const FiniteElement &el, Vector &N, Matrix &dNdu, Matrix3D &d2Ndu2)
{
	PROFILE2("Spline evaluation");
	size_t nBasis = lrspline->getElement(el.iel)->nBasisFunctions();

	std::vector<std::vector<double> > result;
	lrspline->computeBasis(el.u, el.v, el.w, result, 2, el.iel );

	N.resize(nBasis);
	dNdu.resize(nBasis,3);
	d2Ndu2.resize(nBasis,3,3);
	size_t jp, n = 1;
	for (jp = 0; jp < nBasis; jp++, n++) {
		N   (n)       = result[jp][0];
		dNdu (n,1)    = result[jp][1];
		dNdu (n,2)    = result[jp][2];
		dNdu (n,3)    = result[jp][3];
		d2Ndu2(n,1,1) = result[jp][4];
		d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = result[jp][5];
		d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = result[jp][6];
		d2Ndu2(n,2,2) = result[jp][7];
		d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = result[jp][8];
		d2Ndu2(n,3,3) = result[jp][9];
	}
}

void ASMu3D::evaluateBasis (FiniteElement &el, int derivs)
{
	PROFILE2("Spline evaluation");
	size_t nBasis = lrspline->getElement(el.iel)->nBasisFunctions();

	std::vector<std::vector<double> > result;
	lrspline->computeBasis(el.u, el.v, el.w, result, derivs, el.iel );

	el.N.resize(nBasis);
	el.dNdX.resize(nBasis,3);
	el.d2NdX2.resize(nBasis,3,3);
	size_t jp, n = 1;
	for (jp = 0; jp < nBasis; jp++, n++) {
		el.N   (n)      = result[jp][0];
		if(derivs > 0) {  
			el.dNdX (n,1)    = result[jp][1];
			el.dNdX (n,2)    = result[jp][2];
			el.dNdX (n,3)    = result[jp][3];
		}
		if(derivs > 1) {
			el.d2NdX2(n,1,1) = result[jp][4];
			el.d2NdX2(n,1,2) = el.d2NdX2(n,2,1) = result[jp][5];
			el.d2NdX2(n,1,3) = el.d2NdX2(n,3,1) = result[jp][6];
			el.d2NdX2(n,2,2) = result[jp][7];
			el.d2NdX2(n,2,3) = el.d2NdX2(n,3,2) = result[jp][8];
			el.d2NdX2(n,3,3) = result[jp][9];
		}
	}
}

bool ASMu3D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
	std::cout << "YEAH BABY! We arrived at ASMu3D::integrate()" << std::endl;
	if (!lrspline) return true; // silently ignore empty patches

	PROFILE2("ASMu3D::integrate(I)");

	// Get Gaussian quadrature points and weights
	const double* xg = GaussQuadrature::getCoord(nGauss);
	const double* wg = GaussQuadrature::getWeight(nGauss);
	if (!xg || !wg) return false;

	// Get the reduced integration quadrature points, if needed
	const double* xr = 0;
	int nRed = integrand.getReducedIntegration();
	if (nRed < 0)
		nRed = nGauss; // The integrand needs to know nGauss
	else if (nRed > 0 && !(xr = GaussQuadrature::getCoord(nRed)))
		return false;


	// === Assembly loop over all elements in the patch ==========================

	bool ok=true;
	for(LR::Element *el : lrspline->getAllElements()) {
		int iEl = el->getId();
		int nBasis = el->nBasisFunctions();
		FiniteElement fe(nBasis);
		fe.iel = iEl;

		Vector   N;
		Matrix   dNdu, Xnod, Jac;
		Matrix3D d2Ndu2, Hess;
		double   dXidu[3];
		Vec4     X;
		// Get element volume in the parameter space
		double dV = el->volume();
		if (dV < 0.0)
		{
			ok = false; // topology error (probably logic error)
			break;
		}

		// Set up control point (nodal) coordinates for current element
		if (!this->getElementCoordinates(Xnod,iEl+1))
		{
			ok = false;
			break;
		}

		// Compute parameter values of the Gauss points over the whole element
		Vector gpar[3], redpar[3];
		for (int d = 0; d < 3; d++)
		{
			this->getGaussPointParameters(gpar[d],d,nGauss,iEl,xg);
			if (xr)
				this->getGaussPointParameters(redpar[d],d,nRed,iEl,xr);
		}


		if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
			this->getElementCorners(iEl, fe.XC);

		if (integrand.getIntegrandType() & Integrand::G_MATRIX)
		{
			// Element size in parametric space
			dXidu[0] = el->getParmin(0);
			dXidu[1] = el->getParmin(1);
			dXidu[2] = el->getParmin(2);
		}
		else if (integrand.getIntegrandType() & Integrand::AVERAGE)
		{
			// --- Compute average value of basis functions over the element -----

			fe.Navg.resize(nBasis,true);
			double vol = 0.0;
			for (int k = 0; k < nGauss; k++)
				for (int j = 0; j < nGauss; j++)
					for (int i = 0; i < nGauss; i++)
					{
						// Fetch basis function derivatives at current integration point
						evaluateBasis(fe, N, dNdu);

						// Compute Jacobian determinant of coordinate mapping
						// and multiply by weight of current integration point
						double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
						double weight = 0.125*dV*wg[i]*wg[j]*wg[k];

						// Numerical quadrature
						fe.Navg.add(fe.N,detJac*weight);
						vol += detJac*weight;
			}

			// Divide by element volume
			fe.Navg /= vol;
		}

		else if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
		{
			// Compute the element center
			Go::Point X0;
			double u0 = 0.5*(el->getParmin(0) + el->getParmax(0));
			double v0 = 0.5*(el->getParmin(1) + el->getParmax(1));
			double w0 = 0.5*(el->getParmin(2) + el->getParmax(2));
			lrspline->point(X0,u0,v0,w0);
			X = Vec3(X0);
		}

		// Initialize element quantities
		LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
		if (!integrand.initElement(MNPC[iEl],X,nRed*nRed*nRed,*A))
		{
			A->destruct();
			ok = false;
			break;
		}

		if (xr)
		{
			std::cerr << "Haven't really figured out what this part does yet\n";
			exit(42142);
			#if 0
			// --- Selective reduced integration loop ----------------------------

			int ip = (((i3-p3)*nRed*nel2 + i2-p2)*nRed*nel1 + i1-p1)*nRed;
			for (int k = 0; k < nRed; k++, ip += nRed*(nel2-1)*nRed*nel1)
				for (int j = 0; j < nRed; j++, ip += nRed*(nel1-1))
					for (int i = 0; i < nRed; i++, ip++)
					{
						// Local element coordinates of current integration point
						fe.xi   = xr[i];
						fe.eta  = xr[j];
						fe.zeta = xr[k];

						// Parameter values of current integration point
						fe.u = redpar[0](i+1,i1-p1+1);
						fe.v = redpar[1](j+1,i2-p2+1);
						fe.w = redpar[2](k+1,i3-p3+1);

						// Fetch basis function derivatives at current point
						evaluateBasis(fe, 1); 

						// Compute Jacobian inverse and derivatives
						fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

						// Cartesian coordinates of current integration point
						X = Xnod * fe.N;
						X.t = time.t;

						// Compute the reduced integration terms of the integrand
						if (!integrand.reducedInt(*A,fe,X))
							ok = false;
			}
			#endif
		}


		// --- Integration loop over all Gauss points in each direction --------

		fe.iGP = iEl*nGauss*nGauss; // Global integration point counter

		for (int k = 0; k < nGauss; k++) 
			for (int j = 0; j < nGauss; j++) 
				for (int i = 0; i < nGauss; i++, fe.iGP++)
				{
					// Local element coordinates of current integration point
					fe.xi   = xg[i];
					fe.eta  = xg[j];
					fe.zeta = xg[k];

					// Parameter values of current integration point
					fe.u = gpar[0](i+1);
					fe.v = gpar[1](j+1);
					fe.w = gpar[2](k+1);

					// Fetch basis function derivatives at current integration point
					if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
						evaluateBasis(fe, N, dNdu, d2Ndu2);
					else
						evaluateBasis(fe, N, dNdu) ;

					// Compute Jacobian inverse of coordinate mapping and derivatives
					fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
					if (fe.detJxW == 0.0) continue; // skip singular points

					// Compute Hessian of coordinate mapping and 2nd order derivatives
					if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
						if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
							ok = false;

					// Compute G-matrix
					if (integrand.getIntegrandType() & Integrand::G_MATRIX)
						utl::getGmat(Jac,dXidu,fe.G);

#if SP_DEBUG > 4
					std::cout <<"\niel, "<< iEl <<" "
					          <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX << std::endl;
#endif

					// Cartesian coordinates of current integration point
					X   = Xnod * fe.N;
					X.t = time.t;

					// Evaluate the integrand and accumulate element contributions
					fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
					if (!integrand.evalInt(*A,fe,time,X))
						ok = false;

		} // end gauss integrand

		// Finalize the element quantities
		if (ok && !integrand.finalizeElement(*A,time,0))
			ok = false;

		// Assembly of global system integral
		if (ok && !glInt.assemble(A->ref(),fe.iel))
			ok = false;

		A->destruct();
	}

	return ok;
}


bool ASMu3D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
	std::cerr << "ASMu3D::integrate(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	if (!lrspline) return true; // silently ignore empty patches

	PROFILE2("ASMu3D::integrate(B)");

	std::map<char,ThreadGroups>::const_iterator tit;
	if ((tit = threadGroupsFace.find(lIndex)) == threadGroupsFace.end())
	{
		std::cerr <<" *** ASMu3D::integrate: No thread groups for face "<< lIndex
				      << std::endl;
		return false;
	}
	const ThreadGroups& threadGrp = tit->second;

#ifdef USE_OPENMP
	int threads=omp_get_max_threads();
	// GoTools objects cannot be used in parallel sections
	if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
		omp_set_num_threads(1);
#endif

	// Get Gaussian quadrature points and weights
	int nGP = integrand.getBouIntegrationPoints(nGauss);
	const double* xg = GaussQuadrature::getCoord(nGP);
	const double* wg = GaussQuadrature::getWeight(nGP);
	if (!xg || !wg) return false;

	// Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
	const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

	const int t1 = 1 + abs(faceDir)%3; // first tangent direction
	const int t2 = 1 + t1%3;           // second tangent direction

	// Compute parameter values of the Gauss points over the whole patch face
	Matrix gpar[3];
	for (int d = 0; d < 3; d++)
		if (-1-d == faceDir)
		{
			gpar[d].resize(1,1);
			gpar[d].fill(lrspline->startparam(d));
		}
		else if (1+d == faceDir)
		{
			gpar[d].resize(1,1);
			gpar[d].fill(lrspline->endparam(d));
		}
		else
			this->getGaussPointParameters(gpar[d],d,nGP,xg);

	// Evaluate basis function derivatives at all integration points
	std::vector<Go::BasisDerivs> spline;
	{
		PROFILE2("Spline evaluation");
		lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
	}

	const int n1 = lrspline->numCoefs(0);
	const int n2 = lrspline->numCoefs(1);
	const int n3 = lrspline->numCoefs(2);

	const int p1 = lrspline->order(0);
	const int p2 = lrspline->order(1);
	const int p3 = lrspline->order(2);

	const int nel1 = n1 - p1 + 1;
	const int nel2 = n2 - p2 + 1;

	// Integrate the extraordinary elements?
	size_t doXelms = 0;
	if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
		if ((doXelms = (n1-p1+1)*(n2-p2+1)*(n3-p3+1))*2 > MNPC.size())
		{
			std::cerr <<" *** ASMu3D::integrate: Too few XO-elements "
				        << MNPC.size() - doXelms << std::endl;
			return false;
		}

	std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
	size_t firstp = iit == firstBp.end() ? 0 : iit->second;


	// === Assembly loop over all elements on the patch face =====================

	bool ok = true;
	for (size_t g = 0; g < threadGrp.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
		for (size_t t = 0; t < threadGrp[g].size(); ++t) {
			FiniteElement fe(p1*p2*p3);
			fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
			fe.u = gpar[0](1,1);
			fe.v = gpar[1](1,1);
			fe.w = gpar[2](1,1);

			Matrix dNdu, Xnod, Jac;
			Vec4   X;
			Vec3   normal;
			double dXidu[3];

			for (size_t l = 0; l < threadGrp[g][t].size() && ok; ++l)
			{
				int iel = threadGrp[g][t][l];
				fe.iel = MLGE[doXelms+iel];
				if (fe.iel < 1) continue; // zero-volume element

				int i1 = p1 + iel % nel1;
				int i2 = p2 + (iel / nel1) % nel2;
				int i3 = p3 + iel / (nel1*nel2);

				// Get element face area in the parameter space
				double dA = this->getParametricArea(++iel,abs(faceDir));
				if (dA < 0.0) // topology error (probably logic error)
				{
				  ok = false;
				  break;
				}

				// Set up control point coordinates for current element
				if (!this->getElementCoordinates(Xnod,iel))
				{
				  ok = false;
				  break;
				}

				if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
				  this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

				if (integrand.getIntegrandType() & Integrand::G_MATRIX)
				{
				  // Element size in parametric space
				  dXidu[0] = lrspline->knotSpan(0,i1-1);
				  dXidu[1] = lrspline->knotSpan(1,i2-1);
				  dXidu[2] = lrspline->knotSpan(2,i3-1);
				}

				// Initialize element quantities
				LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
				if (!integrand.initElementBou(MNPC[doXelms+iel-1],*A))
				{
				  A->destruct();
				  ok = false;
				  break;
				}

				// Define some loop control variables depending on which face we are on
				int nf1, j1, j2;
				switch (abs(faceDir))
				{
				  case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
				  case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
				  case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
				  default: nf1 = j1 = j2 = 0;
				}


				// --- Integration loop over all Gauss points in each direction --------

				int k1, k2, k3;
				int ip = (j2*nGP*nf1 + j1)*nGP;
				int jp = (j2*nf1 + j1)*nGP*nGP;
				fe.iGP = firstp + jp; // Global integration point counter

				for (int j = 0; j < nGP; j++, ip += nGP*(nf1-1))
				  for (int i = 0; i < nGP; i++, ip++, fe.iGP++)
				  {
				    // Local element coordinates and parameter values
				    // of current integration point
				    switch (abs(faceDir))
				    {
				      case 1: k2 = i+1; k3 = j+1; k1 = 0; break;
				      case 2: k1 = i+1; k3 = j+1; k2 = 0; break;
				      case 3: k1 = i+1; k2 = j+1; k3 = 0; break;
				      default: k1 = k2 = k3 = 0;
				    }
				    if (gpar[0].size() > 1)
				    {
				      fe.xi = xg[k1];
				      fe.u = gpar[0](k1,i1-p1+1);
				    }
				    if (gpar[1].size() > 1)
				    {
				      fe.eta = xg[k2];
				      fe.v = gpar[1](k2,i2-p2+1);
				    }
				    if (gpar[2].size() > 1)
				    {
				      fe.zeta = xg[k3];
				      fe.w = gpar[2](k3,i3-p3+1);
				    }

				    // Fetch basis function derivatives at current integration point
				    evaluateBasis(fe, N, dNdu);

				    // Compute basis function derivatives and the face normal
				    fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
				    if (fe.detJxW == 0.0) continue; // skip singular points

				    if (faceDir < 0) normal *= -1.0;

				    // Compute G-matrix
				    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
				      utl::getGmat(Jac,dXidu,fe.G);

				    // Cartesian coordinates of current integration point
				    X = Xnod * fe.N;
				    X.t = time.t;

				    // Evaluate the integrand and accumulate element contributions
				    fe.detJxW *= 0.25*dA*wg[i]*wg[j];
				    if (!integrand.evalBou(*A,fe,time,X,normal))
				      ok = false;
				  }

				// Assembly of global system integral
				if (ok && !glInt.assemble(A->ref(),fe.iel))
				  ok = false;

				A->destruct();
			}
		}
	}
#ifdef USE_OPENMP
	omp_set_num_threads(threads);
#endif

	return ok;
#endif
	return false;
}


bool ASMu3D::integrateEdge (Integrand& integrand, int lEdge,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
	std::cerr << "ASMu3D::integrateEdge(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	if (!lrspline) return true; // silently ignore empty patches

	PROFILE2("ASMu3D::integrate(E)");

	// Get Gaussian quadrature points and weights
	const double* xg = GaussQuadrature::getCoord(nGauss);
	const double* wg = GaussQuadrature::getWeight(nGauss);
	if (!xg || !wg) return false;

	// Compute parameter values of the Gauss points along the whole patch edge
	Matrix gpar[3];
	for (int d = 0; d < 3; d++)
		if (lEdge < d*4+1 || lEdge >= d*4+5)
		{
			gpar[d].resize(1,1);
			if (lEdge%4 == 1)
	gpar[d].fill(lrspline->startparam(d));
			else if (lEdge%4 == 0)
	gpar[d].fill(lrspline->endparam(d));
			else if (lEdge == 6 || lEdge == 10)
	gpar[d].fill(d == 0 ? lrspline->endparam(d) : lrspline->startparam(d));
			else if (lEdge == 2 || lEdge == 11)
	gpar[d].fill(d == 1 ? lrspline->endparam(d) : lrspline->startparam(d));
			else if (lEdge == 3 || lEdge == 7)
	gpar[d].fill(d == 2 ? lrspline->endparam(d) : lrspline->startparam(d));
		}
		else
		{
			int pm1 = lrspline->order(d) - 1;
			RealArray::const_iterator uit = lrspline->basis(d).begin() + pm1;
			double ucurr, uprev = *(uit++);
			int nCol = lrspline->numCoefs(d) - pm1;
			gpar[d].resize(nGauss,nCol);
			for (int j = 1; j <= nCol; uit++, j++)
			{
	ucurr = *uit;
	for (int i = 1; i <= nGauss; i++)
	  gpar[d](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
	uprev = ucurr;
			}
		}

	// Evaluate basis function derivatives at all integration points
	std::vector<Go::BasisDerivs> spline;
	{
		PROFILE2("Spline evaluation");
		lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
	}

	const int n1 = lrspline->numCoefs(0);
	const int n2 = lrspline->numCoefs(1);
	const int n3 = lrspline->numCoefs(2);

	const int p1 = lrspline->order(0);
	const int p2 = lrspline->order(1);
	const int p3 = lrspline->order(2);

	std::map<char,size_t>::const_iterator iit = firstBp.find(lEdge);
	size_t firstp = iit == firstBp.end() ? 0 : iit->second;

	FiniteElement fe(p1*p2*p3);
	fe.u = gpar[0](1,1);
	fe.v = gpar[1](1,1);
	fe.w = gpar[2](1,1);
	if (gpar[0].size() == 1) fe.xi = fe.u == lrspline->startparam(0) ? -1.0 : 1.0;
	if (gpar[1].size() == 1) fe.eta = fe.v == lrspline->startparam(1) ? -1.0 : 1.0;
	if (gpar[2].size() == 1) fe.zeta = fe.w == lrspline->startparam(2) ? -1.0 : 1.0;

	Matrix dNdu, Xnod, Jac;
	Vec4   X;
	Vec3   tang;


	// === Assembly loop over all elements on the patch edge =====================

	int iel = 1;
	for (int i3 = p3; i3 <= n3; i3++)
		for (int i2 = p2; i2 <= n2; i2++)
			for (int i1 = p1; i1 <= n1; i1++, iel++)
			{
	fe.iel = MLGE[iel-1];
	if (fe.iel < 1) continue; // zero-volume element

	// Skip elements that are not on current boundary edge
	bool skipMe = false;
	switch (lEdge)
	  {
	  case  1: if (i2 > p2 || i3 > p3) skipMe = true; break;
	  case  2: if (i2 < n2 || i3 > p3) skipMe = true; break;
	  case  3: if (i2 > p2 || i3 < n3) skipMe = true; break;
	  case  4: if (i2 < n2 || i3 < n3) skipMe = true; break;
	  case  5: if (i1 > p1 || i3 > p3) skipMe = true; break;
	  case  6: if (i1 < n1 || i3 > p3) skipMe = true; break;
	  case  7: if (i1 > p1 || i3 < n3) skipMe = true; break;
	  case  8: if (i1 < n1 || i3 < n3) skipMe = true; break;
	  case  9: if (i1 > p1 || i2 > p2) skipMe = true; break;
	  case 10: if (i1 < n1 || i2 > p2) skipMe = true; break;
	  case 11: if (i1 > p1 || i2 < n2) skipMe = true; break;
	  case 12: if (i1 < n1 || i2 < n2) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Get element edge length in the parameter space
	double dS = 0.0;
	int ip = MNPC[iel-1].back();
#ifdef INDEX_CHECK
	if (ip < 0 || (size_t)ip >= nodeInd.size()) return false;
#endif
	if (lEdge < 5)
	{
	  dS = lrspline->knotSpan(0,nodeInd[ip].I);
	  ip = (i1-p1)*nGauss;
	}
	else if (lEdge < 9)
	{
	  dS = lrspline->knotSpan(1,nodeInd[ip].J);
	  ip = (i2-p2)*nGauss;
	}
	else if (lEdge < 13)
	{
	  dS = lrspline->knotSpan(2,nodeInd[ip].K);
	  ip = (i3-p3)*nGauss;
	}

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
				LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
				bool ok = integrand.initElementBou(MNPC[iel-1],*A);


	// --- Integration loop over all Gauss points along the edge -----------

	fe.iGP = firstp + ip; // Global integration point counter

	for (int i = 0; i < nGauss && ok; i++, ip++, fe.iGP++)
	{
	  // Parameter values of current integration point
	  if (gpar[0].size() > 1) fe.u = gpar[0](i+1,i1-p1+1);
	  if (gpar[1].size() > 1) fe.v = gpar[1](i+1,i2-p2+1);
	  if (gpar[2].size() > 1) fe.w = gpar[2](i+1,i3-p3+1);

	  // Fetch basis function derivatives at current integration point
	  evaluateBasis(fe, N, dNdu);

	  // Compute basis function derivatives and the edge tang
	  fe.detJxW = utl::Jacobian(Jac,tang,fe.dNdX,Xnod,dNdu,1+(lEdge-1)/4);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= 0.5*dS*wg[i];
				  ok = integrand.evalBou(*A,fe,time,X,tang);
	}

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
	  ok = false;

	A->destruct();

	if (!ok) return false;
			}

	return true;
#endif
	return false;
}


int ASMu3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
	std::cerr << "ASMu3D::evalPoint(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	if (!lrspline) return -3;

	int i;
	for (i = 0; i < 3; i++)
		param[i] = (1.0-xi[i])*lrspline->startparam(i) + xi[i]*lrspline->endparam(i);

	Go::Point X0;
	lrspline->point(X0,param[0],param[1],param[2]);
	for (i = 0; i < 3 && i < lrspline->dimension(); i++)
		X[i] = X0[i];

	// Check if this point matches any of the control points (nodes)
	Vec3 Xnod;
	size_t inod = 1;
	RealArray::const_iterator cit = lrspline->coefs_begin();
	for (i = 0; cit != lrspline->coefs_end(); cit++, i++)
	{
		if (i < 3) Xnod[i] = *cit;
		if (i+1 == lrspline->dimension())
			if (X.equal(Xnod,0.001))
	return inod;
			else
			{
	inod++;
	i = -1;
			}
	}

	return 0;
#endif
	return -1;
}


bool ASMu3D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
	std::cerr << "ASMu3D::getGridParameters not implemented properly yet" << std::endl;
	exit(776654);
	#if 0
	if (!lrspline) return false;

	if (nSegPerSpan < 1)
	{
		std::cerr <<" *** ASMu3D::getGridParameters: Too few knot-span points "
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

	#endif

	return true;
}

bool ASMu3D::tesselate (ElementBlock& grid, const int* npe) const
{
	std::cerr << "ASMu3D::tesselate(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	// Compute parameter values of the nodal points
	RealArray gpar[3];
	for (int dir = 0; dir < 3; dir++)
		if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
			return false;

	// Evaluate the spline volume at all points
	size_t nx = gpar[0].size();
	size_t ny = gpar[1].size();
	size_t nz = gpar[2].size();
	RealArray XYZ(lrspline->dimension()*nx*ny*nz);
	lrspline->gridEvaluator(gpar[0],gpar[1],gpar[2],XYZ);

	// Establish the block grid coordinates
	size_t i, j, k, l;
	grid.resize(nx,ny,nz);
	for (i = j = 0; i < grid.getNoNodes(); i++, j += lrspline->dimension())
		grid.setCoor(i,XYZ[j],XYZ[j+1],XYZ[j+2]);

	// Establish the block grid topology
	int ie, nse1 = npe[0] - 1;
	int je, nse2 = npe[1] - 1;
	int ke, nse3 = npe[2] - 1;
	int nel1 = (nx-1)/nse1;
	int nel2 = (ny-1)/nse2;
	int n[8], ip = 0;
	for (k = ke = 1, n[2] = 0; k < nz; k++)
	{
		for (j = je = 1, n[1] = n[2]; j < ny; j++)
		{
			n[0] = n[1];
			n[1] = n[0] + 1;
			n[2] = n[1] + nx;
			n[3] = n[1] + nx-1;
			n[4] = n[0] + nx*ny;
			n[5] = n[4] + 1;
			n[6] = n[5] + nx;
			n[7] = n[5] + nx-1;
			for (i = ie = 1; i < nx; i++)
			{
	for (l = 0; l < 8; l++)
	  grid.setNode(ip++,n[l]++);
	grid.setElmId(((k-1)*(ny-1)+j-1)*(nx-1)+i,((ke-1)*nel2+je-1)*nel1+ie);
	if (i%nse1 == 0) ie++;
			}
			if (j%nse2 == 0) je++;
		}
		if (k%nse3 == 0) ke++;
	}

	return true;
#endif
	return false;
}


void ASMu3D::scatterInd (int n1, int n2, int n3, int p1, int p2, int p3,
			 const int* start, IntVec& index)
{
	std::cerr << "ASMu3D::scatterInd(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	index.reserve(p1*p2*p3);
	int ip = ((start[2]-p3+1)*n2 + (start[1]-p2+1))*n1 + (start[0]-p1+1);
	for (int i3 = 0; i3 < p3; i3++, ip += n1*(n2-p2))
		for (int i2 = 0; i2 < p2; i2++, ip += n1-p1)
			for (int i1 = 0; i1 < p1; i1++, ip++)
				index.push_back(ip);
			#endif
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
			   const int* npe) const
{
	std::cerr << "ASMu3D::evalSolution(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	// Compute parameter values of the result sampling points
	RealArray gpar[3];
	for (int dir = 0; dir < 3; dir++)
		if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
			return false;

	// Evaluate the primary solution at all sampling points
	return this->evalSolution(sField,locSol,gpar);
#endif
	return false;
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
			   const RealArray* gpar, bool regular) const
{
	std::cerr << "ASMu3D::evalSolution(...) is not properly implemented yet :(" << std::endl;
	exit(776672);
#if 0
	// Evaluate the basis functions at all points
	std::vector<Go::BasisPts> spline;
	if (regular)
	{
		PROFILE2("Spline evaluation");
		lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
	}
	else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
	{
		PROFILE2("Spline evaluation");
		spline.resize(gpar[0].size());
		for (size_t i = 0; i < spline.size(); i++)
			lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline[i]);
	}
	else
		return false;

	const int p1 = lrspline->order(0);
	const int p2 = lrspline->order(1);
	const int p3 = lrspline->order(2);
	const int n1 = lrspline->numCoefs(0);
	const int n2 = lrspline->numCoefs(1);
	const int n3 = lrspline->numCoefs(2);
//  size_t nComp = locSol.size() / this->getNoNodes(-1);
//  if (nComp*this->getNoNodes(-1) != locSol.size())
//    return false;
	size_t nComp = locSol.size() / (n1*n2*n3);

	Matrix Xtmp;
	Vector Ytmp;

	// Evaluate the primary solution field at each point
	size_t nPoints = spline.size();
	sField.resize(nComp,nPoints);
	for (size_t i = 0; i < nPoints; i++)
	{
		IntVec ip;
		scatterInd(n1,n2,n3,p1,p2,p3,spline[i].left_idx,ip);

		utl::gather(ip,nComp,locSol,Xtmp);
		Xtmp.multiply(spline[i].basisValues,Ytmp);
		sField.fillColumn(1+i,Ytmp);
	}

	return true;
#endif
	return false;
}


bool ASMu3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const int* npe, char project) const
{
	std::cerr << "ASMu3D::evalSolution(...) is not properly implemented yet :(" << std::endl;
	exit(776658);
#if 0
	// Project the secondary solution onto the spline basis
	Go::SplineVolume* v = NULL;
	if (project == 'A')
		v = this->projectSolutionLocalApprox(integrand);
	else if (project == 'L')
		v = this->projectSolutionLocal(integrand);
	else if (project == 'W')
		v = this->projectSolutionLeastSquare(integrand);
	else if (project == 'D' || !npe)
		v = this->projectSolution(integrand);

	if (npe)
	{
		// Compute parameter values of the result sampling points
		RealArray gpar[3];
		if (this->getGridParameters(gpar[0],0,npe[0]-1) &&
	this->getGridParameters(gpar[1],1,npe[1]-1) &&
	this->getGridParameters(gpar[2],2,npe[2]-1))
		{
			if (!project)
				// Evaluate the secondary solution directly at all sampling points
				return this->evalSolution(sField,integrand,gpar);
			else if (v)
			{
				// Evaluate the projected field at the result sampling points
				const Vector& svec = sField; // using utl::matrix cast operator
				sField.resize(v->dimension(),
				              gpar[0].size()*gpar[1].size()*gpar[2].size());
				v->gridEvaluator(gpar[0],gpar[1],gpar[2],const_cast<Vector&>(svec));
				delete v;
				return true;
			}
		}
		else if (v)
			delete v;
	}
	else if (v)
	{
		// Extract control point values from the spline object
		sField.resize(v->dimension(),
				          v->numCoefs(0)*v->numCoefs(1)*v->numCoefs(2));
		sField.fill(&(*v->coefs_begin()));
		delete v;
		return true;
	}

	std::cerr <<" *** ASMu3D::evalSolution: Failure!";
	if (project) std::cerr <<" project="<< project;
	std::cerr << std::endl;
	return false;
#endif
	return false;
}


bool ASMu3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const RealArray* gpar, bool regular) const
{
	std::cerr << "ASMu3D::evalSolution(...) is not properly implemented yet :(" << std::endl;
	exit(776655);
#if 0
	sField.resize(0,0);

	// Evaluate the basis functions and their derivatives at all points
	size_t nPoints = gpar[0].size();
	bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
	std::vector<Go::BasisDerivs>  spline1(regular ||  use2ndDer ? 0 : nPoints);
	std::vector<Go::BasisDerivs2> spline2(regular || !use2ndDer ? 0 : nPoints);
	if (regular)
	{
		PROFILE2("Spline evaluation");
		nPoints *= gpar[1].size()*gpar[2].size();
		if (use2ndDer)
			lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
		else
			lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
	}
	else if (nPoints == gpar[1].size() && nPoints == gpar[2].size())
	{
		PROFILE2("Spline evaluation");
		for (size_t i = 0; i < nPoints; i++)
			if (use2ndDer)
				lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
			else
				lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
	}
	else
		return false;

	const int p1 = lrspline->order(0);
	const int p2 = lrspline->order(1);
	const int p3 = lrspline->order(2);
	const int n1 = lrspline->numCoefs(0);
	const int n2 = lrspline->numCoefs(1);
	const int n3 = lrspline->numCoefs(2);

	// Fetch nodal (control point) coordinates
	Matrix Xnod, Xtmp;
	this->getNodalCoordinates(Xnod);

	FiniteElement fe(p1*p2*p3);
	Vector        solPt;
	Matrix        dNdu, Jac;
	Matrix3D      d2Ndu2, Hess;

	// Evaluate the secondary solution field at each point
	for (size_t i = 0; i < nPoints; i++)
	{
		// Fetch indices of the non-zero basis functions at this point
		IntVec ip;
		if (use2ndDer)
		{
			scatterInd(n1,n2,n3,p1,p2,p3,spline2[i].left_idx,ip);
			fe.u = spline2[i].param[0];
			fe.v = spline2[i].param[1];
			fe.w = spline2[i].param[2];
		}
		else
		{
			scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);
			fe.u = spline1[i].param[0];
			fe.v = spline1[i].param[1];
			fe.w = spline1[i].param[2];
		}
		fe.iGP = firstIp + i;

		// Fetch associated control point coordinates
		utl::gather(ip,3,Xnod,Xtmp);

		// Fetch basis function derivatives at current integration point
		if (use2ndDer)
			evaluateBasis(fe, N, dNdu, d2Ndu2);
		else
			evaluateBasis(fe, N, dNdu);

		// Compute the Jacobian inverse and derivatives
		if (utl::Jacobian(Jac,fe.dNdX,Xtmp,dNdu) == 0.0) // Jac = (Xtmp * dNdu)^-1
			continue; // skip singular points

		// Compute Hessian of coordinate mapping and 2nd order derivatives
		if (use2ndDer)
			if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xtmp,d2Ndu2,dNdu))
				continue;

		// Now evaluate the solution field
		if (!integrand.evalSol(solPt,fe,Xtmp*fe.N,ip))
			return false;
		else if (sField.empty())
			sField.resize(solPt.size(),nPoints,true);

		sField.fillColumn(1+i,solPt);
	}

	return true;
#endif
	return false;
}


bool ASMu3D::evaluate (const ASMbase* input, const Vector& locVec, Vector& vec)
{
	std::cerr << "ASMu3D::evaluate(...) is not properly implemented yet :(" << std::endl;
	exit(776654);
#if 0
	ASMu3D* pch = (ASMu3D*)input;
	// Compute parameter values of the result sampling points (Greville points)
	RealArray gpar[3];
	for (int dir = 0; dir < 3; dir++)
		if (!this->getGrevilleParameters(gpar[dir],dir))
			return false;

	Matrix sValues;
	pch->evalSolution(sValues, locVec, gpar, true);

	// Project the results onto the spline basis to find control point
	// values based on the result values evaluated at the Greville points.
	// Note that we here implicitly assume that the number of Greville points
	// equals the number of control points such that we don't have to resize
	// the result array. Think that is always the case, but beware if trying
	// other projection schemes later.

	RealArray weights;
	if (lrspline->rational())
		lrspline->getWeights(weights);

	const Vector& vec2 = sValues;
	Go::SplineVolume* vol_new =
	          Go::VolumeInterpolator::regularInterpolation(lrspline->basis(0),
	                                                       lrspline->basis(1),
	                                                       lrspline->basis(2),
	                                                       gpar[0], gpar[1], gpar[2],
	                                                       const_cast<Vector&>(vec2),
	                                                       sValues.rows(),
	                                                       lrspline->rational(),
	                                                       weights);
	vec.resize(vol_new->coefs_end()-vol_new->coefs_begin());
	std::copy(vol_new->coefs_begin(), vol_new->coefs_end(), vec.begin());
	delete vol_new;

	return true;
#endif
	return false;
}

bool ASMu3D::addXElms (short int dim, short int item, size_t nXn, IntVec& nodes)
{
	std::cerr << "ASMu3D::addXElms not implemented properly yet\n";
	exit(776654);
#if 0
	if (dim != 2)
	{
		std::cerr <<" *** ASMu3D::addXElms: Invalid boundary dimension "<< dim
				      <<", only 2 (face) is allowed."<< std::endl;
		return false;
	}
	else if (!svol || shareFE)
		return false;

	for (size_t i = 0; i < nXn; i++)
	{
		if (nodes.size() == i)
			nodes.push_back(++gNod);
		myMLGN.push_back(nodes[i]);
	}

	const int n1 = svol->numCoefs(0);
	const int n2 = svol->numCoefs(1);
	const int n3 = svol->numCoefs(2);

	const int p1 = svol->order(0);
	const int p2 = svol->order(1);
	const int p3 = svol->order(2);

	nXelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);
	myMNPC.resize(2*nXelm);
	myMLGE.resize(2*nXelm,0);

	int iel = 0;
	bool skipMe = false;
	for (int i3 = p3; i3 <= n3; i3++)
		for (int i2 = p2; i2 <= n2; i2++)
			for (int i1 = p1; i1 <= n1; i1++, iel++)
			{
				if (MLGE[iel] < 1) continue; // Skip zero-volume element

				// Skip elements that are not on current boundary face
				switch (item)
				  {
				  case 1: skipMe = i1 > p1; break;
				  case 2: skipMe = i1 < n1; break;
				  case 3: skipMe = i2 > p2; break;
				  case 4: skipMe = i2 < n2; break;
				  case 5: skipMe = i3 > p3; break;
				  case 6: skipMe = i3 < n3; break;
				  }
				if (skipMe) continue;

				IntVec& mnpc = myMNPC[nXelm+iel];
				if (!mnpc.empty())
				{
				  std::cerr <<" *** ASMu3D::addXElms: Only one X-face allowed."
				            << std::endl;
				  return false;
				}

				mnpc = MNPC[iel]; // Copy the ordinary element nodes

				// Negate node numbers that are not on the boundary face, to flag that
				// they shall not receive any tangent and/or residual contributions
				int lnod = 0;
				for (int j3 = 0; j3 < p3; j3++)
				  for (int j2 = 0; j2 < p2; j2++)
				    for (int j1 = 0; j1 < p1; j1++, lnod++)
				    {
				      switch (item)
				        {
				        case 1: skipMe = j1 > 0;    break;
				        case 2: skipMe = j1 < p1-1; break;
				        case 3: skipMe = j2 > 0;    break;
				        case 4: skipMe = j2 < p2-1; break;
				        case 5: skipMe = j3 > 0;    break;
				        case 6: skipMe = j3 < p3-1; break;
	        }
	      if (skipMe) // Hack for node 0: Using -maxint as flag instead
		mnpc[lnod] = mnpc[lnod] == 0 ? -2147483648 : -mnpc[lnod];
	    }

				// Add connectivity to the extra-ordinary nodes
				for (size_t i = 0; i < nXn; i++)
				  mnpc.push_back(MLGN.size()-nXn+i);

	myMLGE[nXelm+iel] = ++gEl;
			}

#endif

	return true;

}

