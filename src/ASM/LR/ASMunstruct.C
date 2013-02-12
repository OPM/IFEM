// $Id$
//==============================================================================
//!
//! \file ASMunstruct.C
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#include "Profiler.h"
#include "ASMunstruct.h"
#include "GoTools/geometry/GeomObject.h"
#include "LRSpline/LRSplineSurface.h"


int ASMunstruct::gEl = 0;
int ASMunstruct::gNod = 0;


ASMunstruct::ASMunstruct (unsigned char n_p, unsigned char n_s,
			  unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geo = 0;
}


ASMunstruct::ASMunstruct (const ASMunstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  nGauss = patch.nGauss;
  geo = patch.geo;
}

bool ASMunstruct::refine (const std::vector<double>& elementError,
                     const std::vector<int>&    options,
                     const char* fName)
{
	PROFILE2("ASMunstruct::refine()");
	std::cout << "ASMunstruct::refine()" << std::endl;

	if (!geo) return false;
	if (shareFE) return true;

	std::cout << "passed geo and shareFE" << std::endl;

	double                  beta          = (options.size()>0)  ? options[0]/100.0 : 0.10;
	int                     multiplicity  = (options.size()>1)  ? options[1]       : 1;
	enum refinementStrategy strat         = LR_SAFE;
	bool                    linIndepTest  = (options.size()>3)  ? options[3]!=0    : false;
	int                     maxTjoints    = (options.size()>4)  ? options[4]       : -1;
	double                  maxAspectRatio= (options.size()>5)  ? options[5]       : -1;
	bool                    closeGaps     = (options.size()>6)  ? options[6]!=0    : false;

	if(options.size() > 2) {
		if(options[2]==1)      strat  = LR_MINSPAN;
		else if(options[2]==2) strat  = LR_ISOTROPIC_FUNC;
	}

	if (multiplicity > 1)
	{
		int p1 = geo->order(0) - 1;
		int p2 = geo->order(1) - 1;
		multiplicity = options[1];
		if (multiplicity > p1) multiplicity = p1;
		if (multiplicity > p2) multiplicity = p2;
	}

	if (!elementError.empty()) {
		// set refinement parameters
		if(maxTjoints > 0)
			geo->setMaxTjoints(maxTjoints);
		if(maxAspectRatio > 0)
			geo->setMaxAspectRatio(maxAspectRatio);
		geo->setCloseGaps(closeGaps);
		geo->setRefMultiplicity(multiplicity);
		geo->setRefStrat(strat);

		// do actual refinement
		geo->refineByDimensionIncrease(elementError, beta);
	}
	/*
	if (fName)
	{
		char fullFileName[256];

		strcpy(fullFileName, "param_");
		strcat(fullFileName, fName);
		std::ofstream paramMeshFile(fullFileName);

		strcpy(fullFileName, "physical_");
		strcat(fullFileName, fName);
		std::ofstream physicalMeshFile(fullFileName);

		strcpy(fullFileName, "param_dot_");
		strcat(fullFileName, fName);
		std::ofstream paramDotMeshFile(fullFileName);

		strcpy(fullFileName, "physical_dot_");
		strcat(fullFileName, fName);
		std::ofstream physicalDotMeshFile(fullFileName);

		strcpy(fullFileName, "lrspline_");
		strcat(fullFileName, fName);
		std::ofstream lrOut(fullFileName);

		strcpy(fullFileName, "refine_details_");
		strcat(fullFileName, fName);
		std::ofstream refineDetails(fullFileName);

		geo->writePostscriptMesh(paramMeshFile);
		geo->writePostscriptElements(physicalMeshFile);
		geo->writePostscriptFunctionSpace(paramDotMeshFile);
		geo->writePostscriptMeshWithControlPoints(physicalDotMeshFile);
		lrOut << *geo;
		refineDetails  << beta << " " << multiplicity << " " << strat << " "      << std::endl;
		refineDetails  << maxTjoints << " " << maxAspectRatio << " " << closeGaps << std::endl;
		refineDetails  << elementError.size() << std::endl;
		for(size_t i=0; i<elementError.size(); i++) {
			refineDetails << elementError[i] << std::endl;
		}
	}
	*/

	if (!elementError.empty())
		std::cout <<"Refined mesh: "<< geo->nElements()
		          <<" elements "<< geo->nBasisFunctions()
		          <<" nodes."<< std::endl;

	if(linIndepTest)
	{
		std::cout << "Testing for linear independence by overloading " << std::endl;
		bool isLinIndep = geo->isLinearIndepByOverloading(false);
		if(!isLinIndep) {
			std::cout << "Inconclusive..." << std::endl;
			std::cout << "Testing for linear independence by full tensor expansion " << std::endl;
			isLinIndep = geo->isLinearIndepByMappingMatrix(false);
		}
		if(!isLinIndep)
		{
			std::cout << "FAILED!!!\n";
			/*
			std::cerr << std::endl;
			std::cerr << std::endl;
			std::cerr << std::endl;
			std::cerr << "*****************************************************************\n";
			std::cerr << "\nLR B-spline is linear dependant. Continuing analysis anyway\n\n";
			std::cerr << "*****************************************************************\n";
			std::cerr << std::endl;
			std::cerr << std::endl;
			*/
			exit(228);
			// return false;
		}
		else
		{
			std::cout << "...Passed\n";
		}
	}

	// int nBasis    = geo->nBasisFunctions();
	// myMLGN.resize(nBasis);
	// for (int inod = 0; inod < nBasis; inod++)
		// myMLGN[inod] = ++gNod;

	return true;
}

bool ASMunstruct::refine (const std::vector<int>& elements,
                     const std::vector<int>& options,
                     const char* fName)
{
	PROFILE2("ASMunstruct::refine()");
	std::cout << "ASMunstruct::refine()" << std::endl;

	if (!geo) return false;
	if (shareFE) return true;

	std::cout << "passed geo and shareFE" << std::endl;

	double                  beta          = (options.size()>0)  ? options[0]/100.0 : 0.10;
	int                     multiplicity  = (options.size()>1)  ? options[1]       : 1;
	enum refinementStrategy strat         = LR_SAFE;
	bool                    linIndepTest  = (options.size()>3)  ? options[3]!=0    : false;
	int                     maxTjoints    = (options.size()>4)  ? options[4]       : -1;
	double                  maxAspectRatio= (options.size()>5)  ? options[5]       : -1;
	bool                    closeGaps     = (options.size()>6)  ? options[6]!=0    : false;
//        bool                    trueBeta      = (options.size()>7)  ? options[7]!=0    : false;

	if(options.size() > 2) {
		if(options[2]==1)      strat  = LR_MINSPAN;
		else if(options[2]==2) strat  = LR_ISOTROPIC_FUNC;
	}

	if (multiplicity > 1)
	{
		int p1 = geo->order(0) - 1;
		int p2 = geo->order(1) - 1;
		multiplicity = options[1];
		if (multiplicity > p1) multiplicity = p1;
		if (multiplicity > p2) multiplicity = p2;
	}

	if (!elements.empty()) {
		// set refinement parameters
		if(maxTjoints > 0)
			geo->setMaxTjoints(maxTjoints);
		if(maxAspectRatio > 0)
			geo->setMaxAspectRatio(maxAspectRatio);
		geo->setCloseGaps(closeGaps);
		geo->setRefMultiplicity(multiplicity);
		geo->setRefStrat(strat);

		// do actual refinement
		if(strat == LR_ISOTROPIC_FUNC)
			geo->refineBasisFunction(elements);
		else
			geo->refineElement(elements);
// 		if(trueBeta)
// 			geo->refineByDimensionIncrease(elements, beta); 
	}
	/*
	if (fName)
	{
		char fullFileName[256];

		strcpy(fullFileName, "param_");
		strcat(fullFileName, fName);
		std::ofstream paramMeshFile(fullFileName);

		strcpy(fullFileName, "physical_");
		strcat(fullFileName, fName);
		std::ofstream physicalMeshFile(fullFileName);

		strcpy(fullFileName, "param_dot_");
		strcat(fullFileName, fName);
		std::ofstream paramDotMeshFile(fullFileName);

		strcpy(fullFileName, "physical_dot_");
		strcat(fullFileName, fName);
		std::ofstream physicalDotMeshFile(fullFileName);

		strcpy(fullFileName, "lrspline_");
		strcat(fullFileName, fName);
		std::ofstream lrOut(fullFileName);

		strcpy(fullFileName, "refine_details_");
		strcat(fullFileName, fName);
		std::ofstream refineDetails(fullFileName);

		geo->writePostscriptMesh(paramMeshFile);
		geo->writePostscriptElements(physicalMeshFile);
		geo->writePostscriptFunctionSpace(paramDotMeshFile);
		geo->writePostscriptMeshWithControlPoints(physicalDotMeshFile);
		lrOut << *geo;
		refineDetails  << beta << " " << multiplicity << " " << strat << " "      << std::endl;
		refineDetails  << maxTjoints << " " << maxAspectRatio << " " << closeGaps << std::endl;
		refineDetails  << elements.size() << std::endl;
		for(size_t i=0; i<elements.size(); i++) {
			refineDetails << elements[i] << std::endl;
		}
	}
	*/

	if (!elements.empty())
		std::cout <<"Refined mesh: "<< geo->nElements()
		          <<" elements "<< geo->nBasisFunctions()
		          <<" nodes."<< std::endl;

	if(linIndepTest)
	{
		std::cout << "Testing for linear independence by overloading " << std::endl;
		bool isLinIndep = geo->isLinearIndepByOverloading(false);
		if(!isLinIndep) {
			std::cout << "Inconclusive..." << std::endl;
			std::cout << "Testing for linear independence by full tensor expansion " << std::endl;
			isLinIndep = geo->isLinearIndepByMappingMatrix(false);
		}
		if(!isLinIndep)
		{
			std::cout << "FAILED!!!\n";
			/*
			std::cerr << std::endl;
			std::cerr << std::endl;
			std::cerr << std::endl;
			std::cerr << "*****************************************************************\n";
			std::cerr << "\nLR B-spline is linear dependant. Continuing analysis anyway\n\n";
			std::cerr << "*****************************************************************\n";
			std::cerr << std::endl;
			std::cerr << std::endl;
			*/
			exit(228);
			// return false;
		}
		else
		{
			std::cout << "...Passed\n";
		}
	}

	// int nBasis    = geo->nBasisFunctions();
	// myMLGN.resize(nBasis);
	// for (int inod = 0; inod < nBasis; inod++)
		// myMLGN[inod] = ++gNod;

	return true;
}


ASMunstruct::~ASMunstruct ()
{
  if (geo) delete geo;
}

