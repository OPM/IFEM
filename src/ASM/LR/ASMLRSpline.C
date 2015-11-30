// $Id$
//==============================================================================
//!
//! \file ASMLRSpline.C
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#include "ASMunstruct.h"
#include "Profiler.h"
#include "LRSpline/LRSplineSurface.h"
#include <fstream>


ASMunstruct::~ASMunstruct ()
{
	if (geo) delete geo;
}


bool ASMunstruct::refine (const RealArray& elementError,
                          const IntVec& options,
                          Vectors* sol, const char* fName)
{
	PROFILE2("ASMunstruct::refine()");

	if (!geo) return false;
	if (shareFE) return true;

	// to pick up if LR splines get stuck while doing refinement print entry and exit point of this function

	double beta          = (options.size()>0)  ? options[0]/100.0 : 0.10;
	int    multiplicity  = (options.size()>1)  ? options[1]       : 1;
	enum refinementStrategy strat = LR_FULLSPAN;
	if (options.size() > 2)
	  switch (options[2]) {
	  case 1: strat = LR_MINSPAN; break;
	  case 2: strat = LR_STRUCTURED_MESH; break;
	  }
	bool   linIndepTest  = (options.size()>3)  ? options[3]!=0    : false;
	int    maxTjoints    = (options.size()>4)  ? options[4]       : -1;
	double maxAspectRatio= (options.size()>5)  ? options[5]       : -1;
	bool   closeGaps     = (options.size()>6)  ? options[6]!=0    : false;
	bool   isVol         = geo->nVariate()==3;

	if (multiplicity > 1)
	{
		int p1 = geo->order(0) - 1;
		int p2 = geo->order(1) - 1;
		multiplicity = options[1];
		if (multiplicity > p1) multiplicity = p1;
		if (multiplicity > p2) multiplicity = p2;
		if (isVol && multiplicity > geo->order(2)) multiplicity = geo->order(2);
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

	if (fName)
	{
		char fullFileName[256];

		strcpy(fullFileName, "lrspline_");
		strcat(fullFileName, fName);
		std::ofstream lrOut(fullFileName);
		lrOut << *geo;
		lrOut.close();

		if(!isVol) {
			LR::LRSplineSurface *lr = dynamic_cast<LR::LRSplineSurface*>(geo);
			// open files for writing
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

			lr->writePostscriptMesh(paramMeshFile);
			lr->writePostscriptElements(physicalMeshFile);
			lr->writePostscriptFunctionSpace(paramDotMeshFile);
			lr->writePostscriptMeshWithControlPoints(physicalDotMeshFile);

			// close all files
			paramMeshFile.close();
			physicalMeshFile.close();
			paramDotMeshFile.close();
			physicalDotMeshFile.close();

		}
	}

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
#ifdef HAS_BOOST
			std::cout << "Testing for linear independence by full tensor expansion " << std::endl;
			isLinIndep = geo->isLinearIndepByMappingMatrix(false);
#endif
		}
		if(!isLinIndep)
		{
			std::cout << "FAILED!!!" << std::endl;
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
			std::cout << "...Passed" << std::endl;
		}
	}

	return true;
}


bool ASMunstruct::refine (const IntVec& elements,
                          const IntVec& options,
                          Vectors* sol, const char* fName)
{
	PROFILE2("ASMunstruct::refine()");

	if (!geo) return false;
	if (shareFE) return true;

	// to pick up if LR splines get stuck while doing refinement print entry and exit point of this function

	int    multiplicity  = (options.size()>1)  ? options[1]       : 1;
	enum refinementStrategy strat = LR_FULLSPAN;
	if (options.size() > 2)
	  switch (options[2]) {
	  case 1: strat = LR_MINSPAN; break;
	  case 2: strat = LR_STRUCTURED_MESH; break;
	  }
	bool   linIndepTest  = (options.size()>3)  ? options[3]!=0    : false;
	int    maxTjoints    = (options.size()>4)  ? options[4]       : -1;
	double maxAspectRatio= (options.size()>5)  ? options[5]       : -1;
	bool   closeGaps     = (options.size()>6)  ? options[6]!=0    : false;
	bool   isVol         = geo->nVariate()==3;

	if (multiplicity > 1)
	{
		int p1 = geo->order(0) - 1;
		int p2 = geo->order(1) - 1;
		multiplicity = options[1];
		if (multiplicity > p1) multiplicity = p1;
		if (multiplicity > p2) multiplicity = p2;
		if (isVol && multiplicity > geo->order(2)) multiplicity = geo->order(2);
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
		if(strat == LR_STRUCTURED_MESH)
			geo->refineBasisFunction(elements);
		else
			geo->refineElement(elements);
	}
	if (fName)
	{
		char fullFileName[256];

		strcpy(fullFileName, "lrspline_");
		strcat(fullFileName, fName);
		std::ofstream lrOut(fullFileName);
		lrOut << *geo;
		lrOut.close();

		if(!isVol) {
			LR::LRSplineSurface *lr = dynamic_cast<LR::LRSplineSurface*>(geo);
			// open files for writing
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

			lr->writePostscriptMesh(paramMeshFile);
			lr->writePostscriptElements(physicalMeshFile);
			lr->writePostscriptFunctionSpace(paramDotMeshFile);
			lr->writePostscriptMeshWithControlPoints(physicalDotMeshFile);

			// close all files
			paramMeshFile.close();
			physicalMeshFile.close();
			paramDotMeshFile.close();
			physicalDotMeshFile.close();

		}
	}

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
#ifdef HAS_BOOST
			std::cout << "Testing for linear independence by full tensor expansion " << std::endl;
			isLinIndep = geo->isLinearIndepByMappingMatrix(false);
#endif
		}
		if(!isLinIndep)
		{
			std::cout << "FAILED!!!" << std::endl;
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
			std::cout << "...Passed" << std::endl;
		}
	}

	return true;
}

Go::BsplineBasis ASMunstruct::getBezierBasis (int p)
{
	double knot[2*p];
	std::fill(knot,   knot+p,  -1.0);
	std::fill(knot+p, knot+2*p, 1.0);
	return Go::BsplineBasis(p,p,knot);
}
