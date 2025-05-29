// $Id$
//==============================================================================
//!
//! \file SIMsupel.C
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for general superelement FEM analysis.
//!
//==============================================================================

#include "SIMsupel.h"
#include "ASMbase.h"
#include "ASM3D.h"
#include "IntegrandBase.h"
#include "SAM.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMsupel::SIMsupel (const char* hd, char nf) : ncmp(abs(nf))
{
  // Empty integrand class for superelement analysis.
  class NoProblem : public IntegrandBase
  {
  public:
    explicit NoProblem(unsigned char n) : IntegrandBase(n) {}
  };

  if (nf > 0) myProblem = new NoProblem(nf);
  if (hd) myHeading = hd;
}


bool SIMsupel::parse (const tinyxml2::XMLElement* elem)
{
  if (strncasecmp(elem->Value(),"superel",7))
    return this->SIMgeneric::parse(elem);

  bool result = true;
  size_t last = myModel.size();
  SuperElm sup;
  std::string supNodeSet;
  utl::getAttribute(elem,"id",sup.id);
  utl::getAttribute(elem,"nodeset",supNodeSet);
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"mvp") && child->FirstChild())
    {
      sup.MVP.resize(3,4);
      std::stringstream value(child->FirstChild()->Value());
      for (double& v : sup.MVP) value >> v;
      if (!value)
      {
        std::cerr <<" *** SIMsupel::parse: Failed to read transformation matrix"
                  <<" for superelement \""<< sup.id <<"\"."<< std::endl;
        result = false;
      }
      IFEM::cout <<"  Parsing <"<< child->Value() <<">";
#ifdef SP_DEBUG
      IFEM::cout << sup.MVP;
#else
      IFEM::cout << std::endl;
#endif
    }
    else
      result &= this->SIMgeneric::parse(child);

  if (myModel.size() == last)
  {
    std::cerr <<" *** SIMsupel::parse: Missing superelement data for \""
              << sup.id <<"\",\""<< supNodeSet <<"\"."<< std::endl;
    result = false;
  }
  else if (myModel.size() > last+1)
  {
    std::cerr <<" *** SIMsupel::parse: Multiple ("<< myModel.size()-last
              <<") superelement data blocks detected for \""
              << sup.id <<"\",\""<< supNodeSet <<"\"."<< std::endl;
    result = false;
  }

  // Apply the superelement transformation and assign supernode set name
  if (result && !sup.MVP.empty())
    if ((result = myModel.back()->transform(sup.MVP)) && !supNodeSet.empty())
    {
      int topIdx = myModel.back()->parseNodeSet(supNodeSet,nullptr);
      if (topIdx > 0)
        myEntitys[supNodeSet].insert(TopItem(myModel.size(),topIdx,4));
    }

  if (result)
    mySups.push_back(sup);

  return result;
}


ASMbase* SIMsupel::readPatch (std::istream& isp, int pchInd, const CharVec&,
                              const char* whiteSpace) const
{
  ASMbase* pch = ASM3D::create(ASM::SuperElm,ncmp);
  if (pch)
  {
    if (!pch->read(isp) || this->getLocalPatchIndex(pchInd+1) < 1)
    {
      delete pch;
      pch = nullptr;
    }
    else
    {
      if (whiteSpace)
        IFEM::cout << whiteSpace <<"Reading patch "<< pchInd+1 << std::endl;
      pch->idx = myModel.size();
    }
  }

  return pch;
}


bool SIMsupel::createFEMmodel (char resetNumb)
{
  return this->SIMgeneric::createFEMmodel(resetNumb);
}


bool SIMsupel::recoverInternalDOFs (const Vector& glbSol)
{
  size_t pidx = 0;
  for (SuperElm& sup : mySups)
  {
    // Extract superelement solution vector from the global solution vector
    Vector supSol;
    int pchIdx = myModel[pidx]->idx + 1;
    myModel[pidx++]->extractNodalVec(glbSol, supSol, mySam->getMADOF());
#if SP_DEBUG > 2
    std::cout <<"\nSolution vector for superelement "<< pchIdx << supSol;
#endif

    if (sup.sim)
    {
      // Transform to local superelement axes
      bool ok = sup.MVP.empty() ? true : utl::transform(supSol,sup.MVP,true);

      // Recover the internal state of the superelement
      ok &= this->recoverInternalDOFs(myModel[pidx-1], sup, supSol);

      if (!sup.MVP.empty() && ok) // Transform back to global axes
        ok = utl::transform(sup.sol,sup.MVP);

      if (!ok)
      {
        std::cerr <<"\n *** SIMsupel::recoverInternalDOFs: Failed to"
                  <<" recover internal solution for superelement "<< pchIdx
                  << std::endl;
        return false;
      }
    }
    else // No substructure FE model - just use the superelement solution
      sup.sol = supSol;
  }

  return true;
}


bool SIMsupel::recoverInternalDOFs (const ASMbase*, SuperElm& sup,
                                    const Vector& supSol) const
{
  return sup.sim->recoverInternals(supSol,sup.sol);
}
