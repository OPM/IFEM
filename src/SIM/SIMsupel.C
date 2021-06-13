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
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


/*!
  \brief Empty integrand class for superelement analysis.
*/

class NoProblem : public IntegrandBase
{
public:
  //! \brief The constructor forwards to the parent class constructor,
  explicit NoProblem(unsigned char n) : IntegrandBase(n) {}
  //! \brief Empty destructor,
  virtual ~NoProblem() {}
};


SIMsupel::SIMsupel (const char* hd, unsigned char nf) : ncmp(nf)
{
  opt.pSolOnly = true;
  myProblem = new NoProblem(nf);
  if (hd) myHeading = hd;
}


bool SIMsupel::parse (const TiXmlElement* elem)
{
  if (strncasecmp(elem->Value(),"superel",7))
    return this->SIMgeneric::parse(elem);

  bool result = true;
  size_t last = myModel.size();
  SuperElm sup;
  std::string supNodeSet;
  utl::getAttribute(elem,"id",sup.id);
  utl::getAttribute(elem,"nodeset",supNodeSet);
  const TiXmlElement* child = elem->FirstChildElement();
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
      int topIdx = 0;
      myModel.back()->getNodeSet(supNodeSet,topIdx);
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
