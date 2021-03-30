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
#include "IntegrandBase.h"
#include "ASMbase.h"
#include "ASM3D.h"
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
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    result &= this->SIMgeneric::parse(child);

  return true;
}


ASMbase* SIMsupel::readPatch (std::istream& isp, int pchInd, const CharVec&,
                              const char* whiteSpace) const
{
  ASMbase* pch = ASM3D::create(ASM::SuperElm,{ncmp});
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
