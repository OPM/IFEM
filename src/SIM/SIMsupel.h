// $Id$
//==============================================================================
//!
//! \file SIMsupel.h
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for general superelement FEM analysis.
//!
//==============================================================================

#ifndef _SIM_SUPEL_H
#define _SIM_SUPEL_H

#include "SIMdummy.h"
#include "SIMgeneric.h"


/*!
  \brief Solution driver for general superelement FEM analysis.
  \details This class contains methods for setting up a FE analysis
  of a model consisting of superelements resulting from reduced order modeling.
*/

class SIMsupel : public SIMdummy<SIMgeneric>
{
public:
  //! \brief Default constructor.
  //! \param[in] hd Sub-simulator heading
  //! \param[in] nf Dimension of the primary solution field
  explicit SIMsupel(const char* hd = nullptr, unsigned char nf = 6);
  //! \brief Empty destructor.
  virtual ~SIMsupel() {}

  //! \brief Creates the computational FEM model from the spline patches.
  //! \param[in] resetNumb If \e 'y', start element and node numbers from zero
  virtual bool createFEMmodel(char resetNumb);

  //! \brief Returns the name of this simulator.
  virtual std::string getName() const { return "SIMsupel"; }

protected:
  using SIMdummy<SIMgeneric>::parse;
  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] whiteSpace For message formatting
  virtual ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec&,
                             const char* whiteSpace) const;

private:
  unsigned char ncmp; //!< Number of primary solution components per node
};

#endif
