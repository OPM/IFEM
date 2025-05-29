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
  of a model consisting of superelements resulting from static condensation
  or reduced order modeling.
*/

class SIMsupel : public SIMdummy<SIMgeneric>
{
protected:
  //! \brief Struct with superelement data.
  struct SuperElm
  {
    std::string id;  //!< Superelement id
    Matrix      MVP; //!< Local-to-global transformation matrix
    Vector      sol; //!< Recovered primary solution on underlying FE model

    SIMgeneric* sim = nullptr; //!< Pointer to the underlying FE model
  };

public:
  //! \brief Default constructor.
  //! \param[in] hd Sub-simulator heading
  //! \param[in] nf Dimension of the primary solution field
  explicit SIMsupel(const char* hd = nullptr, char nf = 6);

  //! \brief Creates the FE topology of the superelement patches.
  virtual bool createFEMmodel(char resetNumb);

  //! \brief Returns the name of this simulator.
  virtual std::string getName() const { return "SIMsupel"; }

  //! \brief Performs recovery of the internal DOFs for superelements.
  //! \param[in] glbSol Global solution vector
  bool recoverInternalDOFs(const Vector& glbSol);

protected:
  using SIMdummy<SIMgeneric>::parse;
  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] whiteSpace For message formatting
  virtual ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec&,
                             const char* whiteSpace) const;

  //! \brief Recovers the internal DOFs from static condensation.
  //! \param sup The superelement to recover internal DOFs for
  //! \param[in] supSol Supernode solution values
  virtual bool recoverInternalDOFs(const ASMbase*, SuperElm& sup,
                                   const Vector& supSol) const;

private:
  unsigned char ncmp; //!< Number of primary solution components per node

protected:
  std::vector<SuperElm> mySups; //!< Superelement data container
};

#endif
