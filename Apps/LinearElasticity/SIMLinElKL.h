// $Id$
//==============================================================================
//!
//! \file SIMLinElKL.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love plates.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_KL_H
#define _SIM_LIN_EL_KL_H

#include "SIMLinEl2D.h"


/*!
  \brief Driver class for isogeometric FEM analysis of Kirchhoff-Love plates.
*/

class SIMLinElKL : public SIMLinEl2D
{
public:
  //! \brief Default constructor.
  SIMLinElKL();
  //! \brief Empty destructor.
  virtual ~SIMLinElKL() {}

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);

private:
  RealArray tVec; //!< Plate thickness data
};

#endif
