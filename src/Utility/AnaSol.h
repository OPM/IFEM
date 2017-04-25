// $Id$
//==============================================================================
//!
//! \file AnaSol.h
//!
//! \date Dec 7 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Analytical solution fields (primary and secondary).
//!
//==============================================================================

#ifndef _ANA_SOL_H
#define _ANA_SOL_H

#include "Function.h"
#include <iostream>
#include <vector>

class TiXmlElement;


/*!
  \brief Class for analytical solution fields (primary and secondary solution).
*/

class AnaSol
{
protected:
  std::vector<RealFunc*>    scalSol; //!< Primary scalar solution fields
  std::vector<VecFunc*>  scalSecSol; //!< Secondary scalar solution fields (gradient fields)

  VecFunc*     vecSol;     //!< Primary vector solution field
  TensorFunc*  vecSecSol;  //!< Secondary solution field (vector gradient field)
  STensorFunc* stressSol;  //!< Secondary solution field (stress field)

public:
  //! \brief Default constructor initializing all solution fields.
  //! \param[in] s1 Primary scalar solution field
  //! \param[in] s2 Secondary solution field, gradient
  //! \param[in] v1 Primary vector solution field
  //! \param[in] v2 Secondary solution field, gradient
  //! \param[in] v3 Secondary solution field, symmetric stress tensor
  //!
  //! \details It is assumed that all the arguments are pointers to dynamically
  //! allocated objects, as the class destructor will attempt to delete them.
  AnaSol(RealFunc* s1 = 0, VecFunc* s2 = 0,
         VecFunc* v1 = 0, TensorFunc* v2 = 0, STensorFunc* v3 = 0);

  //! \brief Constructor initializing the primary and secondary solution fields.
  //! \param[in] s Primary scalar solution field
  //! \param[in] sigma Symmetric stress tensor field
  AnaSol(RealFunc* s, STensorFunc* sigma)
    : vecSol(0), vecSecSol(0), stressSol(sigma) { if (s) scalSol.push_back(s); }

  //! \brief Constructor initializing the primary and secondary solution fields.
  //! \param[in] v Primary vector solution field
  //! \param[in] sigma Symmetric stress tensor field
  AnaSol(VecFunc* v, STensorFunc* sigma)
    : vecSol(v), vecSecSol(0), stressSol(sigma) {}

  //! \brief Constructor initializing the symmetric stress tensor field only.
  //! \param[in] sigma Symmetric stress tensor field
  AnaSol(STensorFunc* sigma)
    : vecSol(0), vecSecSol(0), stressSol(sigma) {}

  //! \brief Constructor initializing expression functions by parsing a stream.
  //! \param is The file stream to read function definition from
  //! \param[in] nlines Number of lines to read
  //! \param[in] scalarSol If \e true, the primary solution is a scalar field
  AnaSol(std::istream& is, const int nlines, bool scalarSol = true);
  //! \brief Constructor initializing expression functions by parsing XML tags.
  //! \param[in] elem Pointer to XML-element to extract data from
  //! \param[in] scalarSol If \e true, the primary solution is a scalar field
  AnaSol(const TiXmlElement* elem, bool scalarSol = true);

  //! \brief The destructor frees the analytical solution fields.
  virtual ~AnaSol()
  {
    for (auto& it : scalSol)
      delete it;
    for (auto& it : scalSecSol)
      delete it;
    delete vecSol;
    delete vecSecSol;
    delete stressSol;
  }

  //! \brief Checks whether a scalar solution is defined.
  char hasScalarSol() const
  {
    if (stressSol && !vecSecSol && !vecSol)
      return 3;
    else if (!scalSecSol.empty())
      return 2;
    else if (!scalSol.empty())
      return 1;
    else
      return 0;
  }

  //! \brief Checks whether a vector solution is defined.
  char hasVectorSol() const
  {
    if (stressSol)
      return 3;
    else if (vecSecSol)
      return 2;
    else if (vecSol)
      return 1;
    else
      return 0;
  }

  //! \brief Returns the scalar solution, if any.
  RealFunc* getScalarSol(size_t idx = 0) const
  { return scalSol.size() <= idx ? nullptr : scalSol[idx]; }
  //! \brief Returns the secondary scalar solution, if any.
  VecFunc* getScalarSecSol(size_t idx = 0) const
  { return scalSecSol.size() <= idx ? nullptr : scalSecSol[idx]; }

  //! \brief Returns the vector solution, if any.
  VecFunc* getVectorSol() const { return vecSol; }
  //! \brief Returns the secondary vector solution, if any.
  TensorFunc* getVectorSecSol() const { return vecSecSol; }

  //! \brief Returns the stress solution, if any.
  STensorFunc* getStressSol() const { return stressSol; }
};

#endif
