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

#include <iostream>
#include <vector>

class RealFunc;
class VecFunc;
class TensorFunc;
class STensorFunc;
class TiXmlElement;


/*!
  \brief Class for analytical solution fields (primary and secondary solution).
*/

class AnaSol
{
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
  explicit AnaSol(RealFunc* s1 = nullptr, VecFunc* s2 = nullptr,
                  VecFunc* v1 = nullptr, TensorFunc* v2 = nullptr,
                  STensorFunc* v3 = nullptr);

  //! \brief Constructor initializing the primary and secondary solution fields.
  //! \param[in] s Primary scalar solution field
  //! \param[in] sigma Symmetric stress tensor field
  AnaSol(RealFunc* s, STensorFunc* sigma);

  //! \brief Constructor initializing the primary and secondary solution fields.
  //! \param[in] v Primary vector solution field
  //! \param[in] sigma Symmetric stress tensor field
  AnaSol(VecFunc* v, STensorFunc* sigma)
    : vecSol(v), vecSecSol(nullptr), stressSol(sigma) {}

  //! \brief Constructor initializing the symmetric stress tensor field only.
  //! \param[in] sigma Symmetric stress tensor field
  explicit AnaSol(STensorFunc* sigma)
    : vecSol(nullptr), vecSecSol(nullptr), stressSol(sigma) {}

  //! \brief Constructor initializing expression functions by parsing a stream.
  //! \param is The file stream to read function definition from
  //! \param[in] nlines Number of lines to read
  //! \param[in] scalarSol If \e true, the primary solution is a scalar field
  AnaSol(std::istream& is, const int nlines, bool scalarSol = true);

  //! \brief Constructor initializing expression functions by parsing XML tags.
  //! \param[in] elem Pointer to XML-element to extract data from
  //! \param[in] scalarSol If \e true, the primary solution is a scalar field
  explicit AnaSol(const TiXmlElement* elem, bool scalarSol = true);

  //! \brief No copying of this class.
  AnaSol(const AnaSol&) = delete;

  //! \brief The destructor frees the analytical solution fields.
  virtual ~AnaSol();

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

  //! \brief Sets the patch to use.
  void initPatch(size_t pIdx);

  //! \brief Make sure we have a secondary solution.
  //! \details If none is given, we use derivation (automatic or finite difference)
  //!          to obtain one.
  virtual void setupSecondarySolutions();

private:
  //! \brief Parses expression functions from XML definition.
  template<class Scalar>
  void parseExpressionFunctions(const TiXmlElement* elem, bool scalarSol);

  //! \brief Parses field functions from XML definition.
  void parseFieldFunctions(const TiXmlElement* elem, bool scalarSol);

protected:
  bool symmetric = false; //!< True to use symmetric secondary solution

  std::vector<RealFunc*>   scalSol; //!< Primary scalar solution fields
  std::vector<VecFunc*> scalSecSol; //!< Secondary scalar solution fields

  VecFunc*     vecSol;    //!< Primary vector solution field
  TensorFunc*  vecSecSol; //!< Secondary solution field (vector gradient field)
  STensorFunc* stressSol; //!< Secondary solution field (stress field)

};

#endif
