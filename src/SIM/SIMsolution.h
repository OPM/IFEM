// $Id$
//==============================================================================
//!
//! \file SIMsolution.h
//!
//! \date Nov 11 2017
//!
//! \author Knut Morten Okstad
//!
//! \brief General solution vector container for simulator drivers.
//!
//==============================================================================

#ifndef _SIM_SOLUTION_H_
#define _SIM_SOLUTION_H_

#include "MatVec.h"
#include <string>
#include <map>


/*!
  \brief Solution vector container with serialization support.
*/

class SIMsolution
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only
  SIMsolution() = default;
  //! \brief Empty default destructor.
  virtual ~SIMsolution() = default;

  //! \brief Initializes the solution vectors.
  //! \param[in] ndof Number of degrees of freedom
  //! \param[in] nsol Number of solution vectors
  void initSolution(size_t ndof, size_t nsol = 1);

  //! \brief Pushes the solution vector stack.
  //! \param[in] nVecState Number of vectors per state
  void pushSolution(unsigned short int nVecState = 1);

  using SerializeMap = std::map<std::string,std::string>; //!< Convenience type

  //! \brief Writes current solution to a serialization container.
  //! \param data Container for serialized data
  //! \param[in] name Name of simulator the solution belongs to
  bool saveSolution(SerializeMap& data, const std::string& name) const;

  //! \brief Restores the solution from a serialization container.
  //! \param[in] data Container for serialized data
  //! \param[in] name Name of simulator the solution belongs to
  bool restoreSolution(const SerializeMap& data, const std::string& name);

  //! \brief Helper method for serializing a double array into a text string.
  //! \param[in] v Pointer to double array
  //! \param[in] n Length of array
  //! \return The serialized string
  static std::string serialize(const double* v, size_t n);
  //! \brief Helper method for deserializing a double array from a text string.
  //! \param[in] data The serialized string
  //! \param[out] v Pointer to double array (assumed allocated outside)
  //! \param[in] n Length of array
  static void deSerialize(const std::string& data, double* v, size_t n);

public:
  //! \brief Returns a const reference to the solution vectors.
  virtual const Vectors& getSolutions() const { return solution; }
  //! \brief Returns a reference to the solution vectors (for assignment).
  virtual Vectors& theSolutions() { return solution; }

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int ix = 0) const { return solution[ix]; }
  //! \brief Modifies the current solution vector (used by sub-iterations only).
  virtual void setSolution(const RealArray& s, int ix = 0) { solution[ix] = s; }

protected:
  Vectors solution; //!< Stack of solution vectors
};

#endif
