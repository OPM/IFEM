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
public:
  //! \brief Default constructor.
  SIMsolution() {}
  //! \brief Empty destructor.
  virtual ~SIMsolution() {}

  //! \brief Initializes the solution vectors.
  //! \param[in] ndof Number of degrees of freedom
  //! \param[in] nsol Number of solution vectors
  bool initSolution(size_t ndof, size_t nsol = 1);

protected:
  //! \brief Pushes the solution vector stack.
  //! \param[in] nsol Number of vectors to push (0 = all)
  void pushSolution(size_t nsol = 0);

  typedef std::map<std::string,std::string> SerializeMap; //!< Convenience type

  //! \brief Writes current solution to a serialization container.
  //! \param data Container for serialized data
  //! \param[in] name Name of simulator the solution belongs to
  bool saveSolution(SerializeMap& data, const std::string& name) const;

  //! \brief Restores the solution from a serialization container.
  //! \param[in] data Container for serialized data
  //! \param[in] name Name of simulator the solution belongs to
  bool restoreSolution(const SerializeMap& data, const std::string& name);

public:
  //! \brief Returns a const reference to the solution vectors.
  virtual const Vectors& getSolutions() const { return solution; }
  //! \brief Returns a reference to the solution vectors (for assignment).
  virtual Vectors& theSolutions() { return solution; }

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int idx = 0) const { return solution[idx]; }
  //! \brief Modifies the current solution vector (used by sub-iterations only).
  virtual void setSolution(const Vector& s, int idx = 0) { solution[idx] = s; }

protected:
  Vectors solution; //!< Stack of solution vectors
};

#endif
