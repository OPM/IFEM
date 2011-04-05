// $Id$
//==============================================================================
//!
//! \file LinAlgInit.h
//!
//! \date May 06 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Initializer for linear algebra packages.
//!
//==============================================================================

#ifndef _LINALG_INIT_H
#define _LINALG_INIT_H


/*!
  \brief Initializer for linear algebra packages.

  \details This class is a reference-counted singleton. We need control over the
  destruction order since the destructor calls \a PetscFinalize() which shuts
  down MPI. Thus, the LinAlgInit object must be destroyed after all other
  objects that do MPI calls in their destructor.
  This cannot be guaranteed without the reference-counting since this class is
  instanciated as a local variable in the \a main() scope of the applications.
*/

class LinAlgInit
{
  //! \brief The constructor uses the command-line arguments to set up things.
  LinAlgInit(int argc, char** argv);
public:
  //! \brief The destructor finalizes the linear algebra packages.
  ~LinAlgInit();

  //! \brief Increments the reference counter.
  static void increfs() { ++refs; }
  //! \brief Decrements the reference counter.
  static void decrefs() { if (--refs == 0) delete instance; }

  //! \brief Instanciates the singleton LinAlgInit object.
  static LinAlgInit& Init(int argc, char** argv);

  int myPid; //!< Processor ID in parallel simulations

private:
  static int refs; //!< Reference counter

  static LinAlgInit* instance; //!< Pointer to the singleton object
};

#endif
