//==============================================================================
//!
//! \file ProcessAdm.h
//!
//! \date Sep 30, 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for administration MPI processes used by the IFEM library,
//! in particular parallel linear algebra using PETSc.
//!
//==============================================================================

#ifndef _PROCESS_ADM_H
#define _PROCESS_ADM_H

#ifdef HAS_PETSC
#include "petscsys.h"
#elif defined(HAVE_MPI)
#include "mpi.h"
#endif
#include <vector>

#include "LogStream.h"


/*!
  \brief Class for administration of MPI processes in IFEM library.
*/

class ProcessAdm
{
  int myPid;       //!< Process id
  int nProc;       //!< Number of processes in communicator
  bool parallel;   //!< If Processor is parallel

#if defined(HAS_PETSC) || defined(HAVE_MPI)
  MPI_Comm comm;   //!< MPI communicator
#endif

public:
  mutable utl::LogStream cout; //!< Combined standard out for this process group

  //! \brief Construct an empty (serial) process administrator.
  ProcessAdm();
#if defined(HAS_PETSC) || defined(HAVE_MPI)
  //! \brief Construct a parallel process administrator.
  ProcessAdm(MPI_Comm& mpi_comm);
#endif

#ifdef HAVE_MPI
  //! \brief Construct a parallel process administrator.
  //! \details This overload is necessary due to MPI_COMM_WORLD being .. ickily.
  ProcessAdm(bool hack);
#endif

  //! \brief The destructor releases the process administrator.
  ~ProcessAdm();

  //! \brief Return process id
  int getProcId() const  { return myPid; }
  //! \brief Return rank
  int getNoProcs() const { return nProc; }
  //! \brief Return if parallel
  int isParallel() const { return parallel; }

#if defined(HAS_PETSC) || defined(HAVE_MPI)
  //! \brief Return MPI communicator
  MPI_Comm* getCommunicator() { return &comm; }
  //! \brief Return MPI communicator
  const MPI_Comm* getCommunicator() const { return &comm; }

  //! \brief Set MPI communicator
  //! \param comm New communicator
  void setCommunicator(const MPI_Comm* comm2);

  //! \brief Assignment operator
  ProcessAdm& operator=(const ProcessAdm&);

  //! \brief Blocking send of an integer
  //! \param[in] value Integer to be sent
  //! \param[in] dest  Process id for destination
  void send(int value, int dest) const;
  //! \brief Blocking send of an integer vector
  //! \param[in] ivec  Vector to be sent
  //! \param[in] dest  Process id for destination
  void send(std::vector<int>& ivec, int dest) const;
  //! \brief Blocking send of a double
  //! \param[in] value Double to be sent
  //! \param[in] dest  Process id for destination
  void send(double value, int dest) const;
  //! \brief Blocking send of a double vector
  //! \param[in] rvec  Vector to be sent
  //! \param[in] dest  Process id for destination
  void send(std::vector<double>& rvec, int dest) const;

  //! \brief Blocking receive of an integer
  //! \param[out] value  Received value
  //! \param[in] source Process id for source
  void receive(int& value, int source) const;
  //! \brief Blocking receive of an integer vector
  //! \param[out] ivec   Integer vector to receive
  //! \param[in]  source Process id for source
  void receive(std::vector<int>& ivec, int source) const;
  //! \brief Blocking receive of a double
  //! \param[out] value  Received value
  //! \param[in]  source Process id for source
  virtual void receive(double& value, int source) const;
  //! \brief Blocking receive of a double vector
  //! \param[out] rvec   Double vector to receive
  //! \param[in]  source Process id for source
  void receive(std::vector<double>& rvec, int source) const;

  //! \brief AllReduce for integer value.
  //! \param[in] value Integer to be reduced
  //! \param[in] oper MPI operator (MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, etc.)
  int allReduce(int value, MPI_Op oper) const;
  //! \brief AllReduce for integer vector.
  //! \param vec Integer array to be reduced
  //! \param[in] oper MPI operator (MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, etc.)
  void allReduce(std::vector<int>& vec, MPI_Op oper) const;

  //! \brief AllReduce for double value.
  //! \param[in] value Double to be reduced
  //! \param[in] oper MPI operator (MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, etc.)
  double allReduce(double value, MPI_Op oper) const;
  //! \brief AllReduce for double vector.
  //! \param vec Double array to be reduced
  //! \param[in] oper MPI operator (MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, etc.)
  void allReduce(std::vector<double>& vec, MPI_Op oper) const;
#endif

  //! \brief AllReduce with MPI_SUM for a vector.
  template<class T> void allReduceAsSum(std::vector<T>& vec) const
  {
#ifdef HAS_PETSC
    if (parallel)
      this->allReduce(vec,MPI_SUM);
#endif
  }
};

#endif
