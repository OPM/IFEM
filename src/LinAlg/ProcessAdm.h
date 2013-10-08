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
#include "mpi.h"
#include <vector>
#include "petscsys.h"
#endif

/*!
  \brief Class for administration of MPI processes in IFEM library.
*/

class ProcessAdm
{
 protected:
  int myPid;       //!< Process id
  int nProc;       //!< Number of processes in communicator
  bool parallel;   //!< If Processor is parallel

#ifdef HAS_PETSC
  MPI_Comm comm;   //!< MPI communicator
#endif

 public:
  //! \brief Construct an empty (serial) process administrator.
  ProcessAdm();
#ifdef HAS_PETSC
  //! \brief Construct a parallel process administrator.
  ProcessAdm(MPI_Comm& mpi_comm);
#endif
  //! \brief Destructor
  virtual ~ProcessAdm();

  //! \brief Return process id
  int getProcId() const  { return myPid; }
  //! \brief Return rank
  int getNoProcs() const { return nProc; }
  //! \brief Return if parallel
  int isParallel() const { return parallel; }

#ifdef HAS_PETSC
  //! \brief Return MPI communicator
  const MPI_Comm* getCommunicator() const;

  //! \brief Blocking send of an integer
  //! \param[in] value Integer to be sent
  //! \param[in] dest  Process id for destination
  virtual void send(int value, int dest) const;
  //! \brief Blocking send of an integer vector
  //! \param[in] ivec  Vector to be sent
  //! \param[in] dest  Process id for destination
  virtual void send(std::vector<int>& ivec, int dest) const;
  //! \brief Blocking send of a double
  //! \param[in] value Double to be sent
  //! \param[in] dest  Process id for destination
  virtual void send(double value, int dest) const;
  //! \brief Blocking send of a double vector
  //! \param[in] rvec  Vector to be sent
  //! \param[in] dest  Process id for destination
  virtual void send(std::vector<double>& rvec, int dest) const;

  //! \brief Blocking receive of an integer
  //! \param[out] value  Received value
  //! \param[in] source Process id for source
  virtual void receive(int& value, int source) const;
  //! \brief Blocking receive of an integer vector
  //! \param[out] ivec   Integer vector to receive
  //! \param[in]  source Process id for source
  virtual void receive(std::vector<int>& ivec, int source) const;
  //! \brief Blocking receive of a double
  //! \param[out] value  Received value
  //! \param[in]  source Process id for source
  virtual void receive(double& value, int source) const;
  //! \brief Blocking receive of a double vector
  //! \param[out] rvec   Double vector to receive
  //! \param[in]  source Process id for source
  virtual void receive(std::vector<double>& rvec, int source) const;

  //! \brief AllReduce for integer
  //! \param[in] value Integer to be reduced
  //! \param[in] oper  MPI operator (ex. MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, ..)
  virtual int allReduce(int value, MPI_Op oper) const;
  //! \brief AllReduce for integer vector
  //! \param[inout] ivec  Integer array to be reduced
  //! \param[in]    oper  MPI operator (ex. MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, ..)
  virtual void allReduce(std::vector<int>& ivec, MPI_Op oper) const;
  //! \brief AllReduce for double
  //! \param[in] value Double to be reduced
  //! \param[in] oper  MPI operator (ex. MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, ..)
  virtual double allReduce(double value, MPI_Op oper) const;
  //! \brief AllReduce for double vector
  //! \param[inout] rvec  Double array to be reduced
  //! \param[in]    oper  MPI operator (ex. MPI_MIN, MPI_MAX, MPI_SUM, MPI_PROD, ..)
  virtual void allReduce(std::vector<double>& ivec, MPI_Op oper) const;
#endif
};

#endif
