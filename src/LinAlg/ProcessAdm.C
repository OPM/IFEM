//==============================================================================
//!
//! \file ProcessAdm.C
//!
//! \date Sep 30, 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for administration MPI processes used by the IFEM library,
//! in particular parallel linear algebra using PETSc.
//!
//==============================================================================

#include "IFEM.h"
#include "ProcessAdm.h"
#include "LinAlgInit.h"
#include <iostream>

ProcessAdm::ProcessAdm() : cout(std::cout)
{
  LinAlgInit::increfs();
#ifdef HAS_PETSC
  MPI_Comm_dup(PETSC_COMM_SELF,&comm);
#endif
  myPid = 0;
  nProc = 1;
  parallel = false;
}

#if defined(HAS_PETSC) || defined(HAVE_MPI)
ProcessAdm::ProcessAdm(MPI_Comm& mpi_comm) : cout(std::cout)
{
  LinAlgInit::increfs();
#ifdef HAVE_MPI
  MPI_Comm_dup(mpi_comm,&comm);
  MPI_Comm_rank(comm,&myPid);
  MPI_Comm_size(comm,&nProc);
  cout = IFEM::cout;
  parallel = nProc > 1;
#else
  MPI_Comm_dup(PETSC_COMM_SELF,&comm);
  myPid = 0;
  nProc = 1;
  parallel = false;
#endif
}
#endif


#ifdef HAVE_MPI
ProcessAdm::ProcessAdm(bool) : cout(std::cout)
{
  LinAlgInit::increfs();
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  MPI_Comm_rank(comm,&myPid);
  MPI_Comm_size(comm,&nProc);
  cout = IFEM::cout;
  parallel = nProc > 1;
}
#endif


ProcessAdm::~ProcessAdm()
{
  myPid = nProc = 0;
#ifdef HAS_PETSC
  if (parallel)
    MPI_Comm_free(&comm);
#endif
  parallel = false;
  LinAlgInit::decrefs();
}


#if defined(HAS_PETSC) || defined(HAVE_MPI)
void ProcessAdm::setCommunicator(const MPI_Comm* comm2)
{
  if (parallel)
    MPI_Comm_free(&comm);
  MPI_Comm_dup(*comm2, &comm);
  MPI_Comm_rank(comm,&myPid);
  MPI_Comm_size(comm,&nProc);
  parallel = nProc > 1;
}


ProcessAdm& ProcessAdm::operator=(const ProcessAdm& adm2)
{
#if defined(HAS_PETSC) || defined(HAVE_MPI)
  MPI_Comm_dup(adm2.comm,&comm);
#endif
  myPid = adm2.myPid;
  nProc = adm2.nProc;
  parallel = adm2.parallel;
  cout = adm2.cout;

  return *this;
}


void ProcessAdm::send(int value, int dest) const
{
#ifdef HAVE_MPI
  if ((dest >= 0) && (dest < nProc))
    MPI_Send(&value,1,MPI_INT,dest,101,comm);
#endif
}


void ProcessAdm::send(std::vector<int>& ivec, int dest) const
{
#ifdef HAVE_MPI
  int n = ivec.size();
  if ((dest >= 0) && (dest < nProc))
    MPI_Send(&(ivec[0]),n,MPI_INT,dest,102,comm);
#endif
}


void ProcessAdm::send(double value, int dest) const
{
#ifdef HAVE_MPI
  if ((dest >= 0) && (dest < nProc))
    MPI_Send(&value,1,MPI_DOUBLE,dest,103,comm);
#endif
}


void ProcessAdm::send(std::vector<double>& rvec, int dest) const
{
#ifdef HAVE_MPI
  int n = rvec.size();
  if ((dest >= 0) && (dest < nProc))
    MPI_Send(&(rvec[0]),n,MPI_DOUBLE,dest,104,comm);
#endif
}


void ProcessAdm::receive(int& value, int source) const
{
#ifdef HAVE_MPI
  value = 0;
  MPI_Status status;
  if ((source >= 0) && (source < nProc))
    MPI_Recv(&value,1,MPI_INT,source,101,comm,&status);
#endif
}


void ProcessAdm::receive(std::vector<int>& ivec, int source) const
{
#ifdef HAVE_MPI
  int n = ivec.size();
  MPI_Status status;
  if ((source >= 0) && (source < nProc))
    MPI_Recv(&(ivec[0]),n,MPI_INT,source,102,comm,&status);
#endif
}


void ProcessAdm::receive(double& value, int source) const
{
#ifdef HAVE_MPI
  value = 0.0;
  MPI_Status status;
  if ((source >= 0) && (source < nProc))
    MPI_Recv(&value,1,MPI_DOUBLE,source,103,comm,&status);
#endif
}


void ProcessAdm::receive(std::vector<double>& rvec, int source) const
{
#ifdef HAVE_MPI
  int n = rvec.size();
  MPI_Status status;
  if ((source >= 0) && (source < nProc))
    MPI_Recv(&(rvec[0]),n,MPI_DOUBLE,source,104,comm,&status);
#endif
}


int ProcessAdm::allReduce(int value, MPI_Op oper) const
{
  int tmp;
#ifdef HAVE_MPI
  MPI_Allreduce(&value,&tmp,1,MPI_INT,oper,comm);
#else
  tmp = value;
#endif
  return tmp;
}


void ProcessAdm::allReduce(std::vector<int>& ivec, MPI_Op oper) const
{
#ifdef HAVE_MPI
  int n = ivec.size();
  std::vector<int> tmp;
  tmp.resize(n);
  MPI_Allreduce(&(ivec[0]),&(tmp[0]),n,MPI_INT,oper,comm);
  ivec = tmp;
#endif
}

double ProcessAdm::allReduce(double value, MPI_Op oper) const
{
  double tmp;
#ifdef HAVE_MPI
  MPI_Allreduce(&value,&tmp,1,MPI_DOUBLE,oper,comm);
#else
  tmp = value;
#endif
  return tmp;
}


void ProcessAdm::allReduce(std::vector<double>& rvec, MPI_Op oper) const
{
#ifdef HAVE_MPI
  int n = rvec.size();
  std::vector<double> tmp;
  tmp.resize(n);
  MPI_Allreduce(&(rvec[0]),&(tmp[0]),n,MPI_DOUBLE,oper,comm);
  rvec = tmp;
#endif
}
#endif
