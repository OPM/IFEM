// $Id$
//==============================================================================
//!
//! \file HDF5Base.C
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for interacting with HDF5 files.
//!
//==============================================================================

#include "HDF5Base.h"

#ifdef HAS_HDF5
#include <iostream>
#include <unistd.h>
#ifdef HAVE_MPI
#include "ProcessAdm.h"
#include <mpi.h>
#endif
#endif


HDF5Base::HDF5Base (const std::string& name, const ProcessAdm& adm)
  : m_hdf5_name(name)
#ifdef HAVE_MPI
  , m_adm(adm)
#endif
{
#ifdef HAS_HDF5
  m_file = -1;
#endif
}


bool HDF5Base::openFile (unsigned int flags)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    return true;

#ifdef HAVE_MPI
  hid_t acc_tpl;
  if (m_adm.dd.isPartitioned()) {
    if (m_adm.getProcId() != 0)
      return true;
    acc_tpl = H5P_DEFAULT;
  } else {
    MPI_Info info = MPI_INFO_NULL;
    acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(acc_tpl, *m_adm.getCommunicator(), info);
  }
#else
  hid_t acc_tpl = H5P_DEFAULT;
#endif

  if (flags == H5F_ACC_TRUNC)
    m_file = H5Fcreate(m_hdf5_name.c_str(),flags,H5P_DEFAULT,acc_tpl);
  else
    m_file = H5Fopen(m_hdf5_name.c_str(),flags,acc_tpl);

  if (m_file < 0)
  {
    std::cerr <<" *** HDF5Base: Failed to open "<< m_hdf5_name << std::endl;
    return false;
  }

#ifdef HAVE_MPI
  if (!m_adm.dd.isPartitioned())
    H5Pclose(acc_tpl);
#endif
  return true;
#else
  return false;
#endif
}


void HDF5Base::closeFile()
{
#ifdef HAS_HDF5
  if (m_file != -1)
    H5Fclose(m_file);
  m_file = -1;
#endif
}


#ifdef HAS_HDF5
bool HDF5Base::checkGroupExistence (hid_t parent, const char* path)
{
  // turn off errors to avoid cout spew
  H5E_BEGIN_TRY {
#if H5_VERS_MINOR > 8
    return H5Lexists(parent,path,H5P_DEFAULT) == 1;
#else
    return H5Gget_objinfo((hid_t)parent,path,0,nullptr) == 0;
#endif
  } H5E_END_TRY;
  return false;
}
#endif
