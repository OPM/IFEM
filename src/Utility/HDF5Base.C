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


bool HDF5Base::openFile (unsigned int flags, bool latest)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    return true;

  hid_t acc_tpl = H5P_DEFAULT;
  bool acc_tpl_alloc = false;
  if (latest
#ifdef HAVE_MPI
      && m_adm.getProcId() == 0 || !m_adm.dd.isPartitioned()
 #endif
      ) {
    acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
    acc_tpl_alloc = true;
  }

#ifdef HAVE_MPI
  if (m_adm.dd.isPartitioned()) {
    if (m_adm.getProcId() != 0)
      return true;
  } else {
    MPI_Info info = MPI_INFO_NULL;
    if (!acc_tpl_alloc) {
        acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
        acc_tpl_alloc = true;
    }
    H5Pset_fapl_mpio(acc_tpl, *m_adm.getCommunicator(), info);
  }
#endif

  if (latest)
    H5Pset_libver_bounds(acc_tpl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

  if (flags == H5F_ACC_TRUNC)
    m_file = H5Fcreate(m_hdf5_name.c_str(),flags,H5P_DEFAULT,acc_tpl);
  else
    m_file = H5Fopen(m_hdf5_name.c_str(),flags,acc_tpl);

  if (m_file < 0)
  {
    std::cerr <<" *** HDF5Base: Failed to open "<< m_hdf5_name << std::endl;
    if (acc_tpl_alloc)
      H5Pclose(acc_tpl);
    return false;
  }

  if (acc_tpl_alloc)
    H5Pclose(acc_tpl);

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
