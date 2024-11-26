// $Id$
//==============================================================================
//!
//! \file HDF5Restart.C
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of restart data to HDF5 file.
//!
//==============================================================================

#include "HDF5Restart.h"
#include "TimeStep.h"
#include <iostream>
#include <sstream>

#ifdef HAS_HDF5
#include <numeric>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include "ProcessAdm.h"
#endif
#endif


HDF5Restart::HDF5Restart (const std::string& name, const ProcessAdm& adm,
                          int stride, int level0)
  : HDF5Base(name, adm), m_stride(stride), m_level(level0-1), m_last(0)
{
  if (m_hdf5_name.find('.') == std::string::npos)
    m_hdf5_name += ".hdf5";
}


bool HDF5Restart::dumpStep (const TimeStep& tp)
{
  if (tp.step < m_last + m_stride)
    return false;

  m_last = tp.step;
  return true;
}


bool HDF5Restart::writeData (const SerializeData& data)
{
#ifdef HAS_HDF5
  struct stat sb;
  int flag = stat(m_hdf5_name.c_str(),&sb) == 0 ? H5F_ACC_RDWR : H5F_ACC_TRUNC;
  if (!this->openFile(flag))
    return false;

  std::stringstream str;
  str << '/' << ++m_level;
  if (!checkGroupExistence(m_file,str.str().c_str()))
    H5Gclose(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));

#ifdef HAVE_MPI
  int pid  = m_adm.getProcId();
  int ptot = m_adm.getNoProcs();
#else
  int pid  = 0;
  int ptot = 1;
#endif

  // Lambda function for writing data to the HDF5 file.
#ifdef HAVE_MPI
  auto&& writeGroup = [this,pid,ptot]
#else
  auto&& writeGroup = []
#endif
      (hid_t group, const std::string& name,
       int len, const void* data, hid_t type)
  {
#ifdef HAVE_MPI
    std::vector<int> lens(ptot), lens2(ptot);
    std::fill(lens.begin(), lens.end(),len);
    MPI_Alltoall(lens.data(),1,MPI_INT,lens2.data(),1,MPI_INT,*m_adm.getCommunicator());
    hsize_t siz   = std::accumulate(lens2.begin(), lens2.end(), hsize_t(0));
    hsize_t start = std::accumulate(lens2.begin(), lens2.begin() + pid, hsize_t(0));
#else
    hsize_t siz   = static_cast<hsize_t>(len);
    hsize_t start = 0;
#endif
    hid_t space = H5Screate_simple(1,&siz,nullptr);
    hid_t set = H5Dcreate2(group,name.c_str(),
                           type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if (len > 0) {
      hid_t filespace = H5Dget_space(set);
      siz = len;
      hsize_t stride = 1;
      H5Sselect_hyperslab(filespace,H5S_SELECT_SET,&start,&stride,&siz,nullptr);
      hid_t memspace = H5Screate_simple(1,&siz,nullptr);
      H5Dwrite(set,type,memspace,filespace,H5P_DEFAULT,data);
      H5Sclose(memspace);
      H5Sclose(filespace);
    }
    H5Dclose(set);
    H5Sclose(space);
  };

  for (int p = 0; p < ptot; p++)
    for (const std::pair<const std::string,std::string>& it : data) {
      std::stringstream str;
      str << m_level << '/' << p;
      hid_t group;
      if (checkGroupExistence(m_file,str.str().c_str()))
        group = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
      else
        group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
      if (!H5Lexists(group, it.first.c_str(), 0)) {
        int len = pid == p ? it.second.size() : 0;
        writeGroup(group, it.first, len, it.second.data(), H5T_NATIVE_CHAR);
      }
      H5Gclose(group);
    }

  this->closeFile();

#else
  std::cout <<"HDF5Restart: Compiled without HDF5 support, no data written."<< std::endl;
#endif

  return true;
}


//! \brief Convenience type for restart IO.
typedef std::pair<HDF5Restart::SerializeData&,bool> read_restart_ctx;


#ifdef HAS_HDF5
//! \brief Static helper used for reading data.
static herr_t read_restart_data (hid_t group_id, const char* name, void* data)
{
  bool readBasis = strstr(name,"::basis") != nullptr;
  read_restart_ctx* ctx = static_cast<read_restart_ctx*>(data);
  if (ctx->second != readBasis)
    return 0; // Not reading the right data type (basis or solution)

  hid_t setId = H5Dopen2(group_id,name,H5P_DEFAULT);
  if (setId < 0) return setId; // Result set not found
  hsize_t siz = H5Dget_storage_size(setId);
  if (siz < 1) return 0; // Empty result set
  char* cdata = new char[siz];
  H5Dread(setId,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,cdata);
  H5Dclose(setId);
  ctx->first[name] = std::string(cdata,siz);

  return 0;
}
#endif


int HDF5Restart::readData (SerializeData& data, int level, bool basis)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return -1;

  if (level < 0)
  {
    // Lambda function extracting the highest time level among the groups.
    auto&& find_last_level = [] (hid_t, const char* name, void* data) -> herr_t
    {
      int ilevel = atoi(name);
      int* level = static_cast<int*>(data);
      if (ilevel > *level) *level = ilevel;
      return 0;
    };

    int idx = 0;
    H5Giterate(m_file, "/", &idx, find_last_level, &level);
  }

  std::stringstream str;
  str << '/' << level << '/';
#ifdef HAVE_MPI
  str << m_adm.getProcId();
#else
  str << 0;
#endif
  int idx = 0;
  read_restart_ctx ctx(data,basis);
  int it = H5Giterate(m_file, str.str().c_str(), &idx, read_restart_data, &ctx);
  return it < 0 ? it : level;
#else
  std::cout <<"HDF5Writer: Compiled without HDF5 support, no data read."<< std::endl;
  return -1;
#endif
}
