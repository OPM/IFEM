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
#ifdef HAVE_MPI
#include <mpi.h>
#include "ProcessAdm.h"
#endif
#endif


HDF5Restart::HDF5Restart (const std::string& name, const ProcessAdm& adm,
                          int stride)
  : HDF5Base(name, adm), m_stride(stride)
{
}


bool HDF5Restart::dumpStep(const TimeStep& tp)
{
  return (tp.step % m_stride) == 0;
}


void HDF5Restart::readArray(hid_t group, const std::string& name,
                            int& len,  char*& data)
{
#ifdef HAS_HDF5
  hid_t set = H5Dopen2(group,name.c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set);
  len = siz;
  data = new char[siz];
  H5Dread(set,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(set);
#else
  len = 0;
  std::cout << "HDF5Restart: compiled without HDF5 support, no data read" << std::endl;
#endif
}


void HDF5Restart::writeArray(hid_t group, const std::string& name,
                             int len, const void* data, hid_t type)
{
#ifdef HAS_HDF5
#ifdef HAVE_MPI
  int lens[m_adm.getNoProcs()], lens2[m_adm.getNoProcs()];
  std::fill(lens,lens+m_adm.getNoProcs(),len);
  MPI_Alltoall(lens,1,MPI_INT,lens2,1,MPI_INT,*m_adm.getCommunicator());
  hsize_t siz   = (hsize_t)std::accumulate(lens2,lens2+m_adm.getNoProcs(),0);
  hsize_t start = (hsize_t)std::accumulate(lens2,lens2+m_adm.getProcId(),0);
#else
  hsize_t siz   = (hsize_t)len;
  hsize_t start = 0;
#endif
  hid_t space = H5Screate_simple(1,&siz,nullptr);
  hid_t set = H5Dcreate2(group,name.c_str(),
                         type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if (len > 0) {
    hid_t file_space = H5Dget_space(set);
    siz = len;
    hsize_t stride = 1;
    H5Sselect_hyperslab(file_space,H5S_SELECT_SET,&start,&stride,&siz,nullptr);
    hid_t mem_space = H5Screate_simple(1,&siz,nullptr);
    H5Dwrite(set,type,mem_space,file_space,H5P_DEFAULT,data);
    H5Sclose(mem_space);
    H5Sclose(file_space);
  }
  H5Dclose(set);
  H5Sclose(space);
#else
  std::cout << "HDF5Restart: compiled without HDF5 support, no data written" << std::endl;
#endif
}


bool HDF5Restart::writeData(const TimeStep& tp,
                            const SerializeData& data)
{
  int level = tp.step / m_stride;

#ifdef HAS_HDF5
  int flag = H5F_ACC_RDWR;
  struct stat buffer;
  if (stat(m_hdf5_name.c_str(),&buffer) != 0)
    flag = H5F_ACC_TRUNC;

  if (!this->openFile(flag))
    return false;

  std::stringstream str;
  str << '/' << level;
  if (!checkGroupExistence(m_file,str.str().c_str()))
    H5Gclose(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));

  int pid = 0;
#ifdef HAVE_MPI
  pid = m_adm.getProcId();
  int ptot = m_adm.getNoProcs();
#else
  int ptot = 1;
#endif
  for (int p = 0; p < ptot; p++)
    for (const std::pair<std::string,std::string>& it : data) {
      std::stringstream str;
      str << level << '/' << p;
      hid_t group;
      if (checkGroupExistence(m_file,str.str().c_str()))
        group = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
      else
        group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
      if (!H5Lexists(group, it.first.c_str(), 0)) {
        int len = pid == p ? it.second.size() : 0;
        writeArray(group, it.first, len, it.second.data(), H5T_NATIVE_CHAR);
      }
      H5Gclose(group);
    }
  closeFile();
#else
  std::cout << "HDF5Restart: compiled without HDF5 support, no data written" << std::endl;
#endif

  return true;
}


//! \brief Convenience type for restart IO.
typedef std::pair<HDF5Restart*,HDF5Restart::SerializeData*> read_restart_ctx;


#ifdef HAS_HDF5
herr_t HDF5Restart::read_restart_data(hid_t group_id, const char* member_name, void* data)
{
  read_restart_ctx* ctx = static_cast<read_restart_ctx*>(data);

  char* c;
  int len;
  ctx->first->readArray(group_id,member_name,len,c);
  ctx->second->insert(std::make_pair(std::string(member_name),std::string(c,len)));
  return 0;
}
#endif


int HDF5Restart::readData(SerializeData& data, int level)
{
#ifdef HAS_HDF5
  if (!openFile(H5F_ACC_RDONLY))
    return -1;

  if (level == -1) {
    while (true) {
      std::stringstream str;
      str << '/' << level+1;
      if (!checkGroupExistence(m_file, str.str().c_str()))
        break;
      ++level;
    }
  }

  std::stringstream str;
  str << '/' << level << '/';
#ifdef HAVE_MPI
  str << m_adm.getProcId();
#else
  str << 0;
#endif
  int idx = 0;
  read_restart_ctx ctx(this,&data);
  int it = H5Giterate(m_file, str.str().c_str(), &idx, read_restart_data, &ctx);
  return it < 0 ? it : level;
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
  return -1;
#endif
}
