// $Id$
//==============================================================================
//!
//! \file HDF5Writer.C
//!
//! \date Jul 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file.
//!
//==============================================================================

#include "HDF5Writer.h"
#include "GlbForceVec.h"
#include "SIMbase.h"
#include "ASMbase.h"
#include "IntegrandBase.h"
#include "TimeStep.h"
#include "Vec3.h"
#include <sstream>

#ifdef HAS_HDF5
#include <numeric>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <unistd.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#endif

//! \brief If file system has less than this amount free,
//! we bail to avoid corrupting file when a new write is initiated.
#define HDF5_SANITY_LIMIT 10*1024*1024LL // 10MB


HDF5Writer::HDF5Writer (const std::string& name, const ProcessAdm& adm,
                        bool append, bool keepOpen)
  : DataWriter(name,adm,".hdf5"), m_file(-1), m_restart_file(-1),
    m_keepOpen(keepOpen)
#ifdef HAVE_MPI
  , m_adm(adm)
#endif
{
#ifdef HAS_HDF5
  struct stat temp;
  // file already exists - open and find the next group
  if (append && keepOpen)
    m_flag = H5F_ACC_RDONLY;
  else if (append && stat(m_name.c_str(),&temp) == 0)
    m_flag = H5F_ACC_RDWR;
  else
    m_flag = H5F_ACC_TRUNC;
  m_restart_name = name + "_restart.hdf5";
  if (stat(m_restart_name.c_str(),&temp) == 0)
    m_restart_flag = H5F_ACC_RDWR;
  else
    m_restart_flag = H5F_ACC_TRUNC;
#endif
}


int HDF5Writer::getLastTimeLevel ()
{
  int result = 0;
#ifdef HAS_HDF5
  if (m_flag == H5F_ACC_TRUNC)
    return -1;

#ifdef HAVE_MPI
  MPI_Info info = MPI_INFO_NULL;
  hid_t acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl, MPI_COMM_SELF, info);
#else
  hid_t acc_tpl = H5P_DEFAULT;
#endif

  m_file = H5Fopen(m_name.c_str(),m_flag,acc_tpl);
  if (m_file < 0)
  {
    std::cerr <<" *** HDF5Writer: Failed to open "<< m_name << std::endl;
    return -2;
  }

  for (bool ok = true; ok; result++) {
    std::stringstream str;
    str << '/' << result;
    ok = checkGroupExistence(m_file,str.str().c_str());
  }
  H5Fclose(m_file);
  m_file = -1;
#endif
  return result-1;
}


void HDF5Writer::openFile(int level, bool restart)
{
  if ((m_file != -1 && !restart) || (m_restart_file != -1 && restart))
    return;
#ifdef HAS_HDF5
  hid_t file;
#ifdef HAVE_MPI
  MPI_Info info = MPI_INFO_NULL;
  hid_t acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl, *m_adm.getCommunicator(), info);
#else
  hid_t acc_tpl = H5P_DEFAULT;
#endif
  unsigned int flag = restart ? m_restart_flag : m_flag;
  std::string fname = restart ? m_restart_name : m_name;
  if (flag == H5F_ACC_RDWR) {
    struct stat buffer;
    if (stat(fname.c_str(),&buffer) != 0)
      flag = H5F_ACC_TRUNC;
  }

  if (flag == H5F_ACC_TRUNC)
    file = H5Fcreate(fname.c_str(),flag,H5P_DEFAULT,acc_tpl);
  else {
    // check free disk space - to protect against corrupting files
    // due to out of space condition
    if (flag == H5F_ACC_RDWR) {
#ifdef HAVE_GET_CURRENT_DIR_NAME
      char* cwd = get_current_dir_name();
      struct statvfs vfs;
      statvfs(cwd,&vfs);
      if (((int64_t)vfs.f_bavail)*vfs.f_bsize < HDF5_SANITY_LIMIT) {
        std::cerr << "HDF5Writer: Low disk space detected, bailing" << std::endl;
        exit(1);
      }
      free(cwd);
#endif
    }
    file = H5Fopen(fname.c_str(),flag,acc_tpl);
  }
  if (file < 0)
  {
    std::cerr <<" *** HDF5Writer: Failed to open "<< m_name << std::endl;
    return;
  }

  std::stringstream str;
  str << '/' << level;
  if (!checkGroupExistence(file,str.str().c_str()))
    H5Gclose(H5Gcreate2(file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
#ifdef HAVE_MPI
  H5Pclose(acc_tpl);
#endif
#endif
  if (restart)
    m_restart_file = file;
  else
   m_file = file;
}


void HDF5Writer::closeFile(int level, bool force)
{
  if (m_keepOpen && !force)
    return;
#ifdef HAS_HDF5
  if (m_file) {
    H5Fflush(m_file,H5F_SCOPE_GLOBAL);
    H5Fclose(m_file);
  }
  m_flag = H5F_ACC_RDWR;
  m_file = -1;
#endif
}


void HDF5Writer::readArray(hid_t group, const std::string& name,
                           int& len, double*& data)
{
#ifdef HAS_HDF5
  if (!checkGroupExistence(group, name.c_str())) {
    len = 0;
    return;
  }
  hid_t set = H5Dopen2(group,name.c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set) / 8;
  len = siz;
  data = new double[siz];
  H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(set);
#else
  len = 0;
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
#endif
}


void HDF5Writer::readArray(hid_t group, const std::string& name,
                           int& len,  int*& data)
{
#ifdef HAS_HDF5
  hid_t set = H5Dopen2(group,name.c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set) / sizeof(int);
  len = siz;
  data = new int[siz];
  H5Dread(set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(set);
#else
  len = 0;
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
#endif
}


void HDF5Writer::readArray(hid_t group, const std::string& name,
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
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
#endif
}


void HDF5Writer::readString(const std::string& name, std::string& out, bool close)
{
#ifdef HAS_HDF5
  openFile(0);
  hid_t set = H5Dopen2(m_file,name.c_str(),H5P_DEFAULT);
  hid_t space = H5Dget_space(set);
  hsize_t siz = H5Sget_simple_extent_npoints(space);
  char* temp = new char[siz];
  out.resize(siz);
  H5Dread(set,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
  out.assign(temp, siz);
  delete[] temp;
  H5Dclose(set);
  if (close)
    closeFile(0);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
#endif
}


void HDF5Writer::writeArray(hid_t group, const std::string& name,
                            int len, const void* data, int type)
{
#ifdef HAS_HDF5
#ifdef HAVE_MPI
  int lens[m_size], lens2[m_size];
  std::fill(lens,lens+m_size,len);
  MPI_Alltoall(lens,1,MPI_INT,lens2,1,MPI_INT,*m_adm.getCommunicator());
  hsize_t siz   = (hsize_t)std::accumulate(lens2,lens2+m_size,0);
  hsize_t start = (hsize_t)std::accumulate(lens2,lens2+m_rank,0);
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
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


bool HDF5Writer::readVector(int level, const DataEntry& entry)
{
//  readArray(level,entry.first,entry.second.size,entry.second.data);
  return true;
}


void HDF5Writer::writeVector(int level, const DataEntry& entry)
{
  if (!entry.second.enabled)
    return;
#ifdef HAS_HDF5
  int redundant = entry.second.results & DataExporter::REDUNDANT;
  int rank = 0;
#ifdef HAVE_MPI
  if (redundant)
    MPI_Comm_rank(*m_adm.getCommunicator(), &rank);
#endif
  std::stringstream str;
  str << level;
  hid_t group = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  if (entry.second.field == DataExporter::VECTOR) {
    Vector* dvec = (Vector*)entry.second.data;
    int len = !redundant || rank == 0 ? dvec->size() : 0;
    this->writeArray(level,entry.first,len,dvec->data(),H5T_NATIVE_DOUBLE);
  }
  else if (entry.second.field == DataExporter::INTVECTOR) {
    std::vector<int>* ivec = (std::vector<int>*)entry.second.data;
    int len = !redundant || rank == 0 ? ivec->size() : 0;
    this->writeArray(group,entry.first,len,ivec->data(),H5T_NATIVE_INT);
  }
  H5Gclose(group);
#endif
}


bool HDF5Writer::readVector(int level, const std::string& name,
                            int patch, std::vector<double>& vec)
{
  bool ok=true;
  openFile(level);
  vec.clear();
#ifdef HAS_HDF5
  std::stringstream str;
  str << level;
  if (patch > -1) {
    str << '/';
    str << patch;
  }
  hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  double* tmp = nullptr; int siz = 0;
  readArray(group2,name,siz,tmp);
  vec.insert(vec.begin(),tmp,tmp+siz);
  delete[] tmp;
  H5Gclose(group2);
#endif
  closeFile(level);
  return ok;
}


void HDF5Writer::writeBasis (int level, const DataEntry& entry,
                             const std::string& prefix)
{
  if (!entry.second.enabled || !entry.second.data)
    return;

  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);

  this->writeBasis(sim, prefix+sim->getName()+"-1", 1, level,
                   entry.second.results & DataExporter::REDUNDANT);
}


void HDF5Writer::writeSIM (int level, const DataEntry& entry,
                           bool geometryUpdated, const std::string& prefix)
{
  if (!entry.second.enabled || !entry.second.data || !entry.second.data2)
    return;

  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const Vector* sol = static_cast<const Vector*>(entry.second.data2);
  const int results = abs(entry.second.results);

  if (level == 0 || geometryUpdated || (results & DataExporter::GRID)) {
    this->writeBasis(sim,prefix+sim->getName()+"-1",1,level);
    if (sim->mixedProblem())
      for (size_t b=2; b <= sim->getNoBasis(); ++b) {
        std::stringstream str;
        str << sim->getName() << "-" << b;
        this->writeBasis(sim,str.str(),b,level);
      }
  }

  Matrix eNorm;
  Vectors gNorm;
  if (results & DataExporter::NORMS)
    const_cast<SIMbase*>(sim)->solutionNorms(*sol,eNorm,gNorm);

#ifdef HAS_HDF5
  NormBase* norm = sim->getNormIntegrand();
  const IntegrandBase* prob = sim->getProblem();
  bool usedescription = entry.second.results < 0;

  for (int i = 0; i < sim->getNoPatches(); ++i) {
    std::stringstream str;
    str << level;
    str << '/';
    str << i+1;
    hid_t group2;
    if (checkGroupExistence(m_file,str.str().c_str()))
      group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
    else
      group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0 && (!(results & DataExporter::REDUNDANT) ||
                    sim->getGlobalProcessID() == 0)) // we own the patch
    {
      ASMbase* pch = sim->getPatch(loc);
      if (results & DataExporter::PRIMARY) {
        Vector psol;
        int ncmps = entry.second.ncmps;
        if (entry.second.results < 0) { // field assumed to be on basis 1 for now
          size_t ndof1 = sim->extractPatchSolution(*sol,psol,pch,ncmps,1);
          writeArray(group2, entry.second.description,
                     ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
        } else {
          size_t ndof1 = sim->extractPatchSolution(*sol,psol,pch,ncmps);
          if (sim->mixedProblem())
          {
            size_t ofs = 0;
            for (size_t b=1; b <= sim->getNoBasis(); ++b) {
              ndof1 = pch->getNoNodes(b)*pch->getNoFields(b);
              writeArray(group2,prefix+prob->getField1Name(10+b),ndof1,
                         psol.ptr()+ofs,H5T_NATIVE_DOUBLE);
              ofs += ndof1;
            }
          }
          else {
            writeArray(group2, usedescription ? entry.second.description:
                                                prefix+prob->getField1Name(11),
                                         ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
          }
        }
      }

      if (results & DataExporter::SECONDARY) {
        Matrix field;
        if (prefix.empty()) {
          const_cast<SIMbase*>(sim)->setMode(SIM::RECOVERY);
          sim->extractPatchSolution(Vectors(1,*sol), loc-1);
          sim->evalSecondarySolution(field,loc-1);
        }
        else {
          Vector locvec;
          sim->extractPatchSolution(*sol,locvec,pch,prob->getNoFields(2));
          field.resize(prob->getNoFields(2),locvec.size()/prob->getNoFields(2));
          field.fill(locvec.ptr());
        }
        for (size_t j = 0; j < field.rows(); j++)
          writeArray(group2,prefix+prob->getField2Name(j),field.cols(),
                     field.getRow(j+1).ptr(),H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::NORMS && norm) {
        Matrix patchEnorm;
        sim->extractPatchElmRes(eNorm,patchEnorm,loc-1);
        for (size_t j = 1, l = 1; j <= norm->getNoFields(0); j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              writeArray(group2,
                         prefix+norm->getName(j,k,(j>1&&m_prefix?m_prefix[j-2]:0)),
                         patchEnorm.cols(),patchEnorm.getRow(l++).ptr(),
                         H5T_NATIVE_DOUBLE);
      }
      if (results & DataExporter::EIGENMODES) {
        const std::vector<Mode>* vec2 = static_cast<const std::vector<Mode>*>(entry.second.data2);
        const std::vector<Mode>& vec = static_cast<const std::vector<Mode>&>(*vec2);
        H5Gclose(group2);
        for (size_t k = 0; k < vec.size(); ++k) {
          Vector psol;
          size_t ndof1 = sim->extractPatchSolution(vec[k].eigVec,psol,pch);
          std::stringstream name;
          name << entry.second.description << "-" << k+1;
          std::stringstream str;
          str << k;
          if (!checkGroupExistence(m_file,str.str().c_str())) {
            group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
            H5Gclose(group2);
          }
          str << '/';
          str << i+1;
          if (checkGroupExistence(m_file,str.str().c_str()))
            group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
          else
            group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
          writeArray(group2, "eigenmode",
                     ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
          bool isFreq = sim->opt.eig==3 || sim->opt.eig==4 || sim->opt.eig==6;
          if (isFreq)
            writeArray(group2, "eigenfrequency", 1, &vec[k].eigVal, H5T_NATIVE_DOUBLE);
          else
            writeArray(group2, "eigenval", 1, &vec[k].eigVal, H5T_NATIVE_DOUBLE);
          if (k != vec.size()-1)
            H5Gclose(group2);
        }
      }
    }
    else // must write empty dummy records for the other patches
    {
      double dummy=0.0;
      if (results & DataExporter::PRIMARY) {
        if (entry.second.results < 0) {
          writeArray(group2, entry.second.description,
                     0, &dummy, H5T_NATIVE_DOUBLE);
        }
        else if (sim->mixedProblem())
        {
          for (size_t b=1; b <= sim->getNoBasis(); ++b)
            writeArray(group2,prefix+prob->getField1Name(10+b),0,&dummy,H5T_NATIVE_DOUBLE);
        }
        else
          writeArray(group2, usedescription ? entry.second.description:
                                              prefix+prob->getField1Name(11),
                                              0, &dummy, H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::SECONDARY)
        for (size_t j = 0; j < prob->getNoFields(2); j++)
          writeArray(group2,prefix+prob->getField2Name(j),0,&dummy,H5T_NATIVE_DOUBLE);

      if (results & DataExporter::NORMS && norm)
        for (size_t j = 1; j <= norm->getNoFields(0); j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              writeArray(group2,
                         prefix+norm->getName(j,k,(j>1&&m_prefix?m_prefix[j-2]:0)),
                         0,&dummy,H5T_NATIVE_DOUBLE);
    }
    H5Gclose(group2);
  }
  delete norm;
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


void HDF5Writer::writeKnotspan (int level, const DataEntry& entry,
                                const std::string& prefix)
{
  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const Vector* sol = static_cast<const Vector*>(entry.second.data2);
  if (!sim || !sol) return;

  Matrix infield(1,sol->size());
  infield.fillRow(1,sol->ptr());

#ifdef HAS_HDF5
  for (int i = 0; i < sim->getNoPatches(); ++i) {
    std::stringstream str;
    str << level;
    str << '/';
    str << i+1;
    hid_t group2;
    if (checkGroupExistence(m_file,str.str().c_str()))
      group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
    else
      group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0 && (sim->getProcessAdm().isParallel() ||
                    sim->getGlobalProcessID() == 0)) // we own the patch
    {
      Matrix patchEnorm;
      sim->extractPatchElmRes(infield,patchEnorm,loc-1);
      writeArray(group2,prefix+entry.second.description,patchEnorm.cols(),
                 patchEnorm.getRow(1).ptr(),H5T_NATIVE_DOUBLE);
    }
    else { // must write empty dummy records for the other patches
      double dummy=0.0;
      writeArray(group2,prefix+entry.second.description,0,&dummy,H5T_NATIVE_DOUBLE);
    }

    H5Gclose(group2);
  }
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


void HDF5Writer::writeBasis (const SIMbase* sim, const std::string& name,
                             int basis, int level, bool redundant)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << "/" << level << "/basis";
  hid_t group;
  int rank = 0;
#ifdef HAVE_MPI
  if (redundant)
    MPI_Comm_rank(*m_adm.getCommunicator(), &rank);
#endif
  if (checkGroupExistence(m_file,str.str().c_str()))
    group = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  else
    group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
  str << "/" << name;
  if (checkGroupExistence(m_file,str.str().c_str()))
  {
    H5Gclose(group);
    return;
  }
  hid_t group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

  for (int i = 1; i <= sim->getNoPatches(); i++) {
    std::stringstream str, str2;
    int loc = sim->getLocalPatchIndex(i);
    if (loc > 0)
      sim->getPatch(loc)->write(str,basis);
    str2 << i;
    if (!redundant || rank == 0)
      writeArray(group2, str2.str(), str.str().size(), str.str().c_str(),
                 H5T_NATIVE_CHAR);
    if (redundant && rank != 0) {
      char dummy=0;
      writeArray(group2, str2.str(), 0, &dummy, H5T_NATIVE_CHAR);
    }
  }
  H5Gclose(group2);
  H5Gclose(group);
#endif
}


bool HDF5Writer::checkGroupExistence(hid_t parent, const char* path)
{
  bool result = false;
#ifdef HAS_HDF5
  // turn off errors to avoid cout spew
  H5E_BEGIN_TRY {
#if H5_VERS_MINOR > 8
  result = H5Lexists(parent, path, H5P_DEFAULT) == 1;
#else
    result = H5Gget_objinfo((hid_t)parent,path,0,nullptr) == 0;
#endif
  } H5E_END_TRY;
#endif
  return result;
}


// TODO: implement for variable time steps.
// (named time series to allow different timelevels for different fields)
bool HDF5Writer::writeTimeInfo (int level, int interval, const TimeStep& tp)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << "/" << level << "/timeinfo";
  hid_t group;
  if (checkGroupExistence(m_file,str.str().c_str()))
    group = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  else
    group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

  // parallel nodes != 0 write dummy entries
  int toWrite=(m_rank == 0);

  // !TODO: different names
  writeArray(group,"SIMbase-1",toWrite,&tp.time.t,H5T_NATIVE_DOUBLE);
  H5Gclose(group);
#endif
  return true;
}


void HDF5Writer::writeNodalForces(int level, const DataEntry& entry)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << level;
  hid_t group2;
  if (checkGroupExistence(m_file,str.str().c_str()))
    group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  else
    group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

  if (m_rank == 0) {
    SIMbase* sim = static_cast<SIMbase*>(const_cast<void*>(entry.second.data));
    const GlbForceVec* forces = static_cast<const GlbForceVec*>(entry.second.data2);
    Vector results(forces->size()*6);
    Vec3 coord, val;
    for (size_t i = 0; i < forces->size(); i++)
    {
      int inod = forces->getForce(i,val);
      coord = sim->getNodeCoord(inod);
      for (int j=0;j<3;++j) {
        results[i*6+j] = coord[j];
        results[i*6+j+3] = val[j];
      }
    }
    writeArray(group2,entry.first,results.size(),results.data(),H5T_NATIVE_DOUBLE);
  } else {
    double dummy=0.0;
    writeArray(group2,entry.first,0,&dummy,H5T_NATIVE_DOUBLE);
  }
  H5Gclose(group2);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


bool HDF5Writer::hasGeometries(int level, const std::string& basisName)
{
  std::stringstream str;
  str << '/' << level << "/basis";
  if (!basisName.empty())
    str << '/' << basisName;

  bool notOpen = m_file == -1;
  if (notOpen) this->openFile(level);
  bool result = this->checkGroupExistence(m_file,str.str().c_str());
  if (notOpen) this->closeFile(level);

  return result;
}


bool HDF5Writer::readDouble(int level, const std::string& group,
                            const std::string& name, double& data)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << "/" << level << "/" << group;
  if (!checkGroupExistence(m_file,str.str().c_str()))
    return false;

  hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  double* data2; int len=1;
  readArray(group2,name,len,data2);
  H5Gclose(group2);
  if (len > 0) {
    data = data2[0];
    delete[] data2;
    return true;
  }
#endif
  return false;
}


bool HDF5Writer::readVector(int level, const std::string& name,
                            int patch, std::vector<int>& vec)
{
  bool ok=true;
  openFile(level);
#ifdef HAS_HDF5
  std::stringstream str;
  str << level;
  if (patch > -1) {
    str << '/';
    str << patch;
  }
  hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  int* tmp = nullptr; int siz = 0;
  readArray(group2,name,siz,tmp);
  vec.assign(tmp, tmp + siz);
  delete[] tmp;
  H5Gclose(group2);
#endif
  closeFile(level);
  return ok;
}


bool HDF5Writer::writeRestartData(int level, const DataExporter::SerializeData& data)
{
#ifdef HAS_HDF5
  if (level > 0)
    m_restart_flag = H5F_ACC_RDWR;
  openFile(level, true);
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
      if (checkGroupExistence(m_restart_file,str.str().c_str()))
        group = H5Gopen2(m_restart_file,str.str().c_str(),H5P_DEFAULT);
      else
        group = H5Gcreate2(m_restart_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
      if (!H5Lexists(group, it.first.c_str(), 0)) {
        int len = pid == p ? it.second.size() : 0;
        writeArray(group, it.first, len, it.second.data(), H5T_NATIVE_CHAR);
      }
      H5Gclose(group);
    }

  H5Fclose(m_restart_file);
  m_restart_file = -1;
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif

  return true;
}


//! \brief Convenience type for restart IO.
typedef std::pair<HDF5Writer*,DataExporter::SerializeData*> read_restart_ctx;


#ifdef HAS_HDF5
static herr_t read_restart_data(hid_t group_id, const char* member_name, void* data)
{
  read_restart_ctx* ctx = static_cast<read_restart_ctx*>(data);

  char* c;
  int len;
  ctx->first->readArray(group_id,member_name,len,c);
  ctx->second->insert(std::make_pair(std::string(member_name),std::string(c,len)));
  return 0;
}
#endif


int HDF5Writer::readRestartData(DataExporter::SerializeData& data, int level)
{
#ifdef HAS_HDF5
  openFile(0);

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
  closeFile(0);
  return it < 0 ? it : level;
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
  return -1;
#endif
}
