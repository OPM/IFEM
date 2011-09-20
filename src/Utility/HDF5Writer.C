// $Id$
//==============================================================================
//!
//! \file HDF5Writer.C
//!
//! \date Jul 7 2011
//!
//! \author Arne MOrten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file.
//!
//==============================================================================

#include "HDF5Writer.h"
#include "SIMbase.h"
#include "SIMparameters.h"
#include "IntegrandBase.h"
#include <sstream>
#include <sys/stat.h>
#include <numeric>

#ifdef HAS_HDF5
#include <hdf5.h>
#endif

#ifdef PARALLEL_PETSC
#include <mpi.h>
#endif


HDF5Writer::HDF5Writer (const std::string& name, bool append, bool keepOpen)
  : DataWriter(name+".hdf5"), m_file(0), m_keepOpen(keepOpen)
{
#ifdef HAS_HDF5
  struct stat temp;
  // file already exists - open and find the next group
  if (append && stat(m_name.c_str(),&temp) == 0)
    m_flag = H5F_ACC_RDWR;
  else
    m_flag = H5F_ACC_TRUNC;
#endif
}


int HDF5Writer::getLastTimeLevel ()
{
  int result = 0;
#ifdef HAS_HDF5
  if (m_flag == H5F_ACC_TRUNC)
    return -1;

  hid_t acc_tpl = H5P_DEFAULT;
#ifdef PARALLEL_PETSC
  MPI_Info info = MPI_INFO_NULL;
  acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl, MPI_COMM_WORLD, info);
#endif

  m_file = H5Fopen(m_name.c_str(),m_flag,acc_tpl);
  while (1) {
    std::stringstream str;
    str << '/' << result;
    if (!checkGroupExistence(m_file,str.str().c_str()))
      break;
    result++;
  }
  H5Fclose(m_file);
  m_file = 0;
#endif
  return result-1;
}


void HDF5Writer::openFile(int level)
{
  if (m_file)
    return;
#ifdef HAS_HDF5
  hid_t acc_tpl = H5P_DEFAULT;
#ifdef PARALLEL_PETSC
  MPI_Info info = MPI_INFO_NULL;
  acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl, MPI_COMM_WORLD, info);
#endif

  if (m_flag == H5F_ACC_TRUNC)
    m_file = H5Fcreate(m_name.c_str(),m_flag,H5P_DEFAULT,acc_tpl);
  else
    m_file = H5Fopen(m_name.c_str(),m_flag,acc_tpl);

  std::stringstream str;
  str << '/' << level;
  if (!checkGroupExistence(m_file,str.str().c_str()))
    H5Gclose(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
#ifdef PARALLEL_PETSC
  H5Pclose(acc_tpl);
#endif
#endif
}

void HDF5Writer::closeFile(int level, bool force)
{
  if (m_keepOpen && !force)
    return;
#ifdef HAS_HDF5
  H5Fflush(m_file,H5F_SCOPE_GLOBAL);
  H5Fclose(m_file);
  m_flag = H5F_ACC_RDWR;
  m_file = 0;
#endif
}

void HDF5Writer::readArray(int group, const std::string& name,
                           int& len, double*& data)
{
#ifdef HAS_HDF5
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

void HDF5Writer::readString(const std::string& name, std::string& out)
{
#ifdef HAS_HDF5
  openFile(0);
  hid_t set = H5Dopen2(m_file,name.c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set);
  char* temp = new char[siz+1];
  out.resize(siz);
  H5Dread(set,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
  temp[siz] = '\0';
  out = temp;
  delete[] temp;
  H5Dclose(set);
  closeFile(0);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data read" << std::endl;
#endif
}

void HDF5Writer::writeArray(int group, const std::string& name,
                            int len, const void* data, int type)
{
#ifdef HAS_HDF5
#ifdef PARALLEL_PETSC
  int lens[m_size], lens2[m_size];
  std::fill(lens,lens+m_size,len);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(lens,1,MPI_INT,lens2,1,MPI_INT,MPI_COMM_WORLD);
  hsize_t siz   = (hsize_t)std::accumulate(lens2,lens2+m_size,0);
  hsize_t start = (hsize_t)std::accumulate(lens2,lens2+m_rank,0);
#else
  hsize_t siz   = (hsize_t)len;
  hsize_t start = 0;
#endif
  hid_t space = H5Screate_simple(1,&siz,NULL);
  hid_t set = H5Dcreate2(group,name.c_str(),
                         type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if (len > 0) {
    hid_t file_space = H5Dget_space(set);
    siz = len;
    hsize_t stride = 1;
    H5Sselect_hyperslab(file_space,H5S_SELECT_SET,&start,&stride,&siz,NULL);
    hid_t mem_space = H5Screate_simple(1,&siz,NULL);
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
#ifdef HAS_HDF5
  Vector* vector = (Vector*)entry.second.data;
  writeArray(level,entry.first,vector->size(),vector->data(),H5T_NATIVE_DOUBLE);
#endif
}

bool HDF5Writer::readSIM (int level, const DataEntry& entry)
{
  bool ok = true;
#ifdef HAS_HDF5
  SIMbase* sim = static_cast<SIMbase*>(entry.second.data);
  Vector* sol = static_cast<Vector*>(entry.second.data2);
  if (!sol) return false;
  const IntegrandBase* prob = sim->getProblem();

  for (int i = 0; i < sim->getNoPatches() && ok; ++i) {
    std::stringstream str;
    str << '/' << level << '/' << i+1;
    hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0) {
      double* tmp = NULL; int siz = 0;
      if (prob->mixedFormulation())
        readArray(group2,entry.first,siz,tmp);
      else
        readArray(group2,prob->getField1Name(11),siz,tmp);
      ok = sim->injectPatchSolution(*sol,Vector(tmp,siz),
                                    loc-1);
      delete[] tmp;
    }
    H5Gclose(group2);
  }
#endif
  return ok;
}

bool HDF5Writer::readVector(int level, const std::string& name,
                            int patch, Vector& vec)
{
  bool ok=true;
  openFile(level);
#ifdef HAS_HDF5
  std::stringstream str;
  str << level;
  str << '/';
  str << patch;
  hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  double* tmp = NULL; int siz = 0;
  readArray(group2,name,siz,tmp);
  vec.fill(tmp,siz);
  delete[] tmp;
  H5Gclose(group2);
#endif
  closeFile(level);
  return ok;
}

bool HDF5Writer::readField(int level, const std::string& name,
                           Vector& vec, SIMbase* sim, int components)
{
  bool ok = true;
  openFile(level);
#ifdef HAS_HDF5
  vec.resize(sim->getNoNodes()*components);
  for (int i = 0; i < sim->getNoPatches() && ok; ++i) {
    std::stringstream str;
    str << level;
    str << '/';
    str << i+1;
    hid_t group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0) {
      double* tmp = NULL; int siz = 0;
      readArray(group2,name,siz,tmp);
      ok = sim->injectPatchSolution(vec,Vector(tmp,siz),
                                    loc-1,components);
      delete[] tmp;
    }
    H5Gclose(group2);
  }
#endif
  closeFile(level);
  return ok;
}

void HDF5Writer::writeSIM (int level, const DataEntry& entry)
{
#ifdef HAS_HDF5
  SIMbase* sim = static_cast<SIMbase*>(entry.second.data);
  Vector* sol = static_cast<Vector*>(entry.second.data2);
  if (!sol) return;

  const IntegrandBase* prob = sim->getProblem();
  if (level == 0) { // TODO time dependent geometries
    writeBasis(sim,sim->getName()+"-1",1,level);
    if (prob->mixedFormulation())
      writeBasis(sim,sim->getName()+"-2",2,level);
  }

  Matrix eNorm;
  Vector gNorm;
  if (entry.second.results & DataExporter::NORMS)
    sim->solutionNorms(Vectors(1,*sol),eNorm,gNorm);

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
    if (loc > 0) // we own the patch
    {
      size_t ndof1 = sim->extractPatchSolution(*sol,loc-1);
      Vector& psol = const_cast<IntegrandBase*>(prob)->getSolution();
      if (prob->mixedFormulation())
      {
        // Mixed methods: The primary solution vector is referring to two bases
        size_t ndof2 = psol.size() > ndof1 ? psol.size() - ndof1 : 0;
        writeArray(group2,entry.first,psol.size(),psol.ptr(),H5T_NATIVE_DOUBLE);
        writeArray(group2,prob->getField1Name(11),ndof1,psol.ptr(),H5T_NATIVE_DOUBLE);
        writeArray(group2,prob->getField1Name(12),ndof2,psol.ptr()+ndof1,H5T_NATIVE_DOUBLE);
      }
      else
        writeArray(group2,prob->getField1Name(11),psol.size(),psol.ptr(),H5T_NATIVE_DOUBLE);

      if (entry.second.results & DataExporter::SECONDARY) {
        Matrix field;
        sim->evalSecondarySolution(field,loc-1);
        for (size_t j = 0; j < prob->getNoFields(2); ++j)
          writeArray(group2,prob->getField2Name(j),field.cols(),
                     field.getRow(j+1).ptr(),H5T_NATIVE_DOUBLE);
      }
      if (entry.second.results & DataExporter::NORMS) {
        Matrix patchEnorm;
        sim->extractPatchElmRes(eNorm,patchEnorm,loc-1);
        for (size_t j=0;j<eNorm.rows();++j) {
          if (NormBase::hasElementContributions(j)) {
            Vector k;
            writeArray(group2,NormBase::getName(j),patchEnorm.cols(),
                       patchEnorm.getRow(1+j).ptr(),H5T_NATIVE_DOUBLE);
          }
        }
      }
    }
    else // must write empty dummy records for the other patches
    {
      double dummy;
      writeArray(group2,prob->getField1Name(11),0,&dummy,H5T_NATIVE_DOUBLE);
      if (prob->mixedFormulation())
      {
        writeArray(group2,prob->getField1Name(11),0,&dummy,H5T_NATIVE_DOUBLE);
        writeArray(group2,prob->getField1Name(12),0,&dummy,H5T_NATIVE_DOUBLE);
      }
      if (entry.second.results & DataExporter::SECONDARY) {
        for (size_t j = 0; j < prob->getNoFields(2); ++j)
          writeArray(group2,prob->getField2Name(j),0,&dummy,H5T_NATIVE_DOUBLE);
      }
      if (entry.second.results & DataExporter::NORMS) {
        for (size_t j=0;j<eNorm.rows();++j) {
          if (NormBase::hasElementContributions(j))
            writeArray(group2,NormBase::getName(j),0,&dummy,H5T_NATIVE_DOUBLE);
        }
      }
    }
    H5Gclose(group2);
  }
#endif
}

void HDF5Writer::writeBasis(SIMbase* sim, const std::string& name,
                               int basis, int level)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << "/" << level << "/basis";
  int group;
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

  for (int i=0;i<sim->getNoPatches();++i) {
    std::stringstream str;
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0)
      sim->dumpBasis(str,basis,loc);
    std::stringstream str2;
    str2 << i+1;
    writeArray(group2, str2.str(), str.str().size(), str.str().c_str(),
               H5T_NATIVE_CHAR);
  }
  H5Gclose(group2);
  H5Gclose(group);
#endif
}

bool HDF5Writer::checkGroupExistence(int parent, const char* path)
{
  bool result = false;
#ifdef HAS_HDF5
  // turn off errors to avoid cout spew
  H5E_BEGIN_TRY {
    result = H5Gget_objinfo((hid_t)parent,path,0,NULL) == 0;
  } H5E_END_TRY;
#endif
  return result;
}

// TODO: implement for variable timesteps 
// (named time series to allow different timelevels for different fields)
bool HDF5Writer::writeTimeInfo(int level, int order, int interval,
                               SIMparameters& tp)
{
  return true;
}
