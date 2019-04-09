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
                        bool append)
  : DataWriter(name,adm,".hdf5"), HDF5Base(name+".hdf5", adm)
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


void HDF5Writer::openFile(int level)
{
  if (m_file != -1)
    return;

#ifdef HAS_HDF5
  if (m_flag == H5F_ACC_RDWR) {
    struct stat buffer;
    if (stat(m_name.c_str(),&buffer) != 0)
      m_flag = H5F_ACC_TRUNC;
  }

  if (m_flag == H5F_ACC_RDWR) {
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

  if (!HDF5Base::openFile(m_flag))
    return;

  std::stringstream str;
  str << '/' << level;
  if (!checkGroupExistence(m_file,str.str().c_str()))
    H5Gclose(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
#endif
}


void HDF5Writer::closeFile(int level)
{
#ifdef HAS_HDF5
  if (m_file) {
    H5Fflush(m_file,H5F_SCOPE_GLOBAL);
    H5Fclose(m_file);
  }
  m_flag = H5F_ACC_RDWR;
  m_file = -1;
#endif
}


void HDF5Writer::writeArray(hid_t group, const std::string& name, int patch,
                            int len, const void* data, hid_t type)
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
  hid_t space, set;
  hid_t group1 = -1;
  if (patch > -1) {
    if (checkGroupExistence(group, name.c_str()))
      group1 = H5Gopen2(group, name.c_str(),H5P_DEFAULT);
    else
      group1 = H5Gcreate2(group, name.c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    space = H5Screate_simple(1,&siz,nullptr);
    std::stringstream str;
    str << patch;
    set = H5Dcreate2(group1,str.str().c_str(),
                     type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  } else {
    space = H5Screate_simple(1,&siz,nullptr);
    set = H5Dcreate2(group,name.c_str(),
                     type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  }
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
  if (group1 != -1)
    H5Gclose(group1);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
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
    this->writeArray(level,entry.first,1,len,dvec->data(),H5T_NATIVE_DOUBLE);
  }
  else if (entry.second.field == DataExporter::INTVECTOR) {
    std::vector<int>* ivec = (std::vector<int>*)entry.second.data;
    int len = !redundant || rank == 0 ? ivec->size() : 0;
    this->writeArray(group,entry.first,1,len,ivec->data(),H5T_NATIVE_INT);
  }
  H5Gclose(group);
#endif
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
  if (!entry.second.enabled || !entry.second.data || entry.second.data2.empty())
    return;

#ifdef HAS_HDF5
  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const Vector* sol = static_cast<const Vector*>(entry.second.data2.front());
  const Vectors* proj = nullptr;
  if (entry.second.data2.size() > 1 && entry.second.data2[1])
    proj = static_cast<const Vectors*>(entry.second.data2[1]);
  const Matrix* eNorm = nullptr;
  if (entry.second.data2.size() > 2 && entry.second.data2[2])
    eNorm = static_cast<const Matrix*>(entry.second.data2[2]);
  const int results = abs(entry.second.results);

  if (prefix.empty() &&
      (level == 0 || geometryUpdated || (results & DataExporter::GRID)) &&
      entry.second.results > -1) {
    this->writeBasis(sim,sim->getName()+"-1",1,level);
    if (sim->mixedProblem())
      for (size_t b=2; b <= sim->getNoBasis(); ++b) {
        std::stringstream str;
        str << sim->getName() << "-" << b;
        this->writeBasis(sim,str.str(),b,level);
      }

    if (sim->fieldProjections() &&
        proj && !proj->empty() && !proj->front().empty()) {
      std::stringstream str;
      str << sim->getName() << "-proj";
      this->writeBasis(sim,str.str(),-1,level);
    }
  }

  NormBase* norm = sim->getNormIntegrand();
  const IntegrandBase* prob = sim->getProblem();
  bool usedescription = entry.second.results < 0;

  std::vector<hid_t> group(sim->getNoBasis(), -1);
  std::vector<hid_t> egroup;
  for (size_t b = 1; b <= sim->getNoBasis(); ++b) {
    std::stringstream str;
    str << level;
    str << '/' << sim->getName() << "-" << b;
    if (!checkGroupExistence(m_file,str.str().c_str())) {
      hid_t group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
      H5Gclose(group);
    }
    if (abs(results) & DataExporter::NORMS && norm) {
      std::stringstream str2;
      str2 << str.str() << "/knotspan";
      if (checkGroupExistence(m_file,str2.str().c_str()))
        egroup.push_back(H5Gopen2(m_file,str2.str().c_str(),H5P_DEFAULT));
      else
        egroup.push_back(H5Gcreate2(m_file,str2.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
    }
    if (abs(results) & (DataExporter::PRIMARY | DataExporter::SECONDARY) && !sol->empty()) {
      str << "/fields";
      if (checkGroupExistence(m_file,str.str().c_str()))
        group[b-1] = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
      else
        group[b-1] = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    }
  }

  if (sim->fieldProjections() &&
      proj && !proj->empty() && !proj->front().empty()) {
    std::stringstream str;
    str << level;
    str << '/' << sim->getName() << "-proj";
    str << "/fields";
    if (checkGroupExistence(m_file,str.str().c_str()))
      group.push_back(H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT));
    else
      group.push_back(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
  }


  size_t projOfs = 0;
  for (int i = 0; i < sim->getNoPatches(); ++i) {
    int loc = sim->getLocalPatchIndex(i+1);
    if (sim->getProcessAdm().dd.isPartitioned() && sim->getProcessAdm().getProcId() != 0)
      loc = 0;
    if (loc > 0 && (!(results & DataExporter::REDUNDANT) ||
                    sim->getGlobalProcessID() == 0)) // we own the patch
    {
      ASMbase* pch = sim->getPatch(loc);
      if (results & DataExporter::PRIMARY && !sol->empty()) {
        Vector psol;
        int ncmps = entry.second.ncmps;
        if (entry.second.results < 0) { // field assumed to be on basis 1 for now
          size_t ndof1 = sim->extractPatchSolution(*sol,psol,pch,ncmps,1);
          writeArray(group[0], entry.second.description, i+1,
                     ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
        } else {
          size_t ndof1 = sim->extractPatchSolution(*sol,psol,pch,ncmps);
          if (sim->mixedProblem())
          {
            size_t ofs = 0;
            for (size_t b=1; b <= sim->getNoBasis(); ++b) {
              ndof1 = pch->getNoNodes(b)*pch->getNoFields(b);
              writeArray(group[b-1],prefix+prob->getField1Name(10+b),i+1,ndof1,
                         psol.ptr()+ofs,H5T_NATIVE_DOUBLE);
              ofs += ndof1;
            }
          }
          else {
            writeArray(group[0], usedescription ? entry.second.description:
                                                prefix+prob->getField1Name(11), i+1,
                                         ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
          }
        }
      }

      if (results & DataExporter::SECONDARY && !sol->empty()) {
        Matrix field;
        SIM::SolutionMode mode = prob->getMode();
        const_cast<SIMbase*>(sim)->setMode(SIM::RECOVERY);
        sim->extractPatchSolution(Vectors(1,*sol), loc-1);
        sim->evalSecondarySolution(field,loc-1);
        const_cast<SIMbase*>(sim)->setMode(mode);

        for (size_t j = 0; j < field.rows(); j++)
          writeArray(group[0],prefix+prob->getField2Name(j),i+1,field.cols(),
                     field.getRow(j+1).ptr(),H5T_NATIVE_DOUBLE);
      }

      if (proj && !proj->empty()) {
        for (size_t p = 0; p < proj->size(); ++p) {
          if (proj->at(p).empty())
            continue;
          hid_t g = sim->fieldProjections() ? group.back() : group.front();
          Vector locvec;
          if (sim->fieldProjections()) {
            size_t ndof = sim->getPatch(loc)->getNoProjectionNodes() *
                          prob->getNoFields(2);
            locvec.resize(ndof);
            std::copy(proj->at(p).begin()+projOfs,
                      proj->at(p).begin()+projOfs+ndof, locvec.begin());
            if (p == proj->size()-1)
              projOfs += ndof;
          } else
            sim->extractPatchSolution(proj->at(p),locvec,pch,prob->getNoFields(2));

          Matrix field;
          field.resize(prob->getNoFields(2),locvec.size()/prob->getNoFields(2));
          field.fill(locvec.ptr());
          for (size_t j = 0; j < field.rows(); j++)
            writeArray(g,m_prefix[p]+" "+prob->getField2Name(j),i+1,field.cols(),
                       field.getRow(j+1).ptr(),H5T_NATIVE_DOUBLE);
        }
      }

      if (results & DataExporter::NORMS && eNorm) {
        Matrix patchEnorm;
        sim->extractPatchElmRes(*eNorm,patchEnorm,loc-1);
        for (size_t j = 1, l = 1; l < eNorm->rows(); j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              writeArray(egroup[0],
                         prefix+norm->getName(j,k,(j > 1 && !m_prefix.empty() ?
                                                m_prefix[j-2].c_str():0)),i+1,
                         patchEnorm.cols(),patchEnorm.getRow(l++).ptr(),
                         H5T_NATIVE_DOUBLE);
      }
      if (results & DataExporter::EIGENMODES) {
        const std::vector<Mode>* vec2 =
              static_cast<const std::vector<Mode>*>(entry.second.data2.front());
        const std::vector<Mode>& vec = static_cast<const std::vector<Mode>&>(*vec2);
        for (size_t k = 0; k < vec.size(); ++k) {
          Vector psol;
          size_t ndof1 = sim->extractPatchSolution(vec[k].eigVec,psol,pch);
          std::stringstream str;
          str << level;
          str << '/' << sim->getName() << "-1/Eigenmode";
          hid_t group2;
          if (checkGroupExistence(m_file,str.str().c_str()))
            group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
          else
            group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

          std::stringstream str4;
          str4 << k+1;
          writeArray(group2, str4.str(), i+1,
                     ndof1, psol.ptr(), H5T_NATIVE_DOUBLE);
          bool isFreq = sim->opt.eig==3 || sim->opt.eig==4 || sim->opt.eig==6;
          if (isFreq)
            writeArray(group2, str4.str()+"/Frequency", -1, 1, &vec[k].eigVal, H5T_NATIVE_DOUBLE);
          else
            writeArray(group2, str4.str()+"/Value", -1, 1, &vec[k].eigVal, H5T_NATIVE_DOUBLE);
          H5Gclose(group2);
        }
      }
    }
    else // must write empty dummy records for the other patches
    {
      double dummy=0.0;
      if (results & DataExporter::PRIMARY) {
        if (entry.second.results < 0) {
          writeArray(group[0], entry.second.description,i+1,
                     0, &dummy, H5T_NATIVE_DOUBLE);
        }
        else if (sim->mixedProblem())
        {
          for (size_t b=1; b <= sim->getNoBasis(); ++b)
            writeArray(group[b-1],prefix+prob->getField1Name(10+b),i+1,
                       0,&dummy,H5T_NATIVE_DOUBLE);
        }
        else
          writeArray(group[0], usedescription ? entry.second.description:
                                              prefix+prob->getField1Name(11),i+1,
                                              0, &dummy, H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::SECONDARY) {
        for (size_t j = 0; j < prob->getNoFields(2); j++)
          writeArray(group[0],prefix+prob->getField2Name(j),i+1,
                     0,&dummy,H5T_NATIVE_DOUBLE);

      }

      if (proj && !proj->empty()) {
        for (size_t p = 0; p < proj->size(); ++p) {
          for (size_t j = 0; j < prob->getNoFields(2); j++)
            writeArray(group[0],m_prefix[p]+" "+prob->getField2Name(j),i+1,
                       0,&dummy,H5T_NATIVE_DOUBLE);
        }
      }

      if (results & DataExporter::NORMS && eNorm)
        for (size_t j = 1; j <= norm->getNoFields(0); j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              writeArray(egroup[0],
                         prefix+norm->getName(j,k,(j > 1 && !m_prefix.empty() ?
                                                          m_prefix[j-2].c_str():0)),i+1,
                         0,&dummy,H5T_NATIVE_DOUBLE);
    }
  }
  delete norm;
  for (auto& it : group)
    if (it != -1) H5Gclose(it);
  for (auto& it : egroup)
    H5Gclose(it);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


void HDF5Writer::writeKnotspan (int level, const DataEntry& entry,
                                const std::string& prefix)
{
  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const Vector* sol = static_cast<const Vector*>(entry.second.data2.front());
  if (!sim || !sol) return;

  Matrix infield(1,sol->size());
  infield.fillRow(1,sol->ptr());

#ifdef HAS_HDF5
  std::stringstream str;
  str << level << '/' << sim->getName() << "-" << 1;
  if (!checkGroupExistence(m_file,str.str().c_str())) {
    hid_t group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    H5Gclose(group);
  }

  str << "/knotspan";
  hid_t group2;
  if (checkGroupExistence(m_file,str.str().c_str()))
    group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  else
    group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
  for (int i = 0; i < sim->getNoPatches(); ++i) {
    int loc = sim->getLocalPatchIndex(i+1);
    if (loc > 0 && (sim->getProcessAdm().isParallel() ||
                    sim->getGlobalProcessID() == 0)) // we own the patch
    {
      Matrix patchEnorm;
      sim->extractPatchElmRes(infield,patchEnorm,loc-1);
      writeArray(group2,prefix+entry.second.description,i+1,patchEnorm.cols(),
                 patchEnorm.getRow(1).ptr(),H5T_NATIVE_DOUBLE);
    }
    else { // must write empty dummy records for the other patches
      double dummy=0.0;
      writeArray(group2,prefix+entry.second.description,i+1,0,&dummy,H5T_NATIVE_DOUBLE);
    }

  }
  H5Gclose(group2);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


void HDF5Writer::writeBasis (const SIMbase* sim, const std::string& name,
                             int basis, int level, bool redundant)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << "/" << level << '/' << name;
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

  for (int i = 1; i <= sim->getNoPatches(); i++) {
    std::stringstream str;
    int loc = sim->getLocalPatchIndex(i);
    if (loc > 0)
      sim->getPatch(loc)->write(str,basis);
    if (!redundant || rank == 0)
      writeArray(group, "basis", i, str.str().size(), str.str().c_str(),
                 H5T_NATIVE_CHAR);
    if (redundant && rank != 0) {
      char dummy=0;
      writeArray(group, "basis", i, 0, &dummy, H5T_NATIVE_CHAR);
    }
  }
  H5Gclose(group);
#endif
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
  writeArray(group,"level",-1,toWrite,&tp.time.t,H5T_NATIVE_DOUBLE);
  H5Gclose(group);
#endif
  return true;
}


void HDF5Writer::writeNodalForces(int level, const DataEntry& entry)
{
#ifdef HAS_HDF5
  std::stringstream str;
  str << level << "/nodal";

  hid_t group2;
  if (!checkGroupExistence(m_file,str.str().c_str())) {
    hid_t group = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);
    H5Gclose(group);
  }

  str << "/" << entry.first;
  if (checkGroupExistence(m_file,str.str().c_str()))
    group2 = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
  else
    group2 = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

  if (m_rank == 0) {
    SIMbase* sim = static_cast<SIMbase*>(const_cast<void*>(entry.second.data));
    const GlbForceVec* forces = static_cast<const GlbForceVec*>(entry.second.data2.front());
    Vector values(forces->size()*3);
    Vector coords(forces->size()*3);
    Vec3 coord, val;
    for (size_t i = 0; i < forces->size(); i++)
    {
      int inod = forces->getForce(i,val);
      coord = sim->getNodeCoord(inod);
      for (int j=0;j<3;++j) {
        coords[i*3+j] = coord[j];
        values[i*3+j] = val[j];
      }
    }
    writeArray(group2,"values",-1,values.size(),values.data(),H5T_NATIVE_DOUBLE);
    writeArray(group2,"coords",-1,coords.size(),coords.data(),H5T_NATIVE_DOUBLE);
  } else {
    double dummy=0.0;
    writeArray(group2,"values",-1,0,&dummy,H5T_NATIVE_DOUBLE);
    writeArray(group2,"coords",-1,0,&dummy,H5T_NATIVE_DOUBLE);
  }
  H5Gclose(group2);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}
