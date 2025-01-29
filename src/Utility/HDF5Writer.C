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
#include "ASMsupel.h"
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
#ifdef HAS_HDF5
  if (m_file != -1)
    return;

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


#ifdef HAS_HDF5
void HDF5Writer::writeArray(hid_t group, const std::string& name, int patch,
                            int len, const void* data, hid_t type)
{
#if SP_DEBUG > 2
  std::cout <<"HDF5Writer::writeArray: "<< name <<" for patch "<< patch
            <<" size="<< len << std::endl;
#endif
#ifdef HAVE_MPI
  hsize_t siz;
  hsize_t start;
  if (m_adm.dd.isPartitioned()) {
    siz = static_cast<hsize_t>(len);
    start = 0;
  } else {
    std::vector<int> lens(m_size), lens2(m_size);
    std::fill(lens.begin(), lens.end(), len);
    MPI_Alltoall(lens.data(),1,MPI_INT,lens2.data(),1,MPI_INT,*m_adm.getCommunicator());
    siz   = std::accumulate(lens2.begin(), lens2.end(), hsize_t(0));
    start = std::accumulate(lens2.begin(), lens2.begin() + m_rank, hsize_t(0));
  }
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
}
#endif


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

#ifdef HAS_HDF5
  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  this->writeBasis(sim, prefix+sim->getName()+"-1", 1, level,
                   entry.second.results & DataExporter::REDUNDANT,
                   entry.second.results & DataExporter::L2G_NODE);
#else
  std::cout <<"HDF5Writer: Compiled without HDF5 support, no data written."<< std::endl;
#endif
}


void HDF5Writer::writeSIM (int level, const DataEntry& entry,
                           bool geometryUpdated, const std::string& prefix)
{
  if (!entry.second.enabled || !entry.second.data || entry.second.data2.empty())
    return;

#ifdef HAVE_MPI
  if (m_adm.dd.isPartitioned() && m_rank != 0) // only rank 0 write
    return;
#endif

#ifdef HAS_HDF5
  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const Vector*  sol = static_cast<const Vector*>(entry.second.data2.front());
  const Vectors* proj = nullptr;
  const Matrix* eNorm = nullptr;
  if (entry.second.data2.size() > 1 && entry.second.data2[1])
    if ((proj = static_cast<const Vectors*>(entry.second.data2[1])))
      if (proj->empty() || proj->front().empty())
        proj = nullptr;
  if (entry.second.data2.size() > 2 && entry.second.data2[2])
    eNorm = static_cast<const Matrix*>(entry.second.data2[2]);

#if SP_DEBUG > 1
  std::cout <<"\n >>> Writing results to HDF5 for simulator "<< sim->getName()
            <<" at time level "<< level <<" <<<"<< std::endl;
  std::cout <<"\nProvided solution vector:"<< *sol;
  if (proj)
    for (const Vector& psol : *proj)
      std::cout <<"\nProjected solution vector:"<< psol;
  if (eNorm)
    std::cout <<"\nError norms:"<< *eNorm;
#endif

  bool usedescription = entry.second.results < 0;
  const int results = abs(entry.second.results);
  if (!(results & DataExporter::NORMS))
    eNorm = nullptr;

  if (prefix.empty() && !usedescription &&
      (level == 0 || geometryUpdated || (results & DataExporter::GRID)))
  {
    std::string bName = sim->getName() + "-1";
    for (size_t b = 1; b <= sim->getNoBasis(); b++, ++bName.back())
      this->writeBasis(sim,bName,b,level,results & DataExporter::REDUNDANT,results & DataExporter::L2G_NODE);
    if (proj && sim->fieldProjections())
      this->writeBasis(sim, sim->getName() + "-proj", -1, level);
  }

  const IntegrandBase* prob = sim->getProblem();
  NormBase* norm = eNorm ? sim->getNormIntegrand() : nullptr;
  size_t normGrp = norm ? norm->getNoFields(0) : 0;

  // Extract names for the projection methods used
  std::vector<const char*> projPfx;
  if (proj || norm)
  {
    projPfx.reserve(sim->opt.project.size());
    for (const auto& pfx : sim->opt.project)
      projPfx.push_back(pfx.second.c_str());
    if (normGrp > 1+projPfx.size())
      projPfx.resize(normGrp-1,nullptr);
    else
      normGrp = 1+projPfx.size();
  }

  std::vector<hid_t> egroup, group;
  egroup.reserve(sim->getNoBasis());
  group.reserve(sim->getNoBasis()+1);
  for (size_t b = 1; b <= sim->getNoBasis(); ++b) {
    std::stringstream str;
    str << level;
    str << '/' << sim->getName() << "-" << b;
    if (!checkGroupExistence(m_file,str.str().c_str()))
      H5Gclose(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
    if (norm) {
      std::stringstream str2;
      str2 << str.str() << "/knotspan";
      if (checkGroupExistence(m_file,str2.str().c_str()))
        egroup.push_back(H5Gopen2(m_file,str2.str().c_str(),H5P_DEFAULT));
      else
        egroup.push_back(H5Gcreate2(m_file,str2.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
    }
    if (results & (DataExporter::PRIMARY | DataExporter::SECONDARY) && !sol->empty()) {
      str << "/fields";
      if (checkGroupExistence(m_file,str.str().c_str()))
        group.push_back(H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT));
      else
        group.push_back(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
    }
  }

  if (proj && sim->fieldProjections()) {
    std::stringstream str;
    str << level;
    str << '/' << sim->getName() << "-proj";
    str << "/fields";
    if (checkGroupExistence(m_file,str.str().c_str()))
      group.push_back(H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT));
    else
      group.push_back(H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT));
  }

  // Lambda function for extracting projected solution for a patch.
  size_t projOfs = 0;
  auto&& extractProjection=[sim,&projOfs](ASMbase* pch, const Vector& glbVec,
                                          Matrix& field, bool isLastProj = false)
  {
    size_t ncmp = sim->getProblem()->getNoFields(2);
    if (ncmp == 0)
      return;

    Vector pchVec;
    if (sim->fieldProjections())
    {
      pchVec.resize(ncmp * pch->getNoProjectionNodes());
      size_t projEnd = projOfs + pchVec.size();
      std::copy(glbVec.begin()+projOfs, glbVec.begin()+projEnd, pchVec.begin());
      if (isLastProj) projOfs = projEnd;
    }
    else
      sim->extractPatchSolution(glbVec,pchVec,pch,ncmp,1);

    field.resize(ncmp,pchVec.size()/ncmp);
    field.fill(pchVec.ptr());
  };

  for (int idx = 1; idx <= sim->getNoPatches(); idx++) {
    int loc = sim->getLocalPatchIndex(idx);
    if ((!sim->getProcessAdm().dd.isPartitioned() || sim->getProcessAdm().getProcId() == 0) && loc > 0 &&
        (!(results & DataExporter::REDUNDANT) || sim->getGlobalProcessID() == 0)) // we own the patch
    {
      ASMbase* pch = sim->getPatch(loc);
      if (!pch) continue;

      if (results & DataExporter::PRIMARY && !sol->empty()) {
        Vector psol;
        size_t ndof1 = sim->extractPatchSolution(*sol,psol,pch,entry.second.ncmps,
                                                 usedescription ? 1 : 0);
        if (dynamic_cast<ASMsupel*>(pch))
        {
          // Hack for superelement patches: Expand to include the center node
          Matrix sField;
          if (pch->evalSolution(sField,psol,nullptr,0))
          {
            psol = sField;
            ndof1 = psol.size();
          }
        }
        const double* data = psol.ptr();
        if (usedescription)
          // Field assumed to be on basis 1 for now
          this->writeArray(group.front(), entry.second.description,
                           idx, ndof1, data, H5T_NATIVE_DOUBLE);
        else if (sim->mixedProblem())
          for (size_t b = 1; b <= sim->getNoBasis(); b++) {
            ndof1 = pch->getNoNodes(b)*pch->getNoFields(b);
            if (!prob->suppressOutput(10+b,ASM::PRIMARY))
              this->writeArray(group[b-1], prefix+prob->getField1Name(10+b),
                               idx, ndof1, data, H5T_NATIVE_DOUBLE);
            data += ndof1;
          }
        else
          this->writeArray(group.front(), prefix+prob->getField1Name(11),
                           idx, ndof1, data, H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::SECONDARY && !sol->empty()) {
        Matrix field;
        ASM::ResultClass rClass;
        hid_t gid = group.front();
        if (entry.second.description == "projected") {
          rClass = ASM::PROJECTED;
          // The projected solution has been registered as a separate field
          extractProjection(pch, *sol, field);
          if (sim->fieldProjections())
            gid = group.back();
        }
        else {
          rClass = ASM::SECONDARY;
          SIM::SolutionMode mode = prob->getMode();
          const_cast<SIMbase*>(sim)->setMode(SIM::RECOVERY);
          sim->extractPatchSolution({*sol},loc-1);
          sim->evalSecondarySolution(field,loc-1);
          const_cast<SIMbase*>(sim)->setMode(mode);
        }
        for (size_t j = 0; j < field.rows(); j++)
          if (!prob->suppressOutput(j,rClass))
            this->writeArray(gid, prefix+prob->getField2Name(j),
                             idx, field.cols(), field.getRow(j+1).ptr(),
                             H5T_NATIVE_DOUBLE);
      }

      if (proj) {
        hid_t gid = sim->fieldProjections() ? group.back() : group.front();
        for (size_t p = 0; p < proj->size(); p++)
          if (!proj->at(p).empty())
          {
            Matrix field;
            extractProjection(pch, proj->at(p), field, p+1 == proj->size());
            for (size_t j = 0; j < field.rows(); j++)
              if (!prob->suppressOutput(j,ASM::PROJECTED))
                this->writeArray(gid, prob->getField2Name(j,projPfx[p]),
                                 idx, field.cols(), field.getRow(j+1).ptr(),
                                 H5T_NATIVE_DOUBLE);
          }
      }

      if (norm) {
        Matrix patchEnorm;
        sim->extractPatchElmRes(*eNorm,patchEnorm,loc-1);
        for (size_t j = 1, l = 1; j <= normGrp && l < eNorm->rows(); j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              this->writeArray(egroup.front(),
                               prefix+norm->getName(j, k, j > 1 ? projPfx[j-2] : nullptr),
                               idx, patchEnorm.cols(), patchEnorm.getRow(l++).ptr(),
                               H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::EIGENMODES) {
        size_t iMode = 0;
        Vector psol;
        const std::vector<Mode>* modes = static_cast<const std::vector<Mode>*>(entry.second.data2.front());
        for (const Mode& mode : *modes)
        {
          std::stringstream str;
          str << level << '/' << sim->getName() << "-1/Eigenmode";
          hid_t gid;
          if (checkGroupExistence(m_file,str.str().c_str()))
            gid = H5Gopen2(m_file,str.str().c_str(),H5P_DEFAULT);
          else
            gid = H5Gcreate2(m_file,str.str().c_str(),0,H5P_DEFAULT,H5P_DEFAULT);

          std::stringstream str2;
          str2 << ++iMode;
          size_t ndof1 = sim->extractPatchSolution(mode.eigVec,psol,pch);
          this->writeArray(gid, str2.str(), idx, ndof1, psol.ptr(),
                           H5T_NATIVE_DOUBLE);
          if (sim->opt.eig == 3 || sim->opt.eig == 4 || sim->opt.eig == 6)
            str2 << "/Frequency";
          else
            str2 << "/Value";
          this->writeArray(gid, str2.str(), -1, 1, &mode.eigVal,
                           H5T_NATIVE_DOUBLE);
          if (idx == 1) {
            std::stringstream str3;
            str3 << iMode << "/eqn/";
            this->writeArray(gid, str3.str(), idx, mode.eqnVec.size(),
                             mode.eqnVec.ptr(), H5T_NATIVE_DOUBLE);
          }
          H5Gclose(gid);
        }
      }
    }
    else // must write empty dummy records for the other patches
    {
      double dummy=0.0;
      if (results & DataExporter::PRIMARY) {
        if (usedescription)
          writeArray(group.front(), entry.second.description,
                     idx, 0, &dummy, H5T_NATIVE_DOUBLE);
        else if (sim->mixedProblem())
          for (size_t b = 1; b <= sim->getNoBasis(); b++) {
            if (!prob->suppressOutput(10+b,ASM::PRIMARY))
              writeArray(group[b-1], prefix+prob->getField1Name(10+b),
                         idx, 0, &dummy, H5T_NATIVE_DOUBLE);
          }
        else
          writeArray(group.front(), prefix+prob->getField1Name(11),
                     idx, 0, &dummy, H5T_NATIVE_DOUBLE);
      }

      if (results & DataExporter::SECONDARY) {
        hid_t gid = group.front();
        ASM::ResultClass rClass = ASM::SECONDARY;
        if (entry.second.description == "projected") {
          rClass = ASM::PROJECTED;
          if (sim->fieldProjections())
            gid = group.back();
        }
        for (size_t j = 0; j < prob->getNoFields(2); j++)
          if (!prob->suppressOutput(j,rClass))
            writeArray(gid, prefix+prob->getField2Name(j),
                       idx, 0, &dummy,H5T_NATIVE_DOUBLE);
      }

      if (proj) {
        hid_t gid = sim->fieldProjections() ? group.back() : group.front();
        for (size_t p = 0; p < proj->size(); p++)
          if (!proj->at(p).empty())
            for (size_t j = 0; j < prob->getNoFields(2); j++)
              if (!prob->suppressOutput(j,ASM::PROJECTED))
                writeArray(gid, prob->getField2Name(j,projPfx[p]),
                           idx, 0, &dummy,H5T_NATIVE_DOUBLE);
      }

      if (norm)
        for (size_t j = 1; j <= normGrp; j++)
          for (size_t k = 1; k <= norm->getNoFields(j); k++)
            if (norm->hasElementContributions(j,k))
              writeArray(egroup.front(),
                         prefix+norm->getName(j, k, j > 1 ? projPfx[j-2] : nullptr),
                         idx, 0, &dummy,H5T_NATIVE_DOUBLE);

      if (results & DataExporter::EIGENMODES) // TODO (akva?)
        std::cerr <<"  ** HDF5Writer: Oops, eigenmodes not yet supported for distributed patches"
                  <<"\n     The HDF5-file will most likely be inconsistent."<< std::endl;
    }
  }

  delete norm;
  for (hid_t g : group)
    if (g != -1)
      H5Gclose(g);
  for (hid_t g : egroup)
    H5Gclose(g);
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
  for (int idx = 1; idx <= sim->getNoPatches(); idx++) {
    int loc = sim->getLocalPatchIndex(idx);
    if (loc > 0 && (sim->getProcessAdm().isParallel() ||
                    sim->getGlobalProcessID() == 0)) // we own the patch
    {
      Matrix patchEnorm;
      sim->extractPatchElmRes(infield,patchEnorm,loc-1);
      writeArray(group2,prefix+entry.second.description,idx,patchEnorm.cols(),
                 patchEnorm.getRow(1).ptr(),H5T_NATIVE_DOUBLE);
    }
    else { // must write empty dummy records for the other patches
      double dummy=0.0;
      writeArray(group2,prefix+entry.second.description,idx,0,&dummy,H5T_NATIVE_DOUBLE);
    }

  }
  H5Gclose(group2);
#else
  std::cout << "HDF5Writer: compiled without HDF5 support, no data written" << std::endl;
#endif
}


#ifdef HAS_HDF5
void HDF5Writer::writeBasis (const SIMbase* sim, const std::string& name,
                             int basis, int level, bool redundant, bool l2g)
{
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

  std::map<int, int> l2gNode;
  std::map<int, int> prevNode;
  if (l2g && sim->getNoBasis() > 1) {
    int iNode = 1;
    const ProcessAdm& adm = sim->getProcessAdm();
    bool parallel = sim->getProcessAdm().isParallel() && !adm.dd.isPartitioned();
    const std::vector<int>& MLGN = adm.dd.getMLGN();

#if defined(HAS_PETSC) || defined(HAVE_MPI)
    if (parallel && adm.getProcId() > 0) {
      int nNodes;
      adm.receive(nNodes, adm.getProcId()-1);
      std::vector<int> flatMap(nNodes);
      adm.receive(flatMap, adm.getProcId()-1);
      for (size_t n = 0; n < flatMap.size(); n += 2)
        prevNode.insert({flatMap[n], flatMap[n+1]});
      adm.receive(iNode, adm.getProcId() - 1);
    }
#endif

    for (int idx = 1; idx <= sim->getNoPatches(); idx++) {
      const ASMbase* pch = sim->getPatch(idx,true);
      if (!pch) continue;

      size_t ofs = 0;
      for (int b = 1; b < basis; ++b)
        ofs += pch->getNoNodes(b);
      for (size_t n = ofs; n < ofs + pch->getNoNodes(basis); ++n) {
        int node = pch->getNodeID(n+1);
        if (parallel) {
          auto prevIt = prevNode.find(MLGN[node-1]);
          if (prevIt != prevNode.end()) {
            l2gNode[node] = prevIt->second;
            continue;
          }
        }
        if (l2gNode.find(node) == l2gNode.end())
          l2gNode[node] = iNode++;
      }
    }

#if defined(HAS_PETSC) || defined(HAVE_MPI)
    if (parallel && adm.getProcId() < adm.getNoProcs() - 1) {
      std::vector<int> flatMap(2*l2gNode.size() + 2*prevNode.size());
      auto flatIt = flatMap.begin();
      for (const auto& it : l2gNode) {
        *flatIt++ = MLGN[it.first-1];
        *flatIt++ = it.second;
      }
      for (const auto& it : prevNode) {
        *flatIt++ = it.first;
        *flatIt++ = it.second;
      }
      int size = flatMap.size();
      adm.send(size, adm.getProcId() + 1);
      adm.send(flatMap, adm.getProcId() + 1);
      adm.send(iNode, adm.getProcId() + 1);
    }
#endif

  }

  for (int idx = 1; idx <= sim->getNoPatches(); idx++) {
    std::stringstream str;
    ASMbase* pch = sim->getPatch(idx,true);
    int alen = 0;
    if (pch && (!redundant || rank == 0))
    {
      pch->write(str,basis);
      alen = str.str().size();
    }
    this->writeArray(group, "basis", idx, alen , str.str().c_str(),
                     H5T_NATIVE_CHAR);

    if (l2g) {
      std::vector<int> nodeNumsLoc;
      const std::vector<int>* nodeNums = &nodeNumsLoc;
      if (pch) {
        const std::vector<int>& allNodeNums = pch->getGlobalNodeNums();
        if (l2gNode.empty())
          nodeNums = &allNodeNums;
        else {
          size_t start_loc = 0;
          for (int b = 1; b < basis; ++b)
            start_loc += pch->getNoNodes(b);

          nodeNumsLoc.resize(pch->getNoNodes(basis));
          for (size_t n = start_loc; n < start_loc + pch->getNoNodes(basis); ++n)
            nodeNumsLoc[n-start_loc] = l2gNode[allNodeNums[n]];
        }
      }
      this->writeArray(group, "l2g-node", idx, nodeNums->size(),
                       nodeNums->data(), H5T_NATIVE_INT);
    }
  }
  H5Gclose(group);
}
#endif


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


void HDF5Writer::writeNodalForces (int level, const DataEntry& entry)
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


bool HDF5Writer::writeLog (const std::string& data, const std::string& name)
{
#ifdef HAS_HDF5
  hid_t group;
  if (checkGroupExistence(m_file,"/log"))
    group = H5Gopen2(m_file,"/log",H5P_DEFAULT);
  else
    group = H5Gcreate2(m_file,"/log",0,H5P_DEFAULT,H5P_DEFAULT);

  int rank = 0;
  int size = 1;
#ifdef HAVE_MPI
  MPI_Comm_rank(*m_adm.getCommunicator(), &rank);
  MPI_Comm_size(*m_adm.getCommunicator(), &size);
#endif

  for (int i = 0; i < size; i++)
    writeArray(group, name.c_str(), i, rank == i ? data.size() : 0, data.data(), H5T_NATIVE_CHAR);
  H5Gclose(group);
  return true;

#else
  return false;
#endif
}
