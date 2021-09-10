// $Id$
//==============================================================================
//!
//! \file FieldFunctions.C
//!
//! \date Sep 20 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Field function implementations.
//!
//==============================================================================

#include "FieldFunctions.h"
#include "ASMbase.h"
#include "ASM2D.h"
#include "ASM3D.h"
#include "Field.h"
#include "Fields.h"
#include "Vec3.h"
#include "StringUtils.h"
#ifdef HAS_HDF5
#include "HDF5Reader.h"
#include "ProcessAdm.h"
#include <sstream>
#else
class HDF5Reader {};
class ProcessAdm {};
#endif


FieldFuncBase::~FieldFuncBase ()
{
  for (ASMbase* pch : patch) delete pch;
}


bool FieldFuncBase::setPatch (size_t pIdx)
{
  if (pIdx >= npch)
    return false;

  pidx = pIdx;
  return true;
}


FieldFuncHDF5::FieldFuncHDF5 (const std::string& fName)
  : hdf5(nullptr), pAdm(nullptr)
{
  lastLevel = 0;
  lastTime = 0.0;
#ifdef HAS_HDF5
  pAdm = new ProcessAdm();
  hdf5 = new HDF5Reader(fName,*pAdm);
#else
  std::cerr <<"WARNING: Compiled without HDF5 support,"
            <<" field function is not instantiated."<< std::endl;
#endif
}


FieldFuncHDF5::~FieldFuncHDF5 ()
{
  delete hdf5;
  delete pAdm;
}


int FieldFuncHDF5::findClosestLevel (double time) const
{
  if (time == lastTime) return lastLevel;
#ifdef HAS_HDF5
  double t;
  int incLev = time > lastTime ? 1 : -1;
  bool ok = true;
  while (ok)
  {
    std::stringstream str;
    str << lastLevel+incLev << "/timeinfo/SIMbase-1";
    ok = hdf5->readDouble(str.str(),t);
    if (ok && fabs(time-t) >= fabs(time-lastTime))
    {
#ifdef SP_DEBUG
      std::cout <<"FieldFuncHDF5: New time level "<< lastLevel
                <<" at t="<< lastTime <<" (dt="<< time-lastTime
                <<")"<< std::endl;
#endif
      return lastLevel; // lastLevel is the closest to time
    }
    lastTime = t;
    lastLevel += incLev;
  }
#endif
  return -1;
}


bool FieldFuncHDF5::load (const std::vector<std::string>& fieldNames,
                          const std::string& basisName, int level,
                          bool isScalar)
{
  size_t nOK = 0;
  size_t nPatches = 0;
#ifdef HAS_HDF5
  std::stringstream str;
  str << level << "/" << basisName << "/fields/" << fieldNames.front();
  nPatches = hdf5->getFieldSize(str.str());
#endif
  if (nPatches == 0)
    nPatches = patch.size();
  else if (patch.empty())
    patch.resize(nPatches,nullptr);
  this->clearField();

#ifdef HAS_HDF5
  size_t nFldCmp = fieldNames.size();
  size_t nFldC2D = isScalar ? 1 : (nFldCmp < 2 ? 2 : nFldCmp);
  size_t nFldC3D = isScalar ? 1 : (nFldCmp < 3 ? 3 : nFldCmp);
  for (size_t ip = 0; ip < nPatches; ip++)
  {
    std::stringstream str;
    str << level << "/" << basisName << "/basis";
    if (hdf5->getFieldSize(str.str()))
    {
      if (patch[ip])
      {
        // We have an updated basis at this level, replace current one
        if (patch[ip]->getNoParamDim() == 2) nFldC2D = patch[ip]->getNoFields();
        if (patch[ip]->getNoParamDim() == 3) nFldC3D = patch[ip]->getNoFields();
        delete patch[ip];
      }
      std::string g2;
      std::stringstream sbasis;
      sbasis << level << "/" << basisName << "/basis/"<< ip+1;
      hdf5->readString(sbasis.str(),g2);
      if (g2.compare(0,9,"200 1 0 0") == 0)
        patch[ip] = ASM2D::create(ASM::Spline,nFldC2D);
      else if (g2.compare(0,9,"700 1 0 0") == 0)
        patch[ip] = ASM3D::create(ASM::Spline,nFldC3D);
      else if (g2.compare(0,18,"# LRSPLINE SURFACE") == 0)
        patch[ip] = ASM2D::create(ASM::LRSpline,nFldC2D);
      else if (g2.compare(0,17,"# LRSPLINE VOLUME") == 0)
        patch[ip] = ASM3D::create(ASM::LRSpline,nFldC3D);
      else
        patch[ip] = nullptr;

      if (patch[ip])
      {
        std::stringstream strg2(g2);
        patch[ip]->read(strg2);
      }
      else
        std::cerr <<" *** FieldFuncHDF5::load: Undefined basis "<< sbasis.str()
                  <<" ("<< g2.substr(0,9) <<")"<< std::endl;
    }

    if (patch[ip])
    {
      std::vector<RealArray> coefs(nFldCmp);
      for (size_t i = 0; i < nFldCmp; i++)
      {
        std::stringstream str;
        str << level << "/" << basisName << "/fields/" << fieldNames[i] << "/" << ip+1;
        hdf5->readVector(str.str(),coefs[i]);
#if SP_DEBUG > 1
        std::cout <<"FieldFuncHDF5::load: Reading \""<< fieldNames[i]
                  <<"\" ("<< coefs[i].size() <<") for patch "<< ip+1;
        for (size_t j = 0; j < coefs[i].size(); j++)
          std::cout << (j%10 ? ' ' : '\n') << coefs[i][j];
        std::cout << std::endl;
#endif
      }
      if (nFldCmp > 1)
      {
        RealArray coef1;
        coef1.reserve(nFldCmp*coefs.front().size());
        for (size_t i = 0; i < coefs.front().size(); i++)
          for (size_t j = 0; j < nFldCmp; j++)
            coef1.push_back(coefs[j][i]);
        this->addPatchField(patch[ip],coef1);
      }
      else
        this->addPatchField(patch[ip],coefs.front());

      nOK++;
    }
    else
      std::cerr <<" *** FieldFuncHDF5::load: No field function created"
                <<" for patch "<< ip+1 << std::endl;
  }
#endif

  return nOK == nPatches;
}


FieldFunction::FieldFunction (const std::string& fileName,
                              const std::string& basisName,
                              const std::string& fieldName,
                              int level)
  : FieldFuncHDF5(fileName), currentLevel(level),
    fName(fieldName), bName(basisName)
{
  if (level >= 0)
    this->load({fieldName},basisName,level,true);
}


void FieldFunction::clearField ()
{
  for (Field* f : field) delete f;
  field.clear();
  npch = 0;
}


void FieldFunction::addPatchField (ASMbase* pch, const RealArray& coefs)
{
  field.push_back(Field::create(pch,coefs));
  npch = field.size();
}


Real FieldFunction::evaluate (const Vec3& X) const
{
  if (pidx >= field.size() || !field[pidx])
    return Real(0);

  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    return field[pidx]->valueCoor(X);

  int level = this->findClosestLevel(x4->t);
  if (level < 0) return Real(0);
  if (level != currentLevel)
    if (const_cast<FieldFunction*>(this)->load({fName},bName,level))
      currentLevel = level;

  if (x4->idx > 0)
    return field[pidx]->valueNode(x4->idx);
  else
    return field[pidx]->valueCoor(*x4);
}


FieldFuncStream::FieldFuncStream (const std::vector<ASMbase*>& patches,
                                  std::istream& istr)
{
  for (ASMbase* pch : patches)
  {
    RealArray coefs(pch->getNoNodes(1));
    for (size_t i = 0; i < coefs.size() && istr; i++) istr >> coefs[i];
    field.push_back(Field::create(pch,coefs,1,0));
  }
}


FieldFuncStream::~FieldFuncStream ()
{
  for (Field* f : field) delete f;
}


Real FieldFuncStream::evaluate (const Vec3& X) const
{
  if (pidx >= field.size() || !field[pidx])
    return Real(0);

  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    return field[pidx]->valueCoor(X);
  else if (x4->idx > 0)
    return field[pidx]->valueNode(x4->idx);
  else
    return field[pidx]->valueCoor(*x4);
}


FieldsFuncBase::FieldsFuncBase (const std::string& fileName,
                                const std::string& basisName,
                                const std::string& fieldName,
                                int level)
  : FieldFuncHDF5(fileName), currentLevel(level),
    fName(splitString(fieldName,[](int c){ return c == '|' ? 1 : 0; })),
    bName(basisName)
{
  if (level >= 0)
    this->load(fName,basisName,level);
}


void FieldsFuncBase::clearField ()
{
  for (Fields* f : field) delete f;
  field.clear();
  npch = 0;
}


void FieldsFuncBase::addPatchField (ASMbase* pch, const RealArray& coefs)
{
  field.push_back(Fields::create(pch,coefs,1,pch->getNoFields(1)));
  npch = field.size();
}


RealArray FieldsFuncBase::getValues (const Vec3& X)
{
  Vector vals;
  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    field[pidx]->valueCoor(X,vals);
  else
  {
    int level = this->findClosestLevel(x4->t);
    if (level < 0) return std::move(vals);
    if (level != currentLevel)
      if (this->load(fName,bName,level))
        currentLevel = level;

    if (x4->idx > 0)
      field[pidx]->valueNode(x4->idx,vals);
    else
      field[pidx]->valueCoor(*x4,vals);
  }

  return std::move(vals);
}


VecFieldFunction::VecFieldFunction (const std::string& fileName,
                                    const std::string& basisName,
                                    const std::string& fieldName,
                                    int level)
  : FieldsFuncBase(fileName,basisName,fieldName,level)
{
  if (!field.empty())
    ncmp = field.front()->getNoFields();
}


Vec3 VecFieldFunction::evaluate (const Vec3& X) const
{
  if (pidx >= field.size() || !field[pidx])
    return Vec3();

  return Vec3(const_cast<VecFieldFunction*>(this)->getValues(X).data(),ncmp);
}


TensorFieldFunction::TensorFieldFunction (const std::string& fileName,
                                          const std::string& basisName,
                                          const std::string& fieldName,
                                          int level)
  : FieldsFuncBase(fileName,basisName,fieldName,level)
{
  if (!field.empty())
    ncmp = field.front()->getNoFields();
}


Tensor TensorFieldFunction::evaluate (const Vec3& X) const
{
  if (pidx >= field.size() || !field[pidx])
    return Tensor(3);

  return const_cast<TensorFieldFunction*>(this)->getValues(X);
}


STensorFieldFunction::STensorFieldFunction (const std::string& fileName,
                                            const std::string& basisName,
                                            const std::string& fieldName,
                                            int level)
  : FieldsFuncBase(fileName,basisName,fieldName,level)
{
  if (!field.empty())
    ncmp = field.front()->getNoFields();
}


SymmTensor STensorFieldFunction::evaluate (const Vec3& X) const
{
  if (pidx >= field.size() || !field[pidx])
    return SymmTensor(3);

  return const_cast<STensorFieldFunction*>(this)->getValues(X);
}
