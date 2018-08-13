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
#include "HDF5Writer.h"
#include "ProcessAdm.h"
#include <sstream>
#endif


FieldFuncBase::FieldFuncBase (const std::string& fName) :
  hdf5(nullptr), pAdm(nullptr), pidx(0)
{
  lastLevel = 0;
  lastTime = 0.0;
#ifdef HAS_HDF5
  pAdm = new ProcessAdm();
  hdf5 = new HDF5Writer(fName,*pAdm,true,true);
  hdf5->readDouble(lastLevel,"timeinfo","SIMbase-1",lastTime);
#else
  std::cerr <<"WARNING: Compiled without HDF5 support,"
            <<" field function is not instantiated."<< std::endl;
#endif
}


FieldFuncBase::~FieldFuncBase ()
{
  for (ASMbase* pch : patch) delete pch;
#ifdef HAS_HDF5
  delete hdf5;
  delete pAdm;
#endif
}


int FieldFuncBase::findClosestLevel (double time) const
{
  if (time == lastTime) return lastLevel;
#ifdef HAS_HDF5
  double t;
  int incLev = time > lastTime ? 1 : -1;
  while (hdf5->readDouble(lastLevel+incLev,"timeinfo","SIMbase-1",t))
  {
    if (fabs(time-t) >= fabs(time-lastTime))
    {
#ifdef SP_DEBUG
      std::cout <<"FieldFuncBase: New time level "<< lastLevel
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


bool FieldFuncBase::load (const std::vector<std::string>& fieldNames,
                          const std::string& basisName, int level,
                          size_t nPatches, bool isScalar)
{
  size_t nOK = 0;
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
    if (hdf5->hasGeometries(level,basisName))
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
      sbasis << level <<"/basis/"<< basisName <<"/"<< ip+1;
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
        std::cerr <<" *** FieldFuncBase::load: Undefined basis "<< sbasis.str()
                  <<" ("<< g2.substr(0,9) <<")"<< std::endl;
    }

    if (patch[ip])
    {
      RealArrays coefs(nFldCmp);
      for (size_t i = 0; i < nFldCmp; i++)
      {
        hdf5->readVector(level,fieldNames[i],ip+1,coefs[i]);
#if SP_DEBUG > 1
        std::cout <<"FieldFuncBase::load: Reading \""<< fieldNames[i]
                  <<"\" ("<< coefs[i].size() <<") for patch "<< ip+1;
        for (size_t j = 0; j < coefs[i].size(); j++)
          std::cout << (j%10 ? ' ' : '\n') << coefs[i][j];
        std::cout << std::endl;
#endif
      }
      this->addPatchField(patch[ip],coefs);
      nOK++;
    }
    else
      std::cerr <<" *** FieldFuncBase::load: No field function created"
                <<" for patch "<< ip+1 << std::endl;
  }
#endif

  return nOK == nPatches;
}


FieldFunction::FieldFunction (const std::string& fileName,
                              const std::string& basisName,
                              const std::string& fieldName,
                              size_t nPatches, int level)
  : FieldFuncBase(fileName), currentLevel(level),
    fName(fieldName), bName(basisName)
{
  if (level >= 0)
    this->load({fieldName},basisName,level,nPatches,true);
}


void FieldFunction::clearField ()
{
  for (Field* f : field) delete f;
  field.clear();
}


void FieldFunction::addPatchField (ASMbase* pch, const RealArrays& coefs)
{
  field.push_back(Field::create(pch,coefs.front()));
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


FieldsFuncBase::FieldsFuncBase (const std::string& fileName,
                                const std::string& basisName,
                                const std::string& fieldName,
                                size_t nPatches, int level)
  : FieldFuncBase(fileName), currentLevel(level),
    fName(splitString(fieldName,[](int c){ return c == '|' ? 1 : 0; })),
    bName(basisName)
{
  if (level >= 0)
    this->load(fName,basisName,level,nPatches);
}


void FieldsFuncBase::clearField ()
{
  for (Fields* f : field) delete f;
  field.clear();
}


void FieldsFuncBase::addPatchField (ASMbase* pch, const RealArrays& coefs)
{
  if (coefs.size() == 1)
    field.push_back(Fields::create(pch,coefs.front(),1,pch->getNoFields(1)));
  else if (coefs.size() > 1)
  {
    RealArray coef;
    coef.reserve(coefs.size()*coefs.front().size());
    for (size_t i = 0; i < coefs.front().size(); i++)
      for (size_t j = 0; j < coefs.size(); j++)
        coef.push_back(coefs[j][i]);
    field.push_back(Fields::create(pch,coef,1,pch->getNoFields(1)));
  }
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
    if (level < 0) return vals;
    if (level != currentLevel)
      if (this->load(fName,bName,level))
        currentLevel = level;

    if (x4->idx > 0)
      field[pidx]->valueNode(x4->idx,vals);
    else
      field[pidx]->valueCoor(*x4,vals);
  }

  return vals;
}


VecFieldFunction::VecFieldFunction (const std::string& fileName,
                                    const std::string& basisName,
                                    const std::string& fieldName,
                                    size_t nPatches, int level)
  : FieldsFuncBase(fileName,basisName,fieldName,nPatches,level)
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
                                          size_t nPatches, int level)
  : FieldsFuncBase(fileName,basisName,fieldName,nPatches,level)
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
                                            size_t nPatches, int level)
  : FieldsFuncBase(fileName,basisName,fieldName,nPatches,level)
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
