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
#include "ProcessAdm.h"
#include "Vec3.h"
#include "StringUtils.h"
#ifdef HAS_HDF5
#include "HDF5Writer.h"
#endif
#include <sstream>


FieldFunction::FieldFunction (const std::string& fileName,
                              const std::string& basisName,
                              const std::string& fieldName,
                              size_t npatches, int level) : pidx(0)
{
#ifdef HAS_HDF5
  ProcessAdm adm;
  HDF5Writer hdf5(fileName,adm,true,true);

  for (size_t p = 0; p < npatches; ++p) {
    std::string g2;
    std::stringstream str;
    std::stringstream str2;
    str2 << level << "/basis/" << basisName << "/" << p+1;
    hdf5.readString(str2.str(),g2);
    str << g2;
    std::string head = g2.substr(0,9);
    std::string head2 = g2.substr(0,12);
    ASMbase* pch = nullptr;
    if (head == "200 1 0 0")
      pch = ASM2D::create(ASM::Spline);
    else if (head == "700 1 0 0")
      pch = ASM3D::create(ASM::Spline);
    else if (head2 == "# LRSPLINE S")
      pch = ASM2D::create(ASM::LRSpline);
    else if (head2 == "# LRSPLINE V")
      pch = ASM3D::create(ASM::LRSpline);

    if (pch)
    {
      pch->read(str);
      Vector coefs;
      hdf5.readVector(level,fieldName,p+1,coefs);
      field.push_back(Field::create(pch,coefs));
      patch.push_back(pch);
    }
    else {
      std::cerr <<" *** FieldFunction: Unknown basis type, no function created."
                << std::endl;
      return;
    }
  }
#else
  std::cerr <<"WARNING: Compiled without HDF5 support,"
            <<" field function is not instantiated."<< std::endl;
#endif
}


FieldFunction::~FieldFunction ()
{
  for (Field* f : field)
    delete f;
  for (ASMbase* pch : patch)
    delete pch;
}


Real FieldFunction::evaluate (const Vec3& X) const
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
                                size_t npatches, int level) : pidx(0)
{
#ifdef HAS_HDF5
  HDF5Writer hdf5(fileName,ProcessAdm(),true,true);

  std::vector<std::string> fieldNames = splitString(fieldName,
                                                 [](int c) { return c == '|' ? 1 : 0; });

  for (size_t p = 0; p < npatches; ++p) {
    std::string g2;
    std::stringstream str;
    std::stringstream str2;
    str2 << level << "/basis/" << basisName << "/" << p+1;
    hdf5.readString(str2.str(),g2);
    str << g2;
    std::string head = g2.substr(0,9);
    std::string head2 = g2.substr(0,12);
    ASMbase* pch = nullptr;
    if (head == "200 1 0 0")
      pch = ASM2D::create(ASM::Spline, std::max(fieldNames.size(), 2LU));
    else if (head == "700 1 0 0")
      pch = ASM3D::create(ASM::Spline, std::max(fieldNames.size(), 3LU));
    else if (head2 == "# LRSPLINE S")
      pch = ASM2D::create(ASM::LRSpline, std::max(fieldNames.size(), 2LU));
    else if (head2 == "# LRSPLINE V")
      pch = ASM3D::create(ASM::LRSpline, std::max(fieldNames.size(), 3LU));

    if (pch)
    {
      pch->read(str);
      Vector coefs;
      if (fieldNames.size() == 1)
        hdf5.readVector(level,fieldName,p+1,coefs);
      else {
        Vectors coef(fieldNames.size());
        for (size_t i = 0; i < fieldNames.size(); ++i)
          hdf5.readVector(level,fieldNames[i],p+1,coef[i]);
        coefs.reserve(coef.size()*coef.front().size());
        for (size_t i = 0; i < coef.front().size(); ++i)
          for (size_t j = 0; j < coef.size(); ++j)
            coefs.push_back(coef[j][i]);
      }
      field.push_back(Fields::create(pch,coefs,1,pch->getNoFields(1)));
      patch.push_back(pch);
    }
    else {
      std::cerr <<" *** FieldsFuncBase: Unknown basis type, no function created."
                << std::endl;
      return;
    }
  }
#else
  std::cerr <<"WARNING: Compiled without HDF5 support,"
            <<" vector field function is not instantiated."<< std::endl;
#endif
}


FieldsFuncBase::~FieldsFuncBase ()
{
  for (Fields* f : field)
    delete f;
  for (ASMbase* pch : patch)
    delete pch;
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

  Vector vals;
  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    field[pidx]->valueCoor(X, vals);
  else if (x4->idx > 0)
    field[pidx]->valueNode(x4->idx, vals);
  else
    field[pidx]->valueCoor(*x4,vals);

  return Vec3(vals.ptr(), ncmp);
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

  Vector vals;
  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    field[pidx]->valueCoor(X, vals);
  else if (x4->idx > 0)
    field[pidx]->valueNode(x4->idx, vals);
  else
    field[pidx]->valueCoor(*x4, vals);

  return vals;
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

  Vector vals;
  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (!x4)
    field[pidx]->valueCoor(X, vals);
  else if (x4->idx > 0)
    field[pidx]->valueNode(x4->idx, vals);
  else
    field[pidx]->valueCoor(*x4, vals);

  return vals;
}
