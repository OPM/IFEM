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
#include "Field.h"
#include "ProcessAdm.h"
#include "Vec3.h"
#ifdef HAS_HDF5
#include "HDF5Writer.h"
#endif
#include <sstream>


FieldFunction::FieldFunction (const std::string& fileName,
                              const std::string& basisName,
                              const std::string& fieldName)
{
#ifdef HAS_HDF5
  HDF5Writer hdf5(fileName,ProcessAdm(),true,true);

  std::string g2;
  std::stringstream str;
  hdf5.readString("0/basis/"+basisName+"/1",g2);
  str << g2;

  pch = ASM2D::create(ASM::Spline);
  pch->read(str);

  Vector coefs;
  hdf5.readVector(0,fieldName,1,coefs);
  field = Field::create(pch,coefs);
#else
  std::cerr <<"WARNING: Compiled without HDF5 support,"
	    <<" field function is not instanciated."<< std::endl;
  field = NULL;
  pch = NULL;
#endif
}


FieldFunction::~FieldFunction ()
{
  delete field;
  delete pch;
}


Real FieldFunction::evaluate (const Vec3& X) const
{
  if (!field)
    return Real(0);

  const Vec4* x4 = dynamic_cast<const Vec4*>(&X);
  if (x4 && x4->idx > 0)
    return field->valueNode(x4->idx);

  return field->valueCoor(X);
}
