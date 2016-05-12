// $Id$
//==============================================================================
//!
//! \file SIMdummy.h
//!
//! \date May 27 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dummy SIM template class for drivers not associated with a FE model.
//!
//==============================================================================

#ifndef _SIM_DUMMY_H_
#define _SIM_DUMMY_H_

#include "SIMbase.h"


/*!
  \brief Template SIM class with some dummy implementations.
  \details This class only implements dummy versions for the pure virtual
  virtual methods of the base class SIMbase, and can be used as a base for
  simulator drivers that do not require any FE model.
*/

template<class Base> class SIMdummy : public Base
{
public:
  //! \brief Default constructor.
  SIMdummy(IntegrandBase* p = nullptr) : Base(p) {}
  //! \brief Empty destructor.
  virtual ~SIMdummy() {}
  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const { return 0; }
  //! \brief Reads a patch from given input stream.
  virtual ASMbase* readPatch(std::istream&,int, const SIMbase::CharVec&) const
  { return nullptr; }
protected:
  //! \brief Reads patches from given input stream.
  virtual bool readPatches(std::istream&,SIMdependency::PatchVec&,const char*) const
  { return false; }
  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  virtual bool addConstraint(int,int,int,int,int,int&,char)
  { return false; }
  //! \brief Creates a default geometry.
  virtual SIMdependency::PatchVec createDefaultGeometry(const TiXmlElement*) const
  { return SIMdependency::PatchVec(); }
  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints() {}
  //! \brief Creates topology for default geometry
  bool createDefaultTopology(const TiXmlElement*) { return false;}
  //! \brief Creates topology sets for default geometry
  TopologySet createDefaultTopologySets(const TiXmlElement*) const { return TopologySet();}
};

#endif
