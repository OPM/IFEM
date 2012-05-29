// $Id$
//==============================================================================
//!
//! \file SIMdependency.h
//!
//! \date May 22 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration of simulators with dependencies to other simulators.
//!
//==============================================================================

#ifndef _SIM_DEPENDENCY_H
#define _SIM_DEPENDENCY_H

#include "matrix.h"
#include <string>
#include <map>

class ASMbase;
class IntegrandBase;


/*!
  \brief Class administering inter-SIM field dependencies.
*/

class SIMdependency
{
public:
  //! \brief Spline patch container
  typedef std::vector<ASMbase*> PatchVec;

private:
  //! \brief Struct holding information about an inter-SIM dependency.
  struct Dependency
  {
    SIMdependency* sim;            //!< SIM object holding the dependent field
    PatchVec       patches;        //!< Patch geometry the field is defined over
    std::string    name;           //!< Field name
    short int      components;     //!< Number of field components per node
    bool           differentBasis; //!< Toggle usage of an independent basis
    //! \brief Default constructor.
    Dependency() : sim(NULL), components(1), differentBasis(false) {}
  };

  //! \brief SIM dependency container
  typedef std::vector<Dependency> DepVector;
  //! \brief Field name to nodal values map
  typedef std::map<std::string,const utl::vector<double>*> FieldMap;

protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  SIMdependency() {}

public:
  //! \brief Empty destructor.
  virtual ~SIMdependency() {}

  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] patches The geometry the field is defined over
  //! \param[in] diffBasis Different basis for the SIM class and the field
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc, const PatchVec& patches,
			  bool diffBasis = false);
  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc = 1);

  //! \brief Returns the nodal vector of named field in this SIM.
  const utl::vector<double>* getField(const std::string& name);

protected:
  //! \brief Registers a named field with associated nodal vector in this SIM.
  void registerField(const std::string& name, const utl::vector<double>& vec);

protected:
  //! \brief Extracts local solution vector(s) for all dependent fields.
  //! \param problem Object with problem-specific data and methods
  //! \param[in] model Patch geometry of this SIM object
  //! \param[in] pindx Local patch index to extract solution vectors for
  bool extractPatchDependencies(IntegrandBase* problem,
                                const PatchVec& model, size_t pindx);

private:
  FieldMap  myFields;  //!< The named fields of this SIM object
  DepVector depFields; //!< Other fields this SIM objecy depends on
};

#endif
