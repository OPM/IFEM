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

  //! \brief Struct holding information about an initial condition.
  struct ICInfo
  {
    int sim_level;  //!< The time level for the field in the SIM class
    int file_level; //!< The time level for the field in the file
    int geo_level;  //!< The time level for the geometry in the file
    char basis;     //!< The basis to inject field into (for mixed)
    int component;  //!< Component for field (for functions)
    std::string sim_field;  //!< The name of the field in the SIM class
    std::string file_field; //!< The field name in the file or type of function
    std::string function;   //!< Function if given in function form
    //! \brief Default constructor.
    ICInfo() : sim_level(0), file_level(0), geo_level(0) {}
    //! \brief Constructor providing the field name.
    ICInfo(const std::string& f) : sim_level(0), file_level(0), geo_level(0),
                                   basis(1), component(-1),
                                   sim_field(f), file_field(f) {}
  };

  //! \brief Initial condition container
  typedef std::map< std::string,std::vector<ICInfo> > InitialCondMap;

private:
  //! \brief Struct holding information about an inter-SIM dependency.
  struct Dependency
  {
    SIMdependency* sim;            //!< SIM object holding the dependent field
    PatchVec       patches;        //!< Patch geometry the field is defined over
    std::string    name;           //!< Field name
    short int      components;     //!< Number of field components per node
    char           differentBasis; //!< Toggle usage of an independent basis
    //! \brief Default constructor.
    Dependency(SIMdependency* s = nullptr, const std::string& f = "",
               short int n = 1) : sim(s), name(f), components(n),
                                  differentBasis(0) {}
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

  //! \brief Returns a const reference to the initial condition container.
  const InitialCondMap& getICs() const { return myICs; }

  //! \brief Checks whether a named initial condition is present.
  virtual bool hasIC(const std::string& name) const;

  //! \brief Returns the number of spatial dimensions in the model.
  virtual size_t getNoSpaceDim() const = 0;
  //! \brief Returns the given name of this simulator.
  virtual std::string getName() const = 0;

  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] patches The geometry the field is defined over
  //! \param[in] diffBasis If non-null, use diffBasis base from patch vector
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc, const PatchVec& patches,
                                  char diffBasis = 0);
  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc = 1);

  //! \brief Initializes the nodal vector of named field in this SIM.
  bool fillField(const std::string& name, const std::vector<double>& values);
  //! \brief Returns the nodal vector of named field in this SIM.
  virtual utl::vector<double>* getField(const std::string& name);
  //! \brief Returns the nodal vector of named field in this SIM.
  virtual const utl::vector<double>* getField(const std::string& name) const;
  //! \brief Returns the nodal vector of named field in a dependent SIM.
  const utl::vector<double>* getDependentField(const std::string& name) const;
  //! \brief Returns a spline patch associated with a dependent field.
  ASMbase* getDependentPatch(const std::string& name, int pindx) const;
  //! \brief Registers a named field with associated nodal vector in this SIM.
  void registerField(const std::string& name, const utl::vector<double>& vec);

private:
  //! \brief Returns an iterator pointing to a named dependency.
  DepVector::const_iterator getDependency(const std::string& name) const;

protected:
  //! \brief Extracts local solution vector(s) for all dependent fields.
  //! \param problem Object with problem-specific data and methods
  //! \param[in] model Patch geometry of this SIM object
  //! \param[in] pindx Local patch index to extract solution vectors for
  bool extractPatchDependencies(IntegrandBase* problem,
                                const PatchVec& model, size_t pindx) const;

private:
  FieldMap  myFields;  //!< The named fields of this SIM object
  DepVector depFields; //!< Other fields this SIM objecy depends on

protected:
  InitialCondMap myICs; //!< The initial conditions
};

#endif
