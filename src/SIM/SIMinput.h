// $Id$
//==============================================================================
//!
//! \file SIMinput.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for simulators with input parsing functionality.
//!
//==============================================================================

#ifndef _SIM_INPUT_H
#define _SIM_INPUT_H

#include "MatVec.h"
#include <iostream>
#include <string>
#include <map>

class TiXmlElement;
class ASMbase;


/*!
  \brief Base class for NURBS-based FEM simulators with input file parsing.
*/

class SIMinput
{
protected:
  //! \brief The default constructor initializes \a myPid and \a nProc.
  SIMinput();

public:
  //! \brief Empty destructor.
  virtual ~SIMinput() {}

  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);

private:
  //! \brief Reads a flat text input file.
  bool readFlat(const char* fileName);
  //! \brief Reads an XML input file.
  bool readXML(const char* fileName);

protected:
  //! \brief Handles the parsing order for certain XML-tags.
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \param[out] parsed Vector of XML-elements that was parsed
  //!
  //! \details Certain tags need to be parsed before others. This method takes
  //! care of this. It is called by the \a readXML method in order to read the
  //! top level tags in the required order. It can also be called by the
  //! application-specific SIM class prior to parsing its data blocks.
  //! In that case the \a getPrioritizedTags method should be reimplemented
  //! by the sub-class to take care of the application-specific tags.
  bool handlePriorityTags(const TiXmlElement* base,
			  std::vector<const TiXmlElement*>& parsed);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const { return NULL; }

public:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is) = 0;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) = 0;

  static int msgLevel; //!< Controls the amount of console output during solving

  //! \brief Obtain a vector pointing to a named field in this SIM
  //! \param[in] name Name of field
  //! \return Pointer to the vector holding the field, if any
  const Vector* getNamedField(const std::string& name);

  //! \brief Register a dependency on a field from another SIM
  //! \param[in] sim The SIM holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] patches The geometry the field is defined over
  void registerDependency(SIMinput* sim, const std::string& name,
                          short int nvc, const std::vector<ASMbase*>* patches=0);

protected:
  //! \brief Register a named field in this integrand
  void registerNamedField(const std::string& name, const Vector* vec);
  int myPid; //!< Processor ID in parallel simulations
  int nProc; //!< Number of processors in parallel simulations

  //! \brief Struct holding information about a inter-SIM dependency
  struct Dependency {
    SIMinput* sim;
    std::string name;
    short int components;
    std::vector<ASMbase*> patches;
  };
  typedef std::vector<Dependency>             DepVector;
  typedef std::map<std::string,const Vector*> FieldMap;

  FieldMap  myFields;  //!< The fields in this SIM
  DepVector depFields; //!< Fields we depend on
};

#endif
