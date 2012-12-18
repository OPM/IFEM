// $Id$
//==============================================================================
//!
//! \file XMLWriter.h
//!
//! \date May 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of metadata associated with HDF5 to XML.
//!
//==============================================================================

#pragma once

#include "DataExporter.h"

class TiXmlDocument;
class TiXmlNode;
class TimeStep;


/*!
  \brief Write data (metadata) to a XML file.

  \details The XML writer writes metadata (name of fields, description,...)
  in a humanly readable (XML) text format.
*/

class XMLWriter : public DataWriter
{
public:
  //! \brief A structure used when reading info from the file
  struct Entry {
    //! \brief Name of field
    std::string name;
    //! \brief Description of field
    std::string description;
    //! \brief The name of the basis the field is associated with
    std::string basis;
    //! \brief Number of patches in field
    int patches;
    //! \brief Number of components in field
    int components;
    //! \brief The timestep associated with the field
    double timestep;
    //! \brief The temporal order associated with the field
    int order;
    //! \brief The dumping interval for the field
    int interval;
    //! \brief The type of the field
    std::string type;
  };

  //! \brief Default constructor
  //! \param[in] name The name (filename without extension) of data file
  XMLWriter(const std::string& name);

  //! \brief Empty destructor
  virtual ~XMLWriter() {}

  //! \brief Return the last time level stored in file
  virtual int getLastTimeLevel();

  //! \brief Read info from file
  void readInfo();

  //! \brief Returns a const vector to the entries (\sa readInfo)
  const std::vector<Entry>& getEntries() const { return m_entry; }

  //! \brief Open the file at a given time level
  //! \param[in] level The requested time level
  virtual void openFile(int level);

  //! \brief Close the file
  //! \param[in] level Level we just wrote to the file
  //! \param[in] force Ignored
  virtual void closeFile(int level, bool force=false);

  //! \brief Write a vector to file
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry);

  //! \brief Read a vector from file
  //! \param[in] level The time level to read the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual bool readVector(int level, const DataEntry& entry);

  //! \brief Write data from a SIM to file
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  //! \param[in] geometryUpdated Ignored
  //! \param[in] prefix Prefix for field names
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated, const std::string& prefix);

  //! \brief Read data from a file into SIM
  //! \param[in] level The time level to read the data at
  //! \param[in] entry The DataEntry describing the SIM
  virtual bool readSIM(int level, const DataEntry& entry);

  //! \brief Write time stepping info to file (currently a dummy)
  //! \param[in] level The time level to write the info at
  //! \param[in] order The temporal order
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int order, int interval,
                             const TimeStep& tp);
protected:
  //! \brief Internal helper function
  //! \param[in] name The name of the field to add
  //! \param[in] description The description of the field to add
  //! \param[in] geometry The name of the geometry the field uses
  //! \param[in] components Number of components in the field
  //! \param[in] patches Number of patches in the field
  //! \param[in] type The type of the field
  void addField(const std::string& name, const std::string& description,
                const std::string& geometry, int components, int patches,
                const std::string& type="field");

  //! \brief A vector of entries read from file
  std::vector<Entry> m_entry;

  //! \brief Our xml document
  TiXmlDocument* m_doc;
  //! \brief The document's root node
  TiXmlNode* m_node;
  //! \brief The current (constant) timestep
  double m_dt;
  //! \brief The temporal order
  int m_order;
  //! \brief A stride for dumping. We dump at every m_interval'th time level
  int m_interval;
};
