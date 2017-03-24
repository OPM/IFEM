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

#ifndef _XML_WRITER_H
#define _XML_WRITER_H

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
  //! \brief A structure used when reading info from the file.
  struct Entry {
    std::string name;        //!< Field name
    std::string description; //!< Field description
    std::string basis;       //!< Name of the basis associated with the field
    std::string type;        //!< Field type
    int    patches;    //!< Number of patches in field
    int    components; //!< Number of components in field
    double timestep;   //!< The time step associated with the field
    int    interval;   //!< The dumping interval for the field
    bool   once;       //!< If true, field is only stored at first time level
  };

  //! \brief The constructor assigns the file name storing this field.
  //! \param[in] name The name (file name without extension) of data file
  XMLWriter(const std::string& name, const ProcessAdm& adm);

  //! \brief Empty destructor.
  virtual ~XMLWriter() {}

  //! \brief Returns the last time level stored in file.
  virtual int getLastTimeLevel();

  //! \brief Reads information from file.
  void readInfo();

  //! \brief Returns a const vector to the entries (\sa readInfo)
  const std::vector<Entry>& getEntries() const { return m_entry; }

  //! \brief Opens the file at a given time level.
  //! \param[in] level The requested time level
  virtual void openFile(int level);

  //! \brief Closes the file.
  //! \param[in] level Level we just wrote to the file
  virtual void closeFile(int level, bool = false);

  //! \brief Writes a vector to file.
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry);

  //! \brief Reads a vector from file.
  //! \param[in] level The time level to read the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual bool readVector(int level, const DataEntry& entry);

  //! \brief Writes data from a SIM to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The data entry describing the vector
  //! \param[in] prefix Prefix for field names
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool, const std::string& prefix);

  //! \brief Writes nodal forces to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeNodalForces(int level, const DataEntry& entry);

  //! \brief Writes knot span field to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the field
  //! \param[in] prefix Prefix for field
  virtual void writeKnotspan(int level, const DataEntry& entry,
                             const std::string& prefix);

  //! \brief Writes time stepping info to file.
  //! \param[in] level The time level to write the info at
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int interval, const TimeStep& tp);

  //! \brief Write a basis to file
  //! \param[in] level The time level to write the basis at
  //! \param[in] entry The DataEntry describing the basis
  //! \param[in] prefix Prefix for basis
  virtual void writeBasis(int level, const DataEntry& entry,
                          const std::string& prefix) {}

  //! \brief No XML for restart data
  //! \param level Level to write data at
  //! \param data Data to write
  virtual bool writeRestartData(int level,
                                const DataExporter::SerializeData& data)
  { return true; }
protected:
  //! \brief Internal helper function adding an XML-element to the file
  //! \param[in] name The name of the field to add
  //! \param[in] description The description of the field to add
  //! \param[in] geometry The name of the geometry the field uses
  //! \param[in] components Number of components in the field
  //! \param[in] patches Number of patches in the field
  //! \param[in] type The type of the field
  //! \param[in] once If \e true, write field only once
  void addField(const std::string& name, const std::string& description,
                const std::string& geometry, int components, int patches,
                const std::string& type = "field", bool once = false);

private:
  std::vector<Entry> m_entry; //!< A vector of entries read from file

  TiXmlDocument* m_doc;  //!< Our XML document
  TiXmlNode*     m_node; //!< The document's root node

  double m_dt;       //!< The current (constant) time step size
  int    m_interval; //!< Stride for dumping (dump at every m_interval'th level)
};

#endif
