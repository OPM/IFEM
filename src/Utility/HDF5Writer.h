// $Id$

#pragma once

#include "DataExporter.h"
#include "MatVec.h"

class SIMbase;


class HDF5Writer : public DataWriter
{
public:
  HDF5Writer(const std::string& name, bool append = false);
  virtual ~HDF5Writer() {}

  virtual int getLastTimeLevel();

  virtual void openFile(int level);
  virtual void closeFile(int level);

  virtual void writeVector(int level, const DataEntry& entry);
  virtual bool readVector(int level, const DataEntry& entry);
  virtual void writeSIM(int level, const DataEntry& entry);
  virtual bool readSIM(int level, const DataEntry& entry);

  bool readField(int level, const std::string& name,
                 Vector& vec, SIMbase* sim, int components);

protected:
  void writeArray(int group, const std::string& name,
                  int len, void* data);
  void readArray(int group, const std::string& name,
                 int& len, double*& data);
  bool checkGroupExistence(int parent, const char* group);

  int          m_file;
  unsigned int m_flag;
};
