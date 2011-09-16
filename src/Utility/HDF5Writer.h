// $Id$

#pragma once

#include "DataExporter.h"
#include "MatVec.h"

class SIMbase;


class HDF5Writer : public DataWriter
{
public:
  HDF5Writer(const std::string& name, bool append = false, bool keepopen=false);
  virtual ~HDF5Writer() {}

  virtual int getLastTimeLevel();

  virtual void openFile(int level);
  virtual void closeFile(int level, bool force=false);

  virtual void writeVector(int level, const DataEntry& entry);
  virtual bool readVector(int level, const DataEntry& entry);
  virtual void writeSIM(int level, const DataEntry& entry);
  virtual bool readSIM(int level, const DataEntry& entry);
  virtual bool writeTimeInfo(int level, int order, int interval,
                             SIMparameters& tp);

  bool readField(int level, const std::string& name,
                 Vector& vec, SIMbase* sim, int components);
  void readString(const std::string& name, std::string& out);
  bool readVector(int level, const std::string& name,
                  int patch, Vector& vec);

protected:
  void writeArray(int group, const std::string& name,
                  int len, const void* data, int type);
  void writeBasis(SIMbase* SIM, const std::string& name,
                     int basis, int level);
  void readArray(int group, const std::string& name,
                 int& len, double*& data);
  bool checkGroupExistence(int parent, const char* group);

  int          m_file;
  unsigned int m_flag;
  bool m_keepOpen;
};
