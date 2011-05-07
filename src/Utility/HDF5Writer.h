// $Id$

#pragma once

#include "DataExporter.h"


class HDF5Writer : public DataWriter {
  public:

    HDF5Writer(const std::string& name);
    virtual ~HDF5Writer();

    int getLastTimeLevel();

    void openFile(int level);
    void closeFile(int level);

    void writeVector(int level, const DataEntry& entry);
    void readVector(int level, const DataEntry& entry);
    void writeSIM(int level, const DataEntry& entry);
    void readSIM(int level, const DataEntry& entry);

  protected:
    void writeArray(int group, const std::string& name,
                    int len, void* data);
    void readArray(int group, const std::string& name,
                   int& len, double*& data);
    bool checkGroupExistence(int parent, const char* group);

    std::string m_hdf5;
    int m_file;
    unsigned int m_flag;

    int m_rank; // MPI rank
    int m_size; // number of MPI nodes
};
