// $Id$
//==============================================================================
//!
//! \file HDF5WriterSSV.C
//!
//! \date Apr 27 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file for SSV.
//!
//==============================================================================

#include "HDF5WriterSSV.h"

#include "ASMbase.h"
#include "IntegrandBase.h"
#include "MatVec.h"
#include "ProcessAdm.h"
#include "SIMbase.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif


HDF5WriterSSV::HDF5WriterSSV (const std::string& name, const ProcessAdm& adm)
  : HDF5Writer(name+"_ssv",adm,false)
{
}

HDF5WriterSSV::~HDF5WriterSSV()
{
  if (m_file != -1)
    H5Fclose(m_file);
}


int HDF5WriterSSV::getLastTimeLevel ()
{
  return 0;
}


void HDF5WriterSSV::openFile(int level)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    return;

  if (!HDF5Base::openFile(m_flag, true))
    return;

  if (!checkGroupExistence(m_file,"/0"))
    H5Gclose(H5Gcreate2(m_file,"/0",0,H5P_DEFAULT,H5P_DEFAULT));
#endif
}


void HDF5WriterSSV::closeFile(int level)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    H5Fflush(m_file,H5F_SCOPE_GLOBAL);
#endif
}


bool HDF5WriterSSV::prepare (int level, const DataEntry& entry,
                             bool geometryUpdated, const std::string& prefix)
{
  if (geometryUpdated) {
    H5Fclose(m_file);
    m_file = -1;
    openFile(0);
    level = 0;
  }

  DataEntry entry2(entry);
  entry2.second.results = DataExporter::PRIMARY;

  if (level == 0)
    this->writeSIMInt(0, entry2, geometryUpdated, prefix, true);

  return true;
}


void HDF5WriterSSV::writeSIM (int level, const DataEntry& entry,
                              bool geometryUpdated, const std::string& prefix)
{
  H5Fstart_swmr_write(m_file);
  DataEntry entry2(entry);
  entry2.second.results = DataExporter::PRIMARY;
  return this->HDF5Writer::writeSIMInt(0, entry2, geometryUpdated, prefix, false);
}
