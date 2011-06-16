#pragma once

// $Id$
//==============================================================================
//!
//! \file VTU.h
//!
//! \date Jun 14 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basic VTU file writer class
//!
//==============================================================================


#include "MatVec.h"
#include "VTF.h"
#include <string>


class ElementBlock;


class VTU : public VTF {
  public:
    VTU(const char* base, bool single);
    virtual ~VTU();

    bool writeGrid(const ElementBlock* lvb, const char* name);

    bool writeVres(const std::vector<double>& field, int blockID,
                   int geomID, int components);
    bool writeNres(const std::vector<double>& vec, int blockID, int geomID);

    bool writeVblk(const std::vector<int>& vBlockIDs, 
                   const char* resultName, int idBlock, int iStep=1);

    bool writeSblk(const std::vector<int>& sBlockIDs,
                   const char* resultName, int idBlock, int iStep=1);

    bool writeState(int iStep, const char* fmt, real refValue, int refType=0);
  protected:
    std::string m_base;
    std::vector<const ElementBlock*> m_geom;
    struct FieldInfo {
      Vector* data;
      int components;
      int patch;
      std::string name;
    };
    std::map<int,FieldInfo> m_field;
    bool m_single;
};
