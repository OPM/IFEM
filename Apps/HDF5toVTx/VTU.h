// $Id$
//==============================================================================
//!
//! \file VTU.h
//!
//! \date Jun 14 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basic VTU file writer class.
//!
//==============================================================================

#pragma once

#include "VTF.h"
#include <string>
#include <map>


/*!
  \brief Basic VTU file writer class.
*/

class VTU : public VTF {
  public:
    VTU(const char* base, bool single);
    virtual ~VTU();

    void clearGeometryBlocks();

    bool writeGrid(const ElementBlock* lvb, const char* name, int iStep=1);

    using VTF::writeVres;
    bool writeVres(const std::vector<Real>& field, int blockID,
                   int geomID, int components);
    bool writeNres(const std::vector<Real>& vec, int blockID, int geomID);
    bool writeEres(const std::vector<Real>& vec, int blockID, int geomID);

    bool writeVblk(const std::vector<int>& vBlockIDs,
                   const char* resultName, int idBlock, int iStep=1);
    bool writeDblk(const std::vector<int>& vBlockIDs,
                   const char* resultName, int idBlock, int iStep=1);
    bool writeSblk(const std::vector<int>& sBlockIDs,
                   const char* resultName, int idBlock, int iStep=1,
                   bool elementData=false);

    bool writeState(int iStep, const char* fmt, Real refValue, int refType=0);

    const ElementBlock* getBlock(int geomID) const { return m_geom[geomID-1]; }

  protected:
    std::string m_base;
    std::vector<const ElementBlock*> m_geom;
    struct FieldInfo {
      std::vector<Real>* data;
      int components;
      int patch;
      bool cellData;
      std::string name;
    };
    std::map<int,FieldInfo> m_field;
    bool m_single;
};
