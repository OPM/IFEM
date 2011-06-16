// $Id$
//==============================================================================
//!
//! \file VTU.C
//!
//! \date Jun 14 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basic VTU file writer class
//!
//==============================================================================


#include "ElementBlock.h"
#include "VTU.h"

#include <fstream>
#include <iomanip>
#include <sstream>


VTU::VTU(const char* base, int format)
  : VTF(NULL,format), m_base(base)
{
  m_base = m_base.substr(0,m_base.rfind('.'));
}


VTU::~VTU()
{
  for (size_t i=0;i<m_geom.size();++i)
    delete m_geom[i];
  m_geom.clear();
}


bool VTU::writeGrid(const ElementBlock* block, const char* name)
{
  m_geom.push_back(block);
  return true;
}


bool VTU::writeVres(const std::vector<real>& field, int blockID, int geomID,
                    int components)
{
  m_field[blockID].data = new Vector(field.data(),field.size());
  m_field[blockID].components = components;
  m_field[blockID].patch = geomID;
  return true;
}


bool VTU::writeNres(const std::vector<real>& field, int blockID, int geomID)
{
  m_field[blockID].data = new Vector(field.data(),field.size());
  m_field[blockID].components = 1;
  m_field[blockID].patch = geomID;
  return true;
}


bool VTU::writeVblk(const std::vector<int>& vBlockIDs, 
                    const char* resultName, int idBlock, int iStep)
{
  for (size_t i=0;i<vBlockIDs.size();++i)
    m_field[vBlockIDs[i]].name = resultName;
  return true;
}


bool VTU::writeSblk(const std::vector<int>& sBlockIDs, 
                    const char* resultName, int idBlock, int iStep)
{
  for (size_t i=0;i<sBlockIDs.size();++i)
    m_field[sBlockIDs[i]].name = resultName;
  return true;
}


bool VTU::writeState(int iStep, const char* fmt, real refValue, int refType)
{
  std::ofstream file;
  std::stringstream str;
  str << m_base << "-" << std::setfill('0') << std::setw(5) << iStep-1 << ".vtu";
  file.open(str.str().c_str());

  file << "<?xml version=\"1.0\"?>" << std::endl;
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  file << "\t<UnstructuredGrid>" << std::endl;
  for (size_t i=0;i<m_geom.size();++i) {
    file << "\t\t<Piece NumberOfCells=\"" << m_geom[i]->getNoElms()
         << "\" NumberOfPoints=\"" << m_geom[i]->getNoNodes()
         << "\">" << std::endl;

    // dump geometry
    file << "\t\t\t<Points>" << std::endl;
    file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates\""
         << " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    file << "\t\t\t\t\t";
    for (std::vector<Vec3>::const_iterator it  = m_geom[i]->begin_XYZ();
                                           it != m_geom[i]->end_XYZ();++it) {
      it->print(file);
      file << " ";
    }
    file << std::endl << "\t\t\t\t</DataArray>" << std::endl;
    file << "\t\t\t</Points>" << std::endl;

    file << "\t\t\t<Cells>" << std::endl;
    file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\""
         << " NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    file << "\t\t\t\t\t";
    for (size_t j=0;j<m_geom[i]->getNoElms()*m_geom[i]->getNoElmNodes();++j)
      file << m_geom[i]->getElements()[j] << " ";
    file << std::endl << "\t\t\t\t</DataArray>" << std::endl;

    file << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\""
         << " NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    file << "\t\t\t\t\t";
    std::string type = "9";
    if (m_geom[i]->getNoElmNodes() == 8)
      type="12";
    for (size_t k=0;k<m_geom[i]->getNoElms();++k)
      file << type << " ";
    file << std::endl << "\t\t\t\t</DataArray>" << std::endl;

    file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\""
         << " NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    file << "\t\t\t\t\t";
    for (size_t k=0;k<m_geom[i]->getNoElms();++k)
      file << (k+1)*m_geom[i]->getNoElmNodes() << " ";
    file << std::endl << "\t\t\t\t</DataArray>" << std::endl;
    file << "\t\t\t</Cells>" << std::endl;

    // now add point datas
    file << "\t\t\t<PointData Scalars=\"scalars\">" << std::endl;
    for (std::map<int,FieldInfo>::iterator it2  = m_field.begin();
                                           it2 != m_field.end();++it2) {
      if (it2->second.patch == (int)i+1) {
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"" << it2->second.name << "\""
             << " NumberOfComponents=\"" << it2->second.components
             << "\" format=\"ascii\">" << std::endl;
        file << "\t\t\t\t\t";
        for (size_t k=0;k<it2->second.data->size();++k)
          file << (*it2->second.data)[k] << " ";
        delete it2->second.data;
        file << std::endl << "\t\t\t\t</DataArray>" << std::endl;
      }
    }
    file << "\t\t\t</PointData>" << std::endl;
    file << "\t\t</Piece>" << std::endl;
  }
  file << "\t</UnstructuredGrid>" << std::endl;
  file << "</VTKFile>" << std::endl;
  file.close();
  m_field.clear();
  return true;
}
