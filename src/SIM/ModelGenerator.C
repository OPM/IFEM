// $Id$
//==============================================================================
//!
//! \file ModelGenerator.C
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for model generators for NURBS-based FEM simulators.
//!
//==============================================================================


#include "ModelGenerator.h"
#include "IFEM.h"
#include "SIMbase.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "tinyxml.h"


ModelGenerator::ModelGenerator (const TiXmlElement* elem) :
  sets(false), geo(elem)
{
  utl::getAttribute(geo, "sets", sets);
}


std::string DefaultGeometry1D::createG2 (int nsd) const
{
  std::string g2("100 1 0 0\n");
  g2.append(1,'0'+nsd);

  bool rational=false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout << "\t Rational basis\n";
  g2.append(rational?" 1":" 0");
  g2.append("\n2 2\n0 0 1 1\n");

  unsigned char d;
  std::string XYZ;
  if (utl::getAttribute(geo,"X0",XYZ))
  {
    IFEM::cout <<"  X0 = "<< XYZ << std::endl;
    g2.append(XYZ);
  }
  else
  {
    g2.append("0.0");
    for (d = 1; d < nsd; d++)
      g2.append(" 0.0");
  }
  if (rational)
    g2.append(" 1.0");
  g2.append("\n");
  if (utl::getAttribute(geo,"X1",XYZ))
  {
    IFEM::cout <<"\tX1 = "<< XYZ << std::endl;
    g2.append(XYZ);
  }
  else
  {
    XYZ = "1.0";
    if (utl::getAttribute(geo,"L",XYZ))
      IFEM::cout <<"  Length scale: "<< XYZ << std::endl;
    g2.append(XYZ);
    for (d = 1; d < nsd; d++)
      g2.append(" 0.0");
  }
  if (rational)
    g2.append(" 1.0");
  g2.append("\n");

  return g2;
}


SIMdependency::PatchVec DefaultGeometry1D::createGeometry (const SIMbase& sim) const
{
  std::istringstream unitLine(this->createG2(sim.getNoSpaceDim()));
  SIMdependency::PatchVec result;
  sim.readPatches(unitLine,result,"\t");
  return result;
}


TopologySet DefaultGeometry1D::createTopologySets (const SIMbase&) const
{
  if (!sets)
    return TopologySet();

  TopologySet result;
  result["Vertex1"].insert(TopItem(1,1,0));
  result["Vertex2"].insert(TopItem(1,2,0));
  result["Boundary"].insert(TopItem(1,1,0));
  result["Boundary"].insert(TopItem(1,2,0));
  result["Corners"].insert(TopItem(1,1,0));
  result["Corners"].insert(TopItem(1,2,0));

  return result;
}


std::string DefaultGeometry2D::createG2 (int nsd) const
{
  std::string g2("200 1 0 0\n");
  g2.append(nsd > 2 ? "3" : "2");
  bool rational=false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout << "\t Rational basis\n";
  g2.append(rational?" 1":" 0");
  g2.append("\n2 2\n0 0 1 1\n2 2\n0 0 1 1");

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"  Corner: "<< X0 << std::endl;
  }

  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"  Scale: "<< scale << std::endl;

  double Lx = 1.0, Ly = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"  Length in X: "<< Lx << std::endl;
  Lx *= scale;
  if (utl::getAttribute(geo,"Ly",Ly))
    IFEM::cout <<"  Length in Y: "<< Ly << std::endl;
  Ly *= scale;

  std::stringstream str;
  str <<"\n"<< X0.x <<" "<< X0.y;
  if (nsd > 2) str <<" 0.0";
  if (rational) str << " 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x+Lx <<" "<< X0.y;
  if (nsd > 2) str <<" 0.0";
  if (rational) str << " 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x <<" "<< X0.y+Ly;
  if (nsd > 2) str <<" 0.0";
  if (rational) str << " 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x+Lx <<" "<< X0.y+Ly;
  if (nsd > 2) str <<" 0.0";
  if (rational) str << " 1.0";
  g2.append(str.str());
  g2.append("\n");

  return g2;
}


SIMdependency::PatchVec DefaultGeometry2D::createGeometry (const SIMbase& sim) const
{
  std::istringstream unitSquare(this->createG2(sim.getNoSpaceDim()));
  SIMdependency::PatchVec result;
  sim.readPatches(unitSquare,result,"\t");
  return result;
}


TopologySet DefaultGeometry2D::createTopologySets (const SIMbase&) const
{
  if (!sets)
    return TopologySet();

  TopologySet result;
  std::string vert = "Vertex1";
  std::string edge = "Edge1";
  for (size_t i = 1; i <= 4; ++i, ++vert.back(), ++edge.back()) {
    result[vert].insert(TopItem(1,i,0));
    result[edge].insert(TopItem(1,i,1));
    result["Corners"].insert(TopItem(1,i,0));
    result["Boundary"].insert(TopItem(1,i,1));
  }

  return result;
}


std::string DefaultGeometry3D::createG2 (int) const
{
  std::string g2("700 1 0 0\n3 ");

  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout <<"  Rational basis"<< std::endl;

  g2.append(rational ? "1\n" : "0\n");
  g2.append("2 2\n0 0 1 1\n"
            "2 2\n0 0 1 1\n"
            "2 2\n0 0 1 1\n");

  std::array<double,24> nodes =
    {{ 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0,
       0.0, 0.0, 1.0,
       1.0, 0.0, 1.0,
       0.0, 1.0, 1.0,
       1.0, 1.0, 1.0 }};

  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"\tscale = "<< scale << std::endl;

  double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"\tLength in X: "<< Lx << std::endl;
  Lx *= scale;
  if (utl::getAttribute(geo,"Ly",Ly))
    IFEM::cout <<"\tLength in Y: "<< Ly << std::endl;
  Ly *= scale;
  if (utl::getAttribute(geo,"Lz",Lz))
    IFEM::cout <<"\tLength in Z: "<< Lz << std::endl;
  Lz *= scale;

  if (Lx != 1.0)
    nodes[3] = nodes[9] = nodes[15] = nodes[21] = Lx;
  if (Ly != 1.0)
    nodes[7] = nodes[10] = nodes[19] = nodes[22] = Ly;
  if (Lz != 1.0)
    nodes[14] = nodes[17] = nodes[20] = nodes[23] = Lz;

  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner);
    Vec3 X0;
    str >> X0;
    IFEM::cout <<"\tCorner: "<< X0 << std::endl;
    for (size_t i = 0; i < nodes.size(); i += 3)
    {
      nodes[i]   += X0.x;
      nodes[i+1] += X0.y;
      nodes[i+2] += X0.z;
    }
  }

  for (size_t i = 0; i < nodes.size(); i += 3)
  {
    std::stringstream str;
    for (size_t j = 0; j < 3; j++)
      str << nodes[i+j] <<" ";
    g2.append(str.str());
    g2.append(rational ? "1.0\n" : "\n");
  }

  return g2;
}


SIMdependency::PatchVec DefaultGeometry3D::createGeometry (const SIMbase& sim) const
{
  std::istringstream unitCube(this->createG2());
  SIMdependency::PatchVec result;
  sim.readPatches(unitCube,result,"\t");
  return result;
}


TopologySet DefaultGeometry3D::createTopologySets (const SIMbase&) const
{
  if (!sets)
    return TopologySet();

  TopologySet result;

  std::string face = "Face1";
  for (size_t i = 1; i <= 6; ++i, ++face.back()) {
    result[face].insert(TopItem(1,i,2));
    result["Boundary"].insert(TopItem(1,i,2));
  }

  std::string edge = "Edge1";
  for (size_t i = 1; i <= 12; ++i, ++edge.back()) {
    result[edge].insert(TopItem(1,i,1));
    result["Frame"].insert(TopItem(1,i,1));
    if (i == 9)
      edge = "Edge1/"; // '/' + 1 == '0'
  }

  std::string vert = "Vertex1";
  for (size_t i = 1; i <= 8; ++i, ++vert.back()) {
    result[vert].insert(TopItem(1,i,0));
    result["Corners"].insert(TopItem(1,i,0));
  }

  return result;
}
