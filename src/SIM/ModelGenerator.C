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
#include "SIMinput.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <array>


/*!
  \brief Static helper adding knot vectors to a g2 definition.
  \param g2 The string to add knot vectors to
  \param[in] num The number of knot vectors to add
  \param[in] geo XML element describing geometry
*/

static void AddKnots(std::string& g2, size_t num, const TiXmlElement* geo)
{
  std::string min = "umin", max = "umax";
  std::stringstream str;
  for (size_t i = 0; i < num; ++i, ++min[0], ++max[0]) {
    double pmin = 0.0, pmax = 1.0;
    if (utl::getAttribute(geo,min.c_str(),pmin) |
        utl::getAttribute(geo,max.c_str(),pmax))
      IFEM::cout <<"\n\t" << min[0] << " = ["<< pmin << ","<< pmax <<"]";

    str << "\n2 2\n" << pmin << " " << pmin << " " << pmax << " " << pmax;
  }
  g2.append(str.str());
  g2.append("\n");
}


bool ModelGenerator::topologySets () const
{
  bool sets = false;
  return utl::getAttribute(geo,"sets",sets) && sets;
}


SIMdependency::PatchVec ModelGenerator::createGeometry (const SIMinput& m) const
{
  std::istringstream g2(this->createG2(m.getNoSpaceDim()));
  SIMdependency::PatchVec result;
  m.readPatches(g2,result,nullptr);
  return result;
}


std::string DefaultGeometry1D::createG2 (int nsd) const
{
  IFEM::cout <<"  Generating linear geometry on unit parameter domain [0,1]";

  std::string g2("100 1 0 0\n");
  g2.append(1,'0'+nsd);

  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout <<"\n\tRational basis.";
  g2.append(rational ? " 1" : " 0");

  AddKnots(g2, 1, geo);

  unsigned char d;
  std::string X0("0"), X1("1");
  if (utl::getAttribute(geo,"X0",X0))
  {
    IFEM::cout <<"\n\tX0 = "<< X0;
    g2.append(X0);
  }
  else
  {
    g2.append("0");
    for (d = 1; d < nsd; d++)
      g2.append(" 0");
  }
  g2.append(rational ? " 1.0\n" : "\n");
  if (utl::getAttribute(geo,"X1",X1))
  {
    IFEM::cout <<"\n\tX1 = "<< X1;
    g2.append(X1);
  }
  else
  {
    double Xend, L = 1.0;
    std::stringstream is(X0); is >> Xend;
    if (utl::getAttribute(geo,"L",L))
      IFEM::cout <<"\n\tLength = "<< L;
    Xend += L;
    std::stringstream os; os << Xend;
    g2.append(os.str());
    for (d = 1; d < nsd; d++)
      g2.append(" 0");
  }
  g2.append(rational ? " 1.0\n" : "\n");

  IFEM::cout << std::endl;
  return g2;
}


TopologySet DefaultGeometry1D::createTopologySets (const SIMinput&) const
{
  TopologySet result;
  if (this->topologySets())
  {
    result["Vertex1"].insert(TopItem(1,1,0));
    result["Vertex2"].insert(TopItem(1,2,0));
    result["Boundary"].insert(TopItem(1,1,0));
    result["Boundary"].insert(TopItem(1,2,0));
    result["Corners"].insert(TopItem(1,1,0));
    result["Corners"].insert(TopItem(1,2,0));
  }

  return result;
}


std::string DefaultGeometry2D::createG2 (int nsd) const
{
  IFEM::cout <<"  Generating linear geometry on unit parameter domain [0,1]^2";

  std::string g2("200 1 0 0\n");
  g2.append(nsd > 2 ? "3" : "2");

  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout <<"\n\tRational basis.";
  g2.append(rational ? " 1" : " 0");

  AddKnots(g2, 2, geo);

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner))
  {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"\n\tCorner = "<< X0;
  }

  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"\n\tScale = "<< scale;

  double Lx = 1.0, Ly = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"\n\tLength in X = "<< Lx;
  Lx *= scale;
  if (utl::getAttribute(geo,"Ly",Ly))
    IFEM::cout <<"\n\tLength in Y = "<< Ly;
  Ly *= scale;

  std::stringstream str;
  str << X0.x <<" "<< X0.y;
  if (nsd > 2) str <<" 0.0";
  if (rational) str <<" 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x+Lx <<" "<< X0.y;
  if (nsd > 2) str <<" 0.0";
  if (rational) str <<" 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x <<" "<< X0.y+Ly;
  if (nsd > 2) str <<" 0.0";
  if (rational) str <<" 1.0";
  g2.append(str.str());
  str.str("");
  str <<"\n"<< X0.x+Lx <<" "<< X0.y+Ly;
  if (nsd > 2) str <<" 0.0";
  if (rational) str <<" 1.0";
  g2.append(str.str());
  g2.append("\n");

  IFEM::cout << std::endl;
  return g2;
}


TopologySet DefaultGeometry2D::createTopologySets (const SIMinput&) const
{
  TopologySet result;
  if (this->topologySets())
  {
    std::string vert = "Vertex1";
    std::string edge = "Edge1";
    for (size_t i = 1; i <= 4; ++i, ++vert.back(), ++edge.back())
    {
      result[vert].insert(TopItem(1,i,0));
      result[edge].insert(TopItem(1,i,1));
      result["Corners"].insert(TopItem(1,i,0));
      result["Boundary"].insert(TopItem(1,i,1));
    }
  }

  return result;
}


std::string DefaultGeometry3D::createG2 (int) const
{
  IFEM::cout <<"  Generating linear geometry on unit parameter domain [0,1]^3";

  std::string g2("700 1 0 0\n3 ");

  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout <<"\n\tRational basis.";
  g2.append(rational ? "1" : "0");

  AddKnots(g2, 3, geo);

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
    IFEM::cout <<"\n\tScale = "<< scale;

  double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"\n\tLength in X = "<< Lx;
  Lx *= scale;
  if (utl::getAttribute(geo,"Ly",Ly))
    IFEM::cout <<"\n\tLength in Y = "<< Ly;
  Ly *= scale;
  if (utl::getAttribute(geo,"Lz",Lz))
    IFEM::cout <<"\n\tLength in Z = "<< Lz;
  Lz *= scale;

  if (Lx != 1.0)
    nodes[3] = nodes[9] = nodes[15] = nodes[21] = Lx;
  if (Ly != 1.0)
    nodes[7] = nodes[10] = nodes[19] = nodes[22] = Ly;
  if (Lz != 1.0)
    nodes[14] = nodes[17] = nodes[20] = nodes[23] = Lz;

  std::string corner;
  if (utl::getAttribute(geo,"X0",corner))
  {
    std::stringstream str(corner);
    Vec3 X0;
    str >> X0;
    IFEM::cout <<"\n\tCorner = "<< X0;
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

  IFEM::cout << std::endl;
  return g2;
}


TopologySet DefaultGeometry3D::createTopologySets (const SIMinput&) const
{
  TopologySet result;
  if (this->topologySets())
  {
    std::string face = "Face1";
    for (size_t i = 1; i <= 6; ++i, ++face.back())
    {
      result[face].insert(TopItem(1,i,2));
      result["Boundary"].insert(TopItem(1,i,2));
    }

    std::string edge = "Edge1";
    for (size_t i = 1; i <= 12; ++i, ++edge.back())
    {
      result[edge].insert(TopItem(1,i,1));
      result["Frame"].insert(TopItem(1,i,1));
      if (i == 9)
        edge = "Edge1/"; // '/' + 1 == '0'
    }

    std::string vert = "Vertex1";
    for (size_t i = 1; i <= 8; ++i, ++vert.back())
    {
      result[vert].insert(TopItem(1,i,0));
      result["Corners"].insert(TopItem(1,i,0));
    }
  }

  return result;
}
