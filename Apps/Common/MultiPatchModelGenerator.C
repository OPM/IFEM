// $Id$
//==============================================================================
//!
//! \file MultiPatchModelGenerator.C
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Multi-patch model generators for NURBS-based FEM simulators.
//!
//==============================================================================

#include "MultiPatchModelGenerator.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "IFEM.h"
#include "SIMbase.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "tinyxml.h"


MultiPatchModelGenerator2D::MultiPatchModelGenerator2D (const TiXmlElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = 1;
  periodic_x = periodic_y = 0;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
}


std::string MultiPatchModelGenerator2D::createG2 (int nsd) const
{
  bool rational=false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout << "\t Rational basis\n";
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

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"  Corner: "<< X0 << std::endl;
  }

  int nx = 1;
  int ny = 1;
  if (utl::getAttribute(geo,"nx",nx))
    IFEM::cout << "  Split in X: " << nx  << std::endl;
  if (utl::getAttribute(geo,"ny",ny))
    IFEM::cout << "  Split in Y: " << ny << std::endl;

  if (nx > 1)
    Lx /= nx;
  if (ny > 1)
    Ly /= ny;

  std::string g2;
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x) {
      g2.append("200 1 0 0\n");
      g2.append(nsd > 2 ? "3" : "2");
      g2.append(rational?" 1":" 0");
      g2.append("\n2 2\n0 0 1 1\n2 2\n0 0 1 1");

      std::stringstream str;
      str <<"\n"<< X0.x+x*Lx <<" "<< X0.y+y*Ly;
      if (nsd > 2) str <<" 0.0";
      if (rational) str << " 1.0";
      g2.append(str.str());
      str.str("");
      str <<"\n"<< X0.x+(x+1)*Lx <<" "<< X0.y+y*Ly;
      if (nsd > 2) str <<" 0.0";
      if (rational) str << " 1.0";
      g2.append(str.str());
      str.str("");
      str <<"\n"<< X0.x+x*Lx <<" "<< X0.y+(y+1)*Ly;
      if (nsd > 2) str <<" 0.0";
      if (rational) str << " 1.0";
      g2.append(str.str());
      str.str("");
      str <<"\n"<< X0.x+(x+1)*Lx <<" "<< X0.y+(y+1)*Ly;
      if (nsd > 2) str <<" 0.0";
      if (rational) str << " 1.0";
      g2.append(str.str());
      g2.append("\n");
    }
  }

  return g2;
}


SIMdependency::PatchVec
MultiPatchModelGenerator2D::createGeometry (const SIMbase& sim) const
{
  std::istringstream rect(this->createG2(sim.getNoSpaceDim()));
  SIMdependency::PatchVec result;
  sim.readPatches(rect,result,"\t");
  return result;
}


bool MultiPatchModelGenerator2D::createTopology (SIMbase& sim) const
{
  auto&& IJ = [this](int i, int j) { return 1 + j*nx + i; };

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx-1; ++i)
      if (!sim.addConnection(IJ(i,j), IJ(i+1,j), 2, 1, 0))
        return false;

  for (int j = 0; j < ny-1; ++j)
    for (int i = 0; i < nx; ++i)
      if (!sim.addConnection(IJ(i,j), IJ(i,j+1), 4, 3, 0))
        return false;

  if (periodic_x)
    for (int i = 0; i < ny; ++i)
      if (nx > 1) {
        if (!sim.addConnection(IJ(0, i), IJ(nx-1, i), 1, 2, 0, false))
          return false;
      } else {
         IFEM::cout <<"\tPeriodic I-direction P"<< IJ(0,i) << std::endl;
         ASMs2D* pch = dynamic_cast<ASMs2D*>(sim.getPatch(IJ(0,i), true));
         if (pch)
           pch->closeEdges(1);
      }

  if (periodic_y)
    for (int i = 0; i < nx; ++i)
      if (ny > 1)
        if (!sim.addConnection(IJ(i,0), IJ(i,ny-1), 3, 4, 0, false))
          return false;
      else {
         IFEM::cout <<"\tPeriodic J-direction P"<< IJ(i,0)<< std::endl;
         ASMs2D* pch = dynamic_cast<ASMs2D*>(sim.getPatch(IJ(i,0), true));
         if (pch)
           pch->closeEdges(2);
      }

  return true;
}


TopologySet
MultiPatchModelGenerator2D::createTopologySets (const SIMbase& sim) const
{
  TopologySet result;
  TopEntity& e1 = result["Edge1"];
  TopEntity& e2 = result["Edge2"];
  TopEntity& e3 = result["Edge3"];
  TopEntity& e4 = result["Edge4"];
  TopEntity& e5 = result["Boundary"];

  auto&& insertion = [&sim, &e5](TopEntity& e, TopItem top)
  {
    if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
      e.insert(top);
      e5.insert(top);
    }
  };

  for (int i = 0; i < ny; ++i) {
    insertion(e1, TopItem(i*nx+1,1,1));
    insertion(e2, TopItem((i+1)*nx,2,1));
  }
  for (int i = 0; i < nx; ++i) {
    insertion(e3, TopItem(i+1,3,1));
    insertion(e4, TopItem(nx*(ny-1)+1+i,4,1));
  }

  TopEntity& c = result["Corners"];
  auto&& insertionv = [&sim, &c](TopEntity& e, TopItem top)
  {
    if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
      e.insert(top);
      c.insert(top);
    }
  };

  insertionv(result["Vertex1"], TopItem(1,1,0));
  insertionv(result["Vertex2"], TopItem(nx,2,0));
  insertionv(result["Vertex3"], TopItem(nx*(ny-1)+1,3,0));
  insertionv(result["Vertex4"], TopItem(nx*ny,4,0));

  return result;
}


MultiPatchModelGenerator3D::MultiPatchModelGenerator3D (const TiXmlElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = nz = 1;
  periodic_x = periodic_y = periodic_z = 0;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"nz",nz);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
  utl::getAttribute(geo,"periodic_z", periodic_z);
}


std::string MultiPatchModelGenerator3D::createG2 (int) const
{
  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  if (rational)
    IFEM::cout <<"  Rational basis"<< std::endl;

  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"  Scale: "<< scale << std::endl;

  double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"  Length in X: "<< Lx << std::endl;
  Lx *= scale;
  if (utl::getAttribute(geo,"Ly",Ly))
    IFEM::cout <<"  Length in Y: "<< Ly << std::endl;
  Ly *= scale;
  if (utl::getAttribute(geo,"Lz",Lz))
    IFEM::cout <<"  Length in Z: "<< Lz << std::endl;
  Lz *= scale;

  int nx = 1;
  int ny = 1;
  int nz = 1;
  if (utl::getAttribute(geo,"nx",nx))
    IFEM::cout << "  Split in X: " << nx  << std::endl;
  if (utl::getAttribute(geo,"ny",ny))
    IFEM::cout << "  Split in Y: " << ny << std::endl;
  if (utl::getAttribute(geo,"nz",nz))
    IFEM::cout << "  Split in Z: " << nz << std::endl;

  Lx /= nx;
  Ly /= ny;
  Lz /= nz;

  std::string corner;
  Vec3 X0;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner);
    str >> X0;
    IFEM::cout <<"  Corner: "<< X0 << std::endl;
  }

  std::array<double,24> nodes =
    {{ 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0,
       0.0, 0.0, 1.0,
       1.0, 0.0, 1.0,
       0.0, 1.0, 1.0,
       1.0, 1.0, 1.0 }};

  std::string g2;
  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        g2.append("700 1 0 0\n3 ");
        g2.append(rational ? "1\n" : "0\n");
        g2.append("2 2\n0 0 1 1\n"
                  "2 2\n0 0 1 1\n"
                  "2 2\n0 0 1 1\n");

        for (size_t i = 0; i < nodes.size(); i += 3)
        {
          std::stringstream str;
          std::array<int,3> N {x,y,z};
          std::array<double,3> L {Lx,Ly,Lz};
          for (size_t j = 0; j < 3; j++)
            str << (j==0?"":" ") << X0[j]+N[j]*L[j]+nodes[i+j]*L[j];
          g2.append(str.str());
          g2.append(rational ? " 1.0\n" : "\n");
        }
      }
    }
  }

  return g2;
}


SIMdependency::PatchVec
MultiPatchModelGenerator3D::createGeometry (const SIMbase& sim) const
{
  std::istringstream hex(this->createG2(sim.getNoSpaceDim()));
  SIMdependency::PatchVec result;
  sim.readPatches(hex,result,"\t");
  return result;
}


bool MultiPatchModelGenerator3D::createTopology (SIMbase& sim) const
{
  auto&& IJK = [this](int i, int j, int k) { return 1 + (k*ny+j)*nx + i; };

  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx-1; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i+1,j,k), 2, 1, 0))
          return false;

  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny-1; ++j)
      for (int i = 0; i < nx; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i,j+1,k), 4, 3, 0))
          return false;

  for (int k = 0; k < nz-1; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i,j,k+1), 6, 5, 0))
          return false;

  if (periodic_x)
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        if (nx > 1) {
          if (!sim.addConnection(IJK(0,j,k), IJK(nx-1,j,k), 1, 2, 0, false))
            return false;
        } else {
          IFEM::cout <<"\tPeriodic I-direction P"<< IJK(0,j,k) << std::endl;
          ASMs3D* pch = dynamic_cast<ASMs3D*>(sim.getPatch(IJK(0,j,k), true));
          if (pch)
            pch->closeFaces(1);
        }

  if (periodic_y)
    for (int k = 0; k < nz; ++k)
      for (int i = 0; i < nx; ++i)
        if (ny > 1) {
          if (!sim.addConnection(IJK(i,0,k), IJK(i,ny-1,k), 3, 4, 0, false))
            return false;
         } else {
          IFEM::cout <<"\tPeriodic J-direction P"<< IJK(i,0,k) << std::endl;
          ASMs3D* pch = dynamic_cast<ASMs3D*>(sim.getPatch(IJK(i,0,k), true));
          if (pch)
            pch->closeFaces(2);
         }

  if (periodic_z)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        if (nz > 1) {
          if (!sim.addConnection(IJK(i,j,0), IJK(i,j,nz-1), 5, 6, 0, false))
            return false;
        } else {
          IFEM::cout <<"\tPeriodic K-direction P"<< IJK(i,j,0) << std::endl;
          ASMs3D* pch = dynamic_cast<ASMs3D*>(sim.getPatch(IJK(i,j,0), true));
          if (pch)
            pch->closeFaces(3);
        }

  return true;

}

TopologySet
MultiPatchModelGenerator3D::createTopologySets (const SIMbase& sim) const
{
  // 0-based -> 1-based IJK
  auto&& IJK = [this](int i, int j, int k) { return 1 + (k*ny+j)*nx + i; };

  // start/end IJK
  auto&& IJK2 = [this,IJK](int i, int j, int k) { return IJK(i*(nx-1), j*(ny-1), k*(nz-1)); };

  // start/end JK
  auto&& IJKI = [this,IJK](int i, int j, int k) { return IJK(i, j*(ny-1), k*(nz-1)); };
  // start/end IK
  auto&& IJKJ = [this,IJK](int i, int j, int k) { return IJK(i*(nx-1), j, k*(nz-1)); };
  // start/end IJ
  auto&& IJKK = [this,IJK](int i, int j, int k) { return IJK(i*(nx-1), j*(ny-1), k); };

  // start/end I
  auto&& IJK2I = [this,IJK](int i, int j, int k) { return IJK(i*(nx-1), j, k); };
  // start/end J
  auto&& IJK2J = [this,IJK](int i, int j, int k) { return IJK(i, j*(ny-1), k); };
  // start/end K
  auto&& IJK2K = [this,IJK](int i, int j, int k) { return IJK(i, j, k*(nz-1)); };

  TopologySet result;

  // insertion lambda
  auto&& insertion = [&sim,&result](TopItem top,
                                    const std::string& glob,
                                    const std::string& type)
                     {
                       std::stringstream str;
                       str << type << top.item;
                       TopEntity& topI = result[str.str()];
                       TopEntity& globI = result[glob];
                       if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
                         topI.insert(top);
                         globI.insert(top);
                       }
                     };

  size_t r = 1;
  for (int i = 0; i < 2; ++i, ++r)
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        insertion(TopItem(IJK2I(i,j,k),r,2), "Boundary", "Face");

  for (int j = 0; j < 2; ++j, ++r)
    for (int k = 0; k < nz; ++k)
      for (int i = 0; i < nx; ++i)
        insertion(TopItem(IJK2J(i,j,k),r,2), "Boundary", "Face");

  for (int k = 0; k < 2; ++k, ++r)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        insertion(TopItem(IJK2K(i,j,k),r,2), "Boundary", "Face");

  r = 1;
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i, ++r)
        insertion(TopItem(IJK2(i,j,k),r,0), "Corners", "Vertex");

  r = 1;
  for (int k = 0; k < 2; ++k)
    for (int i = 0; i < 2; ++i, ++r)
      for (int j = 0; j < ny; ++j)
        insertion(TopItem(IJKJ(i,j,k),r,1), "Frame", "Edge");

  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 2; ++i, ++r)
      for (int k = 0; k < nz; ++k)
        insertion(TopItem(IJKK(i,j,k),r,1), "Frame", "Edge");

  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j, ++r)
      for (int i = 0; i < nx; ++i)
        insertion(TopItem(IJKI(i,j,k),r,1), "Frame", "Edge");

  return result;
}
