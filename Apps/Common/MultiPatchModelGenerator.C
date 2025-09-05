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
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "IFEM.h"
#include "SIMinput.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "tinyxml2.h"
#include <array>


std::string MultiPatchModelGenerator1D::createG2 (int nsd, bool rational) const
{
  if (rational)
    IFEM::cout <<"\tRational basis.";

  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"\n\tScale: "<< scale;

  double Lx = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx) || utl::getAttribute(geo,"L",Lx))
    IFEM::cout <<"\n\tLength in X: "<< Lx << std::endl;
  Lx *= scale;

  int nx_mp = 1;
  if (nx > 1) {
    IFEM::cout <<"\n\tSplit in X = "<< nx;

    Lx /= nx;

    nx_mp = nx;
  }

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"\n\tCorner: "<< X0;
  }

  std::string g2;
  for (int x = 0; x < nx_mp; ++x) {
    g2.append("100 1 0 0\n");
    g2.append(std::to_string(nsd));
    g2.append(rational?" 1":" 0");
    g2.append("\n2 2\n0 0 1 1");

    std::stringstream str;
    str <<"\n"<< X0.x+x*Lx;
    if (nsd > 1) str <<" 0.0";
    if (nsd > 2) str <<" 0.0";
    if (rational) str << " 1.0";
    g2.append(str.str());
    str.str("");
    str <<"\n"<< X0.x+(x+1)*Lx;
    if (nsd > 1) str <<" 0.0";
    if (nsd > 2) str <<" 0.0";
    if (rational) str << " 1.0";
    g2.append(str.str());
    g2.append("\n");
  }

  return g2;
}


bool MultiPatchModelGenerator1D::createGeometry (SIMinput& sim) const
{
  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  std::istringstream line(this->createG2(sim.getNoSpaceDim(),rational));

  return sim.readPatches(line,"\t");
}


MultiPatchModelGenerator1D::MultiPatchModelGenerator1D (const tinyxml2::XMLElement* geo) :
  ModelGenerator(geo)
{
  nx = 1;
  periodic_x = 0;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"periodic_x", periodic_x);
}


bool MultiPatchModelGenerator1D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  int p1=2;

  for (size_t i = 0; i < nx-1; ++i)
    if (!sim.addConnection(i+1, i+2, 2, 1, true, 0, p1-1))
      return false;

  if (periodic_x)
    if (nx <= 1) {
      IFEM::cout <<"\tPeriodic I-direction P0\n";
      sim.getPatch(1,true)->closeBoundaries();
    }
    else if (!sim.addConnection(1, nx, 1, 2, false, 0, p1-1))
      return false;

  return true;
}


bool MultiPatchModelGenerator1D::createTopologySets (SIMinput& sim) const
{
  if (!this->topologySets())
    return false;

  TopEntity& v1 = sim.topology("Vertex1");
  TopEntity& v2 = sim.topology("Vertex2");
  TopEntity& v3 = sim.topology("Boundary");

  auto&& insertion = [&sim, &v3](TopEntity& v, TopItem top)
  {
    if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
      v.insert(top);
      v3.insert(top);
    }
  };

  insertion(v1, TopItem(1,1,0));
  insertion(v2, TopItem(nx,2,0));

  return true;
}


MultiPatchModelGenerator2D::MultiPatchModelGenerator2D (const tinyxml2::XMLElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = 1;
  periodic_x = periodic_y = 0;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
}


std::string MultiPatchModelGenerator2D::createG2 (int nsd, bool rational) const
{
  if (rational)
    IFEM::cout <<"\tRational basis.";

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

  int nx_mp = 1, ny_mp = 1;
  if (nx*ny > 1) {
    IFEM::cout <<"\n\tSplit in X = "<< nx;
    IFEM::cout <<"\n\tSplit in Y = "<< ny;

    Lx /= nx;
    Ly /= ny;

    nx_mp = nx;
    ny_mp = ny;
  }

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"\n\tCorner = "<< X0;
  }

  std::string g2;
  for (int y = 0; y < ny_mp; ++y) {
    for (int x = 0; x < nx_mp; ++x) {
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

  IFEM::cout << std::endl;
  return g2;
}


bool MultiPatchModelGenerator2D::createGeometry (SIMinput& sim) const
{
  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  std::istringstream rect(this->createG2(sim.getNoSpaceDim(),rational));
  return sim.readPatches(rect,"\t");
}


bool MultiPatchModelGenerator2D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  auto&& IJ = [this](int i, int j) { return 1 + j*nx + i; };

  int p1=2, p2=2;

  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx-1; ++i)
      if (!sim.addConnection(IJ(i,j), IJ(i+1,j), 2, 1, true, 1, p1-1))
        return false;

  for (size_t j = 0; j < ny-1; ++j)
    for (size_t i = 0; i < nx; ++i)
      if (!sim.addConnection(IJ(i,j), IJ(i,j+1), 4, 3, true, 1, p2-1))
        return false;

  if (periodic_x)
    for (size_t j = 0; j < ny; ++j)
      if (nx <= 1) {
        IFEM::cout <<"\tPeriodic I-direction P"<< IJ(0,j) << std::endl;
        sim.getPatch(IJ(0,j),true)->closeBoundaries(1);
      }
      else if (!sim.addConnection(IJ(0,j), IJ(nx-1,j), 1, 2))
        return false;

  if (periodic_y)
    for (size_t i = 0; i < nx; ++i)
      if (ny <= 1) {
        IFEM::cout <<"\tPeriodic J-direction P"<< IJ(i,0)<< std::endl;
        sim.getPatch(IJ(i,0),true)->closeBoundaries(2);
      }
      else if (!sim.addConnection(IJ(i,0), IJ(i,ny-1), 3, 4))
        return false;

  return true;
}


bool MultiPatchModelGenerator2D::createTopologySets (SIMinput& sim) const
{
  if (!this->topologySets())
    return false;

  TopEntity& e1  = sim.topology("Edge1");
  TopEntity& e2  = sim.topology("Edge2");
  TopEntity& e3  = sim.topology("Edge3");
  TopEntity& e4  = sim.topology("Edge4");
  TopEntity& e5  = sim.topology("Boundary");
  TopEntity& e1f = sim.topology("Edge1Patches");
  TopEntity& e2f = sim.topology("Edge2Patches");
  TopEntity& e3f = sim.topology("Edge3Patches");
  TopEntity& e4f = sim.topology("Edge4Patches");
  TopEntity& xb  = sim.topology("BoundaryX");
  TopEntity& yb  = sim.topology("BoundaryY");

  auto&& insertion = [&sim, &e5](TopEntity& e, TopEntity& ef, TopItem top)
  {
    if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
      e.insert(top);
      ef.insert(TopItem(top.patch,0,2));
      e5.insert(top);
    }
  };

  for (size_t i = 0; i < ny; ++i) {
    insertion(e1, e1f, TopItem(i*nx+1,1,1));
    insertion(e2, e2f, TopItem((i+1)*nx,2,1));
  }
  for (const TopItem& e : e1)
    xb.insert(e);
  for (const TopItem& e : e2)
    xb.insert(e);

  for (size_t i = 0; i < nx; ++i) {
    insertion(e3, e3f, TopItem(i+1,3,1));
    insertion(e4, e4f, TopItem(nx*(ny-1)+1+i,4,1));
  }
  for (const TopItem& e : e3)
    yb.insert(e);
  for (const TopItem& e : e4)
    yb.insert(e);

  TopEntity& c = sim.topology("Corners");
  auto&& insertionv = [&sim, &c](TopEntity& e, TopItem top)
  {
    if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
      e.insert(top);
      c.insert(top);
    }
  };

  insertionv(sim.topology("Vertex1"), TopItem(1,1,0));
  insertionv(sim.topology("Vertex2"), TopItem(nx,2,0));
  insertionv(sim.topology("Vertex3"), TopItem(nx*(ny-1)+1,3,0));
  insertionv(sim.topology("Vertex4"), TopItem(nx*ny,4,0));

  std::set<int> used;
  for (size_t i = 1; i <= 4; ++i) {
    std::stringstream str;
    str << "Edge" << i << "Patches";
    for (const TopItem& top : sim.topology(str.str()))
      used.insert(top.patch);
  }

  TopEntity& innerp = sim.topology("InnerPatches");
  for (size_t i = 1; i <= nx*ny; ++i) {
    size_t j = sim.getLocalPatchIndex(i);
    if (j > 0 && used.find(j) == used.end())
      innerp.insert(TopItem(j, 0, 2));
  }

  return true;
}


MultiPatchModelGenerator3D::MultiPatchModelGenerator3D (const tinyxml2::XMLElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = nz = 1;
  periodic_x = periodic_y = periodic_z = 0;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"nz",nz);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
  utl::getAttribute(geo,"periodic_z", periodic_z);
}


std::string MultiPatchModelGenerator3D::createG2 (int, bool rational) const
{
  if (rational)
    IFEM::cout <<"\tRational basis.";

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

  int nx_mp = 1, ny_mp = 1, nz_mp = 1;
  if (nx*ny*nz > 1) {
    IFEM::cout <<"\n\tSplit in X = "<< nx;
    IFEM::cout <<"\n\tSplit in Y = "<< ny;
    IFEM::cout <<"\n\tSplit in Z = "<< nz;

    Lx /= nx;
    Ly /= ny;
    Lz /= nz;

    nx_mp = nx;
    ny_mp = ny;
    nz_mp = nz;
  }

  std::string corner;
  Vec3 X0;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner);
    str >> X0;
    IFEM::cout <<"\n\tCorner = "<< X0;
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
  for (int z = 0; z < nz_mp; ++z) {
    for (int y = 0; y < ny_mp; ++y) {
      for (int x = 0; x < nx_mp; ++x) {
        g2.append("700 1 0 0\n3 ");
        g2.append(rational ? "1\n" : "0\n");
        g2.append("2 2\n0 0 1 1\n"
                  "2 2\n0 0 1 1\n"
                  "2 2\n0 0 1 1\n");

        for (size_t i = 0; i < nodes.size(); i += 3)
        {
          std::stringstream str;
          std::array<int,3> N = {{x,y,z}};
          std::array<double,3> L = {{Lx,Ly,Lz}};
          for (size_t j = 0; j < 3; j++)
            str << (j==0?"":" ") << X0[j]+N[j]*L[j]+nodes[i+j]*L[j];
          g2.append(str.str());
          g2.append(rational ? " 1.0\n" : "\n");
        }
      }
    }
  }

  IFEM::cout << std::endl;
  return g2;
}


bool MultiPatchModelGenerator3D::createGeometry (SIMinput& sim) const
{
  bool rational = false;
  utl::getAttribute(geo,"rational",rational);
  std::istringstream cube(this->createG2(sim.getNoSpaceDim(),rational));
  return sim.readPatches(cube,"\t");
}


bool MultiPatchModelGenerator3D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  auto&& IJK = [this](int i, int j, int k) { return 1 + (k*ny+j)*nx + i; };

  int p1=2, p2=2, p3=2;

  for (size_t k = 0; k < nz; ++k)
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx-1; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i+1,j,k), 2, 1, true, 2, p1-1))
          return false;

  for (size_t k = 0; k < nz; ++k)
    for (size_t j = 0; j < ny-1; ++j)
      for (size_t i = 0; i < nx; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i,j+1,k), 4, 3, true, 2, p2-1))
          return false;

  for (size_t k = 0; k < nz-1; ++k)
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
        if (!sim.addConnection(IJK(i,j,k), IJK(i,j,k+1), 6, 5, true, 2, p3-1))
          return false;

  if (periodic_x)
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        if (nx <= 1) {
          IFEM::cout <<"\tPeriodic I-direction P"<< IJK(0,j,k) << std::endl;
          sim.getPatch(IJK(0,j,k),true)->closeBoundaries(1);
        }
        else if (!sim.addConnection(IJK(0,j,k), IJK(nx-1,j,k), 1, 2, false, 2))
          return false;

  if (periodic_y)
    for (size_t k = 0; k < nz; ++k)
      for (size_t i = 0; i < nx; ++i)
        if (ny <= 1) {
          IFEM::cout <<"\tPeriodic J-direction P"<< IJK(i,0,k) << std::endl;
          sim.getPatch(IJK(i,0,k),true)->closeBoundaries(2);
        }
        else if (!sim.addConnection(IJK(i,0,k), IJK(i,ny-1,k), 3, 4, false, 2))
          return false;

  if (periodic_z)
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
        if (nz <= 1) {
          IFEM::cout <<"\tPeriodic K-direction P"<< IJK(i,j,0) << std::endl;
          sim.getPatch(IJK(i,j,0),true)->closeBoundaries(3);
        }
        else if (!sim.addConnection(IJK(i,j,0), IJK(i,j,nz-1), 5, 6, false, 2))
          return false;

  return true;
}


bool MultiPatchModelGenerator3D::createTopologySets (SIMinput& sim) const
{
  if (!this->topologySets())
    return false;

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

  // insertion lambda
  auto&& insertion = [&sim](TopItem top,
                            const std::string& glob,
                            const std::string& type)
    {
      std::stringstream str;
      str << type << top.item;
      TopEntity& topI = sim.topology(str.str());
      TopEntity* topIf=nullptr;
      if (type == "Face")
        topIf = &sim.topology(str.str()+"Patches");
      TopEntity& globI = sim.topology(glob);
      if ((top.patch = sim.getLocalPatchIndex(top.patch)) > 0) {
        topI.insert(top);
        globI.insert(top);
        if (topIf)
          topIf->insert(TopItem(top.patch, 0, 3));
      }
    };

  size_t r = 1;
  for (size_t i = 0; i < 2; ++i, ++r)
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        insertion(TopItem(IJK2I(i,j,k),r,2), "Boundary", "Face");

  for (size_t j = 0; j < 2; ++j, ++r)
    for (size_t k = 0; k < nz; ++k)
      for (size_t i = 0; i < nx; ++i)
        insertion(TopItem(IJK2J(i,j,k),r,2), "Boundary", "Face");

  for (size_t k = 0; k < 2; ++k, ++r)
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
        insertion(TopItem(IJK2K(i,j,k),r,2), "Boundary", "Face");

  TopEntity& bx  = sim.topology("BoundaryX");
  TopEntity& by  = sim.topology("BoundaryY");
  TopEntity& bz  = sim.topology("BoundaryZ");

  auto&& side_insertion = [&sim](TopEntity& set, int start)
  {
    for (int i = start; i < start+2; ++i) {
      for (const auto& entry : sim.topology("Face" + std::to_string(i)))
        set.insert(entry);
    }
  };

  side_insertion(bx, 1);
  side_insertion(by, 3);
  side_insertion(bz, 5);

  r = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i, ++r)
        insertion(TopItem(IJK2(i,j,k),r,0), "Corners", "Vertex");

  r = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t i = 0; i < 2; ++i, ++r)
      for (size_t j = 0; j < ny; ++j)
        insertion(TopItem(IJKJ(i,j,k),r,1), "Frame", "Edge");

  for (size_t j = 0; j < 2; ++j)
    for (size_t i = 0; i < 2; ++i, ++r)
      for (size_t k = 0; k < nz; ++k)
        insertion(TopItem(IJKK(i,j,k),r,1), "Frame", "Edge");

  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j, ++r)
      for (size_t i = 0; i < nx; ++i)
        insertion(TopItem(IJKI(i,j,k),r,1), "Frame", "Edge");

  std::set<int> used;
  for (size_t i = 1; i <= 6; ++i) {
    std::stringstream str;
    str << "Face" << i << "Patches";
    for (const TopItem& top : sim.topology(str.str()))
      used.insert(top.patch);
  }

  TopEntity& innerp = sim.topology("InnerPatches");
  for (size_t i = 1; i <= nx*ny*nz; ++i) {
    size_t j = sim.getLocalPatchIndex(i);
    if (j > 0 && used.find(j) == used.end())
      innerp.insert(TopItem(j, 0, 3));
  }

  return true;
}
