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
#include "tinyxml.h"
#include <array>


std::string MultiPatchModelGenerator1D::createG2 (int nsd, bool rational) const
{
  if (rational)
    IFEM::cout <<"\tRational basis.";
  double scale = 1.0;
  if (utl::getAttribute(geo,"scale",scale))
    IFEM::cout <<"\n\tScale: "<< scale;

  double Lx = 1.0;
  if (utl::getAttribute(geo,"Lx",Lx))
    IFEM::cout <<"\n\tLength in X: "<< Lx << std::endl;
  Lx *= scale;

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"\n\tCorner: "<< X0;
  }

  int nx_mp = 1;
  if (!subdivision) {
    IFEM::cout <<"\n\tSplit in X = "<< nx;

    Lx /= nx;

    nx_mp = nx;
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

  if (!subdivision)
    return sim.readPatches(line,"\t");

  // for now this consists of a single patch. do refine / raiseorder
  // to obtain knot vector, then split in pieces. written this way
  // so code can be used with .g2 files.

  IFEM::cout << "  Subdivision in X: " << nx << std::endl;
  ASMs1D pch;
  pch.read(line);

  // Parse XML input
  const TiXmlElement* sub = geo->FirstChildElement("subdivision")->FirstChildElement();
  for (; sub; sub = sub->NextSiblingElement())
    if (strcasecmp(sub->Value(),"raiseorder") == 0) {
      int nu;
      utl::getAttribute(sub,"u",nu);
      pch.raiseOrder(nu);
    } else if (strcasecmp(sub->Value(),"refine") == 0) {
      int nu;
      utl::getAttribute(sub,"u",nu);
      pch.uniformRefine(nu);
    }

  // Compute parameters
  const Go::SplineCurve* cur = pch.getCurve();
  size_t p = cur->order() - 1;
  size_t nelems = cur->numCoefs() - p;
  size_t nelems_sub = nelems / nx;
  size_t nelems_rem = nelems % nx;
  std::string header = "100 1 0 0\n";

  // Extract subpatches
  std::stringstream str;
  for (size_t i = 0; i < nx; ++i) {
    size_t ni = nelems_sub + (i < nelems_rem ? 1 : 0);
    size_t i0 = ni*i + (i < nelems_rem ? 0 : nelems_rem);
    size_t di = ni + p;

    Go::SplineCurve subcur = getSubPatch(cur, i0, di, p+1);
    IFEM::cout << "  Number of knot spans in patch " << i << ": "
               << subcur.numCoefs()-subcur.order()+1 << std::endl;
    str << header << subcur;
  }

  return sim.readPatches(str,"\t");
}


Go::SplineCurve
MultiPatchModelGenerator1D::getSubPatch(const Go::SplineCurve* cur,
                                        const size_t startu,
                                        const size_t numcoefsu, const int orderu)
{
  size_t d = cur->dimension() + (cur->rational()?1:0); // #coefs per control point
  bool rat = cur->rational();
  std::vector<double>::const_iterator i0 = cur->basis().begin() + startu;
  std::vector<double>::const_iterator j0 = rat ? cur->rcoefs_begin() : cur->coefs_begin();
  return Go::SplineCurve(numcoefsu, orderu, i0, j0 + d*startu, cur->dimension(), rat);
}


MultiPatchModelGenerator1D::MultiPatchModelGenerator1D (const TiXmlElement* geo) :
  ModelGenerator(geo)
{
  nx = 1;
  periodic_x = 0;
  subdivision = false;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"subdivision",subdivision);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  const TiXmlElement* subd = geo->FirstChildElement("subdivision");
  if (subd) {
    subdivision = true;
    utl::getAttribute(subd,"nx",nx);
  }
}


bool MultiPatchModelGenerator1D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  int p1=2;
  if (subdivision)
    sim.getPatch(1)->getOrder(p1,p1,p1);

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


MultiPatchModelGenerator2D::MultiPatchModelGenerator2D (const TiXmlElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = 1;
  periodic_x = periodic_y = 0;
  subdivision = false;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"subdivision",subdivision);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
  const TiXmlElement* subd = geo->FirstChildElement("subdivision");
  if (subd) {
    subdivision = true;
    utl::getAttribute(subd,"nx",nx);
    utl::getAttribute(subd,"ny",ny);
  }
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

  Vec3 X0;
  std::string corner;
  if (utl::getAttribute(geo,"X0",corner)) {
    std::stringstream str(corner); str >> X0;
    IFEM::cout <<"\n\tCorner = "<< X0;
  }

  int nx_mp = 1, ny_mp = 1;
  if (!subdivision) {
    IFEM::cout <<"\n\tSplit in X = "<< nx;
    IFEM::cout <<"\n\tSplit in Y = "<< ny;

    Lx /= nx;
    Ly /= ny;

    nx_mp = nx;
    ny_mp = ny;
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
  if (!subdivision)
    return sim.readPatches(rect,"\t");

  // for now this consists of a single patch. do refine / raiseorder
  // to obtain knot vector, then split in pieces. written this way
  // so code can be used with .g2 files.
  IFEM::cout << "  Subdivision in X: " << nx << std::endl;
  IFEM::cout << "  Subdivision in Y: " << ny << std::endl;
  ASMs2D pch;
  pch.read(rect);

  // Parse XML input
  const TiXmlElement* sub = geo->FirstChildElement("subdivision")->FirstChildElement();
  for (; sub; sub = sub->NextSiblingElement())
    if (strcasecmp(sub->Value(),"raiseorder") == 0) {
      int nu, nv;
      utl::getAttribute(sub,"u",nu);
      utl::getAttribute(sub,"v",nv);
      pch.raiseOrder(nu,nv);
    } else if (strcasecmp(sub->Value(),"refine") == 0) {
      int nu, nv;
      utl::getAttribute(sub,"u",nu);
      utl::getAttribute(sub,"v",nv);
      pch.uniformRefine(0,nu);
      pch.uniformRefine(1,nv);
    }

  // Compute parameters
  const Go::SplineSurface* srf = pch.getBasis(ASM::GEOMETRY_BASIS);
  size_t px = srf->order_u()-1;
  size_t py = srf->order_v()-1;
  size_t nelemsx = srf->numCoefs_u() - px;
  size_t nelemsy = srf->numCoefs_v() - py;
  size_t nelemsx_sub = nelemsx / nx;
  size_t nelemsy_sub = nelemsy / ny;
  size_t nelemsx_rem = nelemsx % nx;
  size_t nelemsy_rem = nelemsy % ny;
  std::string header = "200 1 0 0\n";

  // Extract subpatches
  std::stringstream str;
  for (size_t j = 0; j < ny; ++j) {
    size_t nj = nelemsy_sub + (j < nelemsy_rem ? 1 : 0);
    size_t j0 = nj*j + (j < nelemsy_rem ? 0 : nelemsy_rem);
    size_t dj = nj + py;
    for (size_t i = 0; i < nx; ++i) {
      size_t ni = nelemsx_sub + (i < nelemsx_rem ? 1 : 0);
      size_t i0 = ni*i + (i < nelemsx_rem ? 0 : nelemsx_rem);
      size_t di = ni + px;

      Go::SplineSurface subsrf = this->getSubPatch(srf, {{i0, j0}}, {{di, dj}},
                                                   {{px+1, py+1}});
      IFEM::cout << "  Number of knot spans in patch (" << i << ", " << j << "): "
                 << subsrf.numCoefs_u()-subsrf.order_u()+1 << "x"
                 << subsrf.numCoefs_v()-subsrf.order_v()+1 << std::endl;
      str << header << subsrf;
    }
  }

  return sim.readPatches(str,"\t");
}


Go::SplineSurface
MultiPatchModelGenerator2D::getSubPatch(const Go::SplineSurface* srf,
                                        const std::array<size_t,2>& start,
                                        const std::array<size_t,2>& numcoefs,
                                        const std::array<size_t,2>& order)
{
  size_t d = srf->dimension() + (srf->rational()?1:0); // #coefs per control point
  bool rat = srf->rational();
  size_t I = srf->numCoefs_u();
  std::vector<double>::const_iterator i0 = srf->basis_u().begin() + start[0];
  std::vector<double>::const_iterator j0 = srf->basis_v().begin() + start[1];
  std::vector<double>::const_iterator k0 = rat ? srf->rcoefs_begin() : srf->coefs_begin();

  std::vector<double> subcoefs;
  for (size_t j = start[1]; j != start[1] + numcoefs[1]; j++)
    for (size_t i = start[0]; i != start[0] + numcoefs[0]; i++)
      for (size_t k = 0;  k != d; k++)
        subcoefs.push_back(*(k0 + k + d*i + d*I*j));
  return Go::SplineSurface(numcoefs[0], numcoefs[1], order[0], order[1], i0, j0,
                           subcoefs.begin(), srf->dimension(), rat);
}


bool MultiPatchModelGenerator2D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  auto&& IJ = [this](int i, int j) { return 1 + j*nx + i; };

  int p1=2, p2=2;
  if (subdivision)
    sim.getPatch(1)->getOrder(p1,p2,p2);

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
  for (size_t i = 0; i < nx; ++i) {
    insertion(e3, e3f, TopItem(i+1,3,1));
    insertion(e4, e4f, TopItem(nx*(ny-1)+1+i,4,1));
  }

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


MultiPatchModelGenerator3D::MultiPatchModelGenerator3D (const TiXmlElement* geo) :
  ModelGenerator(geo)
{
  nx = ny = nz = 1;
  periodic_x = periodic_y = periodic_z = 0;
  subdivision = false;
  if (!geo) return;
  utl::getAttribute(geo,"nx",nx);
  utl::getAttribute(geo,"ny",ny);
  utl::getAttribute(geo,"nz",nz);
  utl::getAttribute(geo,"periodic_x", periodic_x);
  utl::getAttribute(geo,"periodic_y", periodic_y);
  utl::getAttribute(geo,"periodic_z", periodic_z);
  utl::getAttribute(geo,"subdivision", subdivision);
  const TiXmlElement* subd = geo->FirstChildElement("subdivision");
  if (subd) {
    subdivision = true;
    utl::getAttribute(subd,"nx",nx);
    utl::getAttribute(subd,"ny",ny);
    utl::getAttribute(subd,"nz",nz);
  }
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
  if (!subdivision) {
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
  if (!subdivision)
    return sim.readPatches(cube,"\t");

  // for now this consists of a single patch. do refine / raiseorder
  // to obtain knot vector, then split in pieces. written this way
  // so code can be used with .g2 files.
  IFEM::cout << "  Subdivision in X: " << nx << std::endl;
  IFEM::cout << "  Subdivision in Y: " << ny << std::endl;
  IFEM::cout << "  Subdivision in Z: " << nz << std::endl;
  ASMs3D pch;
  pch.read(cube);

  // Parse XML input
  const TiXmlElement* sub = geo->FirstChildElement("subdivision")->FirstChildElement();
  for (; sub; sub = sub->NextSiblingElement())
    if (strcasecmp(sub->Value(),"raiseorder") == 0) {
      int nu, nv, nw;
      utl::getAttribute(sub,"u",nu);
      utl::getAttribute(sub,"v",nv);
      utl::getAttribute(sub,"w",nw);
      pch.raiseOrder(nu,nv,nw);
    } else if (strcasecmp(sub->Value(),"refine") == 0) {
      int nu, nv, nw;
      utl::getAttribute(sub,"u",nu);
      utl::getAttribute(sub,"v",nv);
      utl::getAttribute(sub,"w",nw);
      pch.uniformRefine(0,nu);
      pch.uniformRefine(1,nv);
      pch.uniformRefine(2,nw);
    }

  // Compute parameters
  const Go::SplineVolume* vol = pch.getBasis(ASM::GEOMETRY_BASIS);
  size_t px = vol->order(0)-1;
  size_t py = vol->order(1)-1;
  size_t pz = vol->order(2)-1;
  size_t nelemsx = vol->numCoefs(0) - px;
  size_t nelemsy = vol->numCoefs(1) - py;
  size_t nelemsz = vol->numCoefs(2) - pz;
  size_t nelemsx_sub = nelemsx / nx;
  size_t nelemsy_sub = nelemsy / ny;
  size_t nelemsz_sub = nelemsz / nz;
  size_t nelemsx_rem = nelemsx % nx;
  size_t nelemsy_rem = nelemsy % ny;
  size_t nelemsz_rem = nelemsz % nz;
  std::string header = "700 1 0 0\n";

  // Extract subpatches
  std::stringstream str;
  for (size_t k = 0; k < nz; ++k) {
    size_t nk = nelemsz_sub + (k < nelemsz_rem ? 1 : 0);
    size_t k0 = nk*k + (k < nelemsz_rem ? 0 : nelemsz_rem);
    size_t dk = nk + pz;
    for (size_t j = 0; j < ny; ++j) {
      size_t nj = nelemsy_sub + (j < nelemsy_rem ? 1 : 0);
      size_t j0 = nj*j + (j < nelemsy_rem ? 0 : nelemsy_rem);
      size_t dj = nj + py;
      for (size_t i = 0; i < nx; ++i) {
        size_t ni = nelemsx_sub + (i < nelemsx_rem ? 1 : 0);
        size_t i0 = ni*i + (i < nelemsx_rem ? 0 : nelemsx_rem);
        size_t di = ni + px;

        Go::SplineVolume subvol = getSubPatch(vol, {{i0, j0, k0}},
                                              {{di, dj, dk}},
                                              {{px+1, py+1, pz+1}});
        IFEM::cout << "  Number of knot spans in patch (" << i << ", " << j << ", " << k << "): "
                   << subvol.numCoefs(0)-subvol.order(0)+1 << "x"
                   << subvol.numCoefs(1)-subvol.order(1)+1 << "x"
                   << subvol.numCoefs(2)-subvol.order(2)+1 << std::endl;
        str << header << subvol;
      }
    }
  }

  return sim.readPatches(str,"\t");
}


Go::SplineVolume
MultiPatchModelGenerator3D::getSubPatch(const Go::SplineVolume* vol,
                                        const std::array<size_t,3>& start,
                                        const std::array<size_t,3>& numcoefs,
                                        const std::array<size_t,3>& order)
{
  bool rat = vol->rational();
  size_t I = vol->numCoefs(0);
  size_t J = vol->numCoefs(1);
  size_t N = vol->dimension() + (rat?1:0); // #coefs per control point
  std::vector<double>::const_iterator i0 = vol->basis(0).begin() + start[0];
  std::vector<double>::const_iterator j0 = vol->basis(1).begin() + start[1];
  std::vector<double>::const_iterator k0 = vol->basis(2).begin() + start[2];
  std::vector<double>::const_iterator n0 = rat ? vol->rcoefs_begin() : vol->coefs_begin();

  std::vector<double> subcoefs;
  for (size_t k = start[2]; k != start[2] + numcoefs[2]; k++)
    for (size_t j = start[1]; j != start[1] + numcoefs[1]; j++)
      for (size_t i = start[0]; i != start[0] + numcoefs[0]; i++)
        for (size_t n=0; n!=N; n++)
          subcoefs.push_back(*(n0 + n + N*i + N*I*j + N*I*J*k));

  return Go::SplineVolume(numcoefs[0], numcoefs[1], numcoefs[2],
                          order[0], order[1], order[2],
                          i0, j0, k0, subcoefs.begin(), vol->dimension(), rat);
}


bool MultiPatchModelGenerator3D::createTopology (SIMinput& sim) const
{
  if (!sim.createFEMmodel())
    return false;

  auto&& IJK = [this](int i, int j, int k) { return 1 + (k*ny+j)*nx + i; };

  int p1=2, p2=2, p3=2;
  if (subdivision)
    sim.getPatch(1)->getOrder(p1,p2,p3);

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
