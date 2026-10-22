// $Id$
//==============================================================================
//!
//! \file SIM3D.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based FEM analysis.
//!
//==============================================================================

#include "SIM3D.h"
#include "ModelGenerator.h"
#include "ASMs3D.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <map>
#include <utility>
#include <vector>


SIM3D::SIM3D (unsigned char n1, bool check)
{
  nf.push_back(n1);
  checkRHSys = check;
}


SIM3D::SIM3D (const CharVec& fields, bool check) : nf(fields)
{
  checkRHSys = check;
}


SIM3D::SIM3D (IntegrandBase* itg, unsigned char n, bool check) : SIMgeneric(itg)
{
  nf.push_back(n);
  checkRHSys = check;
}


bool SIM3D::connectPatches (const ASM::Interface& ifc, bool coordCheck)
{
  if (ifc.master == ifc.slave ||
      ifc.master < 1 || ifc.master > nGlPatches ||
      ifc.slave  < 1 || ifc.slave  > nGlPatches)
  {
    std::cerr <<" *** SIM3D::connectPatches: Invalid patch indices "
              << ifc.master <<" "<< ifc.slave << std::endl;
    return false;
  }

  if (ifc.orient < 0 || ifc.orient > 7)
  {
    std::cerr <<" *** SIM3D::connectPatches: Invalid orientation flag "
              << ifc.orient <<"."<< std::endl;
    return false;
  }

  int lmaster = this->getLocalPatchIndex(ifc.master);
  int lslave  = this->getLocalPatchIndex(ifc.slave);
  if (lmaster > 0 && lslave > 0)
  {
    if (ifc.dim < 2) return true; // ignored in serial

    IFEM::cout <<"\tConnecting P"<< ifc.slave <<" F"<< ifc.sidx
               <<" to P"<< ifc.master <<" F"<< ifc.midx
               <<" orient "<< ifc.orient << std::endl;

    ASM3D* spch = dynamic_cast<ASM3D*>(myModel[lslave-1]);
    ASM3D* mpch = dynamic_cast<ASM3D*>(myModel[lmaster-1]);
    if (spch && mpch)
    {
      std::set<int> bases;
      if (ifc.basis == 0)
      {
        // Only the continuous bases, since a basis that is discontinuous
        // across the interface has nothing to tie there
        for (size_t b = 1; b <= myModel[lslave-1]->getNoBasis(); b++)
          if (myModel[lslave-1]->isContinuousBasis(b))
            bases.insert(b);
      }
      else
        bases = utl::getDigits(ifc.basis);

      for (int b : bases)
        if (!spch->connectPatch(ifc.sidx,*mpch,ifc.midx,ifc.orient,b,
                                coordCheck,ifc.thick))
        {
          std::cerr <<" *** SIM3D::connectPatches: Failed to connect basis "
                    << b <<"."<< std::endl;
          return false;
        }
    }

    myInterfaces.push_back(ifc);
  }
  else
  {
    // The basis selection above does not reach the distributed topology:
    // DomainDecomposition reads a zero basis as every basis, so an interface
    // handed over with one would tie a discontinuous basis across the ranks
    // and leave the space, and the answer, depending on how the model was
    // partitioned. Until the selection is carried through, a model with such
    // a basis is refused here rather than solved differently on every
    // partition.
    const ASMbase* pch = lslave > 0 ? myModel[lslave-1]
                       : lmaster > 0 ? myModel[lmaster-1] : nullptr;
    if (pch)
      for (size_t b = 1; b <= pch->getNoBasis(); b++)
      {
        // A basis whose degrees of freedom follow the parametrization is
        // tied across a turned face with a sign, and the distributed
        // topology equates node numbers and carries no sign, so it cannot
        // say what this interface needs. Unlike the selection below, this
        // holds whether the bases were chosen here or named in the
        // interface, so it is asked of both.
        if (ifc.orient != 0 && pch->dofsFollowParametrization(b) &&
            (ifc.basis == 0 || utl::getDigits(ifc.basis).count(b)))
        {
          std::cerr <<" *** SIM3D::connectPatches: The degrees of freedom of"
                    <<" basis "<< b <<" follow the\n     parametrization,"
                    <<" and the face between P"<< ifc.master <<" and P"
                    << ifc.slave <<" is turned over\n     and shared with"
                    <<" another process, whose topology carries no sign."
                    <<" Run this\n     model on one process."<< std::endl;
          return false;
        }

        // A basis this interface leaves untied, and which is meant to
        // stay that way, is two things the distributed topology cannot
        // carry. The selection itself is one: a zero reaches
        // DomainDecomposition as every basis, so it would tie the basis
        // after all. The other is what is left undetermined by not tying
        // it -- the mode around a loop of patches, which connectCrossPoints
        // walks for and ties. That walk reads the interfaces of this
        // process, and an interface shared with another is not among them,
        // so a loop which crosses the partition is neither seen nor tied.
        if (!pch->isContinuousBasis(b) &&
            (ifc.basis == 0 || !utl::getDigits(ifc.basis).count(b)))
        {
          std::cerr <<" *** SIM3D::connectPatches: Basis "<< b <<" is"
                    <<" discontinuous across an interface,\n     and the "
                    "face between P"<< ifc.master <<" and P"<< ifc.slave
                    <<" is shared with another\n     process. Neither which"
                    <<" bases are tied nor the mode left untied where"
                    <<"\n     patches close a loop survives into the"
                    <<" topology of a run on several\n     processes. Run"
                    <<" this model on one."<< std::endl;
          return false;
        }
      }

    adm.dd.ghostConnections.insert(ifc);
  }

  return true;
}


/*!
  Where more than two patches meet around a common edge, a basis which is not
  tied across the interfaces has a mode that alternates in sign around the
  edge, exactly as it does around a cross point in two dimensions -- extruding
  the four-patch grid SIM2D::connectCrossPoints is written for produces four
  patches around an edge and the very same mode. Tying the basis between two
  of the patches meeting there removes it, and that is what the two
  dimensional case does; nothing here does it yet.

  Rather than hand the solver a model with an undetermined mode in it, such a
  model is refused. The patches meeting around an edge are found from the
  geometry, by the midpoint of each of the twelve edges of each patch: an edge
  shared by more than two patches is a cross edge.
*/

bool SIM3D::connectCrossPoints ()
{
  // The bases which are left discontinuous at the interfaces
  std::vector<size_t> dBases;
  for (const ASMbase* pch : myModel)
    if (!pch->empty())
    {
      for (size_t b = 1; b <= pch->getNoBasis(); b++)
        if (!pch->isContinuousBasis(b))
          dBases.push_back(b);
      break;
    }

  if (dBases.empty() || myInterfaces.empty())
    return true;

  // The two ends of each edge, and which patch it belongs to. Two patches
  // are connected when their coordinates agree to the tolerance the
  // connection itself accepts, so the ends are matched by that same
  // tolerance rather than rounded onto a grid: a grid separates two points
  // which fall either side of a line between cells however close they are,
  // and would then miss the very edge it is looking for. There are twelve
  // edges to a patch, so comparing them all against each other is nothing.
  //
  // It takes both ends to say which edge this is. A midpoint alone is the
  // same for every edge through it, so two edges of different patches
  // crossing at their middles -- which a symmetric model has plenty of --
  // would be counted as one and the model refused for a cross edge it does
  // not have.
  //! \brief A patch meeting an edge, and the two faces of it the edge is on.
  struct Incident { size_t patch; int face[2]; };

  const double xtol = 1.0e-4;
  std::vector<std::pair<std::array<Vec3,2>,std::vector<Incident>>> edges;
  for (size_t p = 0; p < myModel.size(); p++)
  {
    const ASM3D* pch = dynamic_cast<const ASM3D*>(myModel[p]);
    if (!pch || myModel[p]->empty())
      continue;

    for (int dir = 0; dir < 3; dir++)
      for (int a = -1; a <= 1; a += 2)
        for (int b = -1; b <= 1; b += 2)
        {
          // The two ends of one edge: the direction it runs along takes
          // both values, the other two are fixed at a corner
          int ijk[2][3];
          for (int e = 0; e < 2; e++)
          {
            ijk[e][dir] = e ? 1 : -1;
            ijk[e][(dir+1)%3] = a;
            ijk[e][(dir+2)%3] = b;
          }

          const int n1 = pch->getCorner(ijk[0][0],ijk[0][1],ijk[0][2],1);
          const int n2 = pch->getCorner(ijk[1][0],ijk[1][1],ijk[1][2],1);
          if (n1 < 1 || n2 < 1)
            continue;

          const Vec3 X1 = myModel[p]->getCoord(n1);
          const Vec3 X2 = myModel[p]->getCoord(n2);

          // The two patches sharing an edge may run it either way round
          auto&& same = [&X1,&X2,xtol](const std::array<Vec3,2>& known)
          {
            return (known[0].equal(X1,xtol) && known[1].equal(X2,xtol)) ||
                   (known[0].equal(X2,xtol) && known[1].equal(X1,xtol));
          };

          // The faces this edge lies on: the one across each of the two
          // directions it does not run along, at the side it sits at.
          // Faces are numbered umin, umax, vmin, vmax, wmin, wmax.
          const Incident here = {p, {2*((dir+1)%3) + (a < 0 ? 1 : 2),
                                     2*((dir+2)%3) + (b < 0 ? 1 : 2)}};

          auto it = std::find_if(edges.begin(), edges.end(),
                                 [&same](const auto& known)
                                 { return same(known.first); });
          if (it == edges.end())
            edges.push_back({{X1,X2},{here}});
          else
            it->second.push_back(here);
        }
  }

  for (const auto& known : edges)
  {
    // Named rather than bound structurally: the lambdas below
    // capture this, and capturing a structured binding is C++20
    const std::vector<Incident>& incident = known.second;

    // Two patches meeting is the interior of an interface, and there the
    // discontinuous bases are meant to stay discontinuous
    if (incident.size() < 3)
      continue;

    // Walk the interfaces meeting along this edge while keeping track of
    // which patches they have joined up, as SIM2D::connectCrossPoints does
    // at a point. It is around a loop that the alternating mode lives, and
    // a loop is closed by an interface joining two patches which are
    // connected already. Patches which merely touch along the edge, or
    // which are joined in a chain, carry no such mode and are left alone.
    std::map<size_t,size_t> group;
    for (const Incident& i : incident)
      group[i.patch] = i.patch;

    auto&& onEdge = [&incident](int lpatch, int face)
    {
      for (const Incident& i : incident)
        if (i.patch + 1 == static_cast<size_t>(lpatch) &&
            (i.face[0] == face || i.face[1] == face))
          return true;
      return false;
    };

    // One basis at a time. An interface which names the bases it ties can
    // name a discontinuous one, and connectPatches then connects it: that
    // basis is tied along that interface and the loop is broken for it,
    // whatever the others do. So the walk is per basis, over the
    // interfaces which leave that basis free.
    for (size_t dBasis : dBases)
    {
      std::map<size_t,size_t> group0(group);
      auto&& root = [&group0](size_t i)
      {
        while (group0[i] != i) i = group0[i];
        return i;
      };

      for (const ASM::Interface& ifc : myInterfaces)
      {
        if (ifc.dim < 2) continue; // joins less than a whole face

        // Named bases are tied as named; a zero leaves out the ones which
        // are discontinuous, which is this one
        if (ifc.basis != 0 && utl::getDigits(ifc.basis).count(dBasis))
          continue;

        const int lmaster = this->getLocalPatchIndex(ifc.master);
        const int lslave  = this->getLocalPatchIndex(ifc.slave);
        if (lmaster < 1 || lslave < 1)
          continue;
        if (!onEdge(lmaster,ifc.midx) || !onEdge(lslave,ifc.sidx))
          continue;

        const size_t rm = root(lmaster-1);
        const size_t rs = root(lslave-1);
        if (rm != rs)
        {
          group0[rm] = rs;
          continue;
        }

        std::cerr <<" *** SIM3D::connectCrossPoints: "<< incident.size()
                  <<" patches close a loop around one edge,\n     and basis "
                  << dBasis <<" is not tied across their interfaces,"
                  <<" which leaves a mode\n     alternating in sign around"
                  <<" that edge with nothing to determine it. Tying it\n"
                  <<"     there is what SIM2D::connectCrossPoints does in two"
                  <<" dimensions, and it is\n     not written for three."
                  << std::endl;
        return false;
      }
    }
  }

  return true;
}


bool SIM3D::parseGeometryTag (const tinyxml2::XMLElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"refine") && !isRefined)
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0, addv = 0, addw = 0;
      utl::getAttribute(elem,"u",addu);
      utl::getAttribute(elem,"v",addv);
      utl::getAttribute(elem,"w",addw);
      for (int j : patches)
      {
        IFEM::cout <<"\tRefining P"<< j
                   <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
        ASM3D* pch = dynamic_cast<ASM3D*>(this->getPatch(j,true));
        if (pch)
        {
          pch->uniformRefine(0,addu);
          pch->uniformRefine(1,addv);
          pch->uniformRefine(2,addw);
        }
      }
    }
    else
    {
      // Non-uniform (graded) refinement
      int dir = 1;
      utl::getAttribute(elem,"dir",dir);
      const char* refdata = elem->FirstChild()->Value();
      for (int j : patches)
      {
        IFEM::cout <<"\tRefining P"<< j <<" dir="<< dir;
        if (refdata && isalpha(refdata[0]))
          IFEM::cout <<" with grading "<< refdata <<":";
        for (size_t i = 0; i < xi.size(); i++)
          IFEM::cout << (i%10 || xi.size() < 11 ? " " : "\n\t") << xi[i];
        IFEM::cout << std::endl;
        ASM3D* pch = dynamic_cast<ASM3D*>(this->getPatch(j,true));
        if (pch) pch->refine(dir-1,xi);
      }
    }
  }

  else if (!strcasecmp(elem->Value(),"raiseorder") && !isRefined)
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    bool setOrder = false;
    int addu = 0, addv = 0, addw = 0;
    utl::getAttribute(elem,"u",addu);
    utl::getAttribute(elem,"v",addv);
    utl::getAttribute(elem,"w",addw);
    utl::getAttribute(elem,"setTo",setOrder);
    for (int j : patches)
    {
      IFEM::cout << (setOrder ? "\tSetting":"\tRaising") <<" order of P"<< j
                 <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
      ASM3D* pch = dynamic_cast<ASM3D*>(this->getPatch(j,true));
      if (pch) pch->raiseOrder(addu,addv,addw,setOrder);
    }
  }

  else if (!strcasecmp(elem->Value(),"topology"))
  {
    if (!this->createFEMmodel()) return false;

    int offset = 0;
    utl::getAttribute(elem,"offset",offset);

    const tinyxml2::XMLElement* child = elem->FirstChildElement("connection");
    for (; child; child = child->NextSiblingElement("connection"))
    {
      ASM::Interface ifc;
      bool periodic = false;
      if (utl::getAttribute(child,"master",ifc.master))
        ifc.master += offset;
      if (!utl::getAttribute(child,"midx",ifc.midx))
        utl::getAttribute(child,"mface",ifc.midx);
      if (utl::getAttribute(child,"slave",ifc.slave))
        ifc.slave += offset;
      if (!utl::getAttribute(child,"sidx",ifc.sidx))
        utl::getAttribute(child,"sface",ifc.sidx);
      utl::getAttribute(child,"orient",ifc.orient);
      utl::getAttribute(child,"basis",ifc.basis);
      utl::getAttribute(child,"periodic",periodic);
      if (!utl::getAttribute(child,"dim",ifc.dim))
        ifc.dim = 2;

      if (!this->connectPatches(ifc,!periodic))
        return false;
    }
  }

  else if (!strcasecmp(elem->Value(),"periodic"))
    return this->parsePeriodic(elem);

  else if (!strcasecmp(elem->Value(),"collapse"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, face = 1, edge = 0;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"face",face);
    utl::getAttribute(elem,"edge",edge);

    if (patch < 1 || patch > nGlPatches)
    {
      std::cerr <<" *** SIM3D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }

    IFEM::cout <<"\tCollapsed face P"<< patch <<" F"<< face;
    if (edge > 0) IFEM::cout <<" on to edge "<< edge;
    IFEM::cout << std::endl;
    ASMs3D* pch = dynamic_cast<ASMs3D*>(this->getPatch(patch,true));
    if (pch) return pch->collapseFace(face,edge);
  }

  else if (!strcasecmp(elem->Value(),"projection") && !isRefined)
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    if (child && !strncasecmp(child->Value(),"patch",5) && child->FirstChild())
    {
      // Read projection basis from file
      const char* patch = child->FirstChild()->Value();
      std::istream* isp = getPatchStream(child->Value(),patch);
      if (!isp) return false;

      for (ASMbase* pch : myModel)
        pch->createProjectionBasis(false);

      bool ok = true;
      for (int pid = 1; isp->good() && ok; pid++)
      {
        IFEM::cout <<"\tReading projection basis for patch "<< pid << std::endl;
        ASMbase* pch = this->getPatch(pid,true);
        if (pch)
          ok = pch->read(*isp);
        else if ((pch = ASM3D::create(opt.discretization,nf)))
        {
          // Skip this patch
          ok = pch->read(*isp);
          delete pch;
        }
      }

      delete isp;
      if (!ok) return false;
    }
    else // Generate separate projection basis from current geometry basis
      for (ASMbase* pch : myModel)
        pch->createProjectionBasis(true);

    // Apply refine and/raise-order commands to define the projection basis
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"refine") ||
          !strcasecmp(child->Value(),"raiseorder"))
        if (!this->parseGeometryTag(child))
          return false;

    for (ASMbase* pch : myModel)
      if (!pch->createProjectionBasis(false))
      {
        std::cerr <<" *** SIM3D::parseGeometryTag: Failed to create projection"
                  <<" basis, check patch file specification."<< std::endl;
        return false;
      }
  }

  return true;
}


bool SIM3D::parseBCTag (const tinyxml2::XMLElement* elem)
{
  if (ignoreDirichlet)
    return true; // Ignore all boundary conditions

  else if (!strcasecmp(elem->Value(),"fixpoint"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0;
    utl::getAttribute(elem,"patch",patch);
    int pid = this->getLocalPatchIndex(patch);
    if (pid < 1) return pid == 0;

    ASM3D* pch = dynamic_cast<ASM3D*>(myModel[pid-1]);
    if (!pch) return false;

    int code = 123;
    double rx = 0.0, ry = 0.0, rz = 0.0;
    utl::getAttribute(elem,"code",code);
    utl::getAttribute(elem,"rx",rx);
    utl::getAttribute(elem,"ry",ry);
    utl::getAttribute(elem,"rz",rz);
    IFEM::cout <<"\tConstraining P"<< patch
               <<" point at "<< rx <<" "<< ry <<" "<< rz
               <<" with code "<< code << std::endl;
    pch->constrainNode(rx,ry,rz,code);
  }

  else if (!strcasecmp(elem->Value(),"fixline"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0;
    utl::getAttribute(elem,"patch",patch);
    if ((patch = this->getLocalPatchIndex(patch)) < 1)
      return patch == 0;

    int line = 0, code = 123;
    double xi = 0.0;
    utl::getAttribute(elem,"line",line);
    utl::getAttribute(elem,"code",code);
    utl::getAttribute(elem,"xi",xi);
    return this->addLineConstraint(patch,line%10,line/10,xi,code);
  }

  return true;
}


bool SIM3D::parse (const tinyxml2::XMLElement* elem)
{
  bool result = this->SIMgeneric::parse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"geometry"))
      result &= this->parseGeometryTag(child);
    else if (!strcasecmp(elem->Value(),"boundaryconditions"))
      result &= this->parseBCTag(child);

  if (myGen && result)
    result = myGen->createTopology(*this);

  delete myGen;
  myGen = nullptr;

  return result;
}


bool SIM3D::parse (char* keyWord, std::istream& is)
{
  char* cline = nullptr;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM3D* pch = nullptr;
      IFEM::cout <<"\nNumber of patch refinements: "<< nref << std::endl;
      for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
      {
        bool uniform = !strchr(cline,'.');
        int patch = atoi(strtok(cline," "));
        if (patch == 0 || abs(patch) > (int)myModel.size())
        {
          std::cerr <<" *** SIM3D::parse: Invalid patch index "
                    << patch << std::endl;
          return false;
        }
        int ipatch = patch-1;
        if (patch < 0)
        {
          ipatch = 0;
          patch = -patch;
        }
        if (uniform)
        {
          int addu = atoi(strtok(nullptr," "));
          int addv = atoi(strtok(nullptr," "));
          int addw = atoi(strtok(nullptr," "));
          for (int j = ipatch; j < patch; j++)
            if ((pch = dynamic_cast<ASM3D*>(myModel[j])))
            {
              IFEM::cout <<"\tRefining P"<< j+1
                         <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
              pch->uniformRefine(0,addu);
              pch->uniformRefine(1,addv);
              pch->uniformRefine(2,addw);
            }
        }
        else
        {
          RealArray xi;
          int dir = atoi(strtok(nullptr," "));
          if (utl::parseKnots(xi))
            for (int j = ipatch; j < patch; j++)
              if ((pch = dynamic_cast<ASM3D*>(myModel[j])))
              {
                IFEM::cout <<"\tRefining P"<< j+1 <<" dir="<< dir;
                for (double u : xi) IFEM::cout <<" "<< u;
                IFEM::cout << std::endl;
                pch->refine(dir-1,xi);
              }
        }
      }
    }
  }

  else if (!strncasecmp(keyWord,"RAISEORDER",10))
  {
    int nref = atoi(keyWord+10);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM3D* pch = nullptr;
      IFEM::cout <<"\nNumber of order raise: "<< nref << std::endl;
      for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
      {
        int patch = atoi(strtok(cline," "));
        int addu  = atoi(strtok(nullptr," "));
        int addv  = atoi(strtok(nullptr," "));
        int addw  = atoi(strtok(nullptr," "));
        if (patch == 0 || abs(patch) > (int)myModel.size())
        {
          std::cerr <<" *** SIM3D::parse: Invalid patch index "
                    << patch << std::endl;
          return false;
        }
        int ipatch = patch-1;
        if (patch < 0)
        {
          ipatch = 0;
          patch = -patch;
        }
        for (int j = ipatch; j < patch; j++)
          if ((pch = dynamic_cast<ASM3D*>(myModel[j])))
          {
            IFEM::cout <<"\tRaising order of P"<< j+1
                       <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
            pch->raiseOrder(addu,addv,addw);
          }
      }
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGYFILE",12))
  {
    if (!this->createFEMmodel()) return false;

    size_t i = 12; while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::ifstream ist(keyWord+i);
    if (ist)
      IFEM::cout <<"\nReading data file "<< keyWord+i << std::endl;
    else
    {
      std::cerr <<" *** SIM3D::parse: Failure opening input file "
                << std::string(keyWord+i) << std::endl;
      return false;
    }

    while ((cline = utl::readLine(ist)))
    {
      ASM::Interface ifc;
      ifc.master = atoi(strtok(cline," "));
      ifc.midx   = atoi(strtok(nullptr," "));
      ifc.slave  = atoi(strtok(nullptr," "));
      ifc.sidx   = atoi(strtok(nullptr," "));
      int swapd  = atoi(strtok(nullptr," "));
      int rev_u  = atoi(strtok(nullptr," "));
      int rev_v  = atoi(strtok(nullptr," "));
      ifc.orient = 4*swapd+2*rev_u+rev_v;
      ifc.dim    = 2;
      if (!this->connectPatches(ifc))
        return false;
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGY",8))
  {
    if (!this->createFEMmodel()) return false;

    int ntop = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
    for (int i = 0; i < ntop && (cline = utl::readLine(is)); i++)
    {
      ASM::Interface ifc;
      ifc.master = atoi(strtok(cline," "));
      ifc.midx   = atoi(strtok(nullptr," "));
      ifc.slave  = atoi(strtok(nullptr," "));
      ifc.sidx   = atoi(strtok(nullptr," "));
      ifc.orient = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      ifc.dim    = 2;
      if (!this->connectPatches(ifc))
        return false;
    }
  }

  else if (!strncasecmp(keyWord,"CONSTRAINTS",11))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int ngno = 0;
    int ncon = atoi(keyWord+11);
    IFEM::cout <<"\nNumber of constraints: "<< ncon << std::endl;
    for (int i = 0; i < ncon && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pface = atoi(strtok(nullptr," "));
      int bcode = atoi(strtok(nullptr," "));
      double pd = (cline = strtok(nullptr," ")) ? atof(cline) : 0.0;

      patch = this->getLocalPatchIndex(patch);
      if (patch < 1) continue;

      int ldim = pface < 0 ? 0 : 2;
      if (pface < 0) pface = -pface;

      if (pface > 10)
      {
        if (!this->addLineConstraint(patch,pface%10,pface/10,pd,bcode))
          return false;
      }
      else if (pd == 0.0)
      {
        if (!this->addConstraint(patch,pface,ldim,bcode%1000000,0,ngno))
          return false;
      }
      else
      {
        int code = 1000000 + bcode;
        while (myScalars.find(code) != myScalars.end())
          code += 1000000;

        if (!this->addConstraint(patch,pface,ldim,bcode%1000000,-code,ngno))
          return false;

        IFEM::cout <<" ";
        cline = strtok(nullptr," ");
        myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,pd));
      }
      if (pface < 10) IFEM::cout << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"FIXPOINTS",9))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int nfix = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
    for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      double rx = atof(strtok(nullptr," "));
      double ry = atof(strtok(nullptr," "));
      double rz = atof(strtok(nullptr," "));
      int bcode = (cline = strtok(nullptr," ")) ? atoi(cline) : 123;

      ASM3D* pch = dynamic_cast<ASM3D*>(this->getPatch(patch,true));
      if (pch)
      {
        IFEM::cout <<"\tConstraining P"<< patch
                   <<" point at "<< rx <<" "<< ry <<" "<< rz
                   <<" with code "<< bcode << std::endl;
        pch->constrainNode(rx,ry,rz,bcode);
      }
    }
  }

  else
    return this->SIMgeneric::parse(keyWord,is);

  return true;
}


bool SIM3D::addConstraint (int patch, int lndx, int ldim, int dirs, int code,
                           int& ngnod, char basis, bool ovrD)
{
  // Lambda function for error message generation
  auto&& error = [](const char* message, int idx, bool coutNL = false)
  {
    if (coutNL) IFEM::cout << std::endl;
    std::cerr <<" *** SIM3D::addConstraint: Invalid "
              << message <<" ("<< idx <<")."<< std::endl;
    return false;
  };

  if (patch < 1 || patch > (int)myModel.size())
    return error("patch index",patch);

  int aldim = abs(ldim);
  if (this->SIMbase::addConstraint(patch,lndx,aldim,dirs,code,ngnod,basis,ovrD))
    return true;

  bool open = ldim < 0; // open means without face edges or edge ends
  bool project = lndx < -10;
  if (project) lndx += 10;
  if (lndx < 0 && aldim > 3) aldim = 2; // local tangent direction is indicated

  // Must dynamic cast here, since ASM3D is not derived from ASMbase
  ASM3D* pch = dynamic_cast<ASM3D*>(myModel[patch-1]);
  if (!pch) return error("3D patch",myModel[patch-1]->idx+1);

  switch (aldim)
    {
    case 0: // Vertex constraints
      switch (lndx)
        {
        case 1: pch->constrainCorner(-1,-1,-1,dirs,abs(code),basis); break;
        case 2: pch->constrainCorner( 1,-1,-1,dirs,abs(code),basis); break;
        case 3: pch->constrainCorner(-1, 1,-1,dirs,abs(code),basis); break;
        case 4: pch->constrainCorner( 1, 1,-1,dirs,abs(code),basis); break;
        case 5: pch->constrainCorner(-1,-1, 1,dirs,abs(code),basis); break;
        case 6: pch->constrainCorner( 1,-1, 1,dirs,abs(code),basis); break;
        case 7: pch->constrainCorner(-1, 1, 1,dirs,abs(code),basis); break;
        case 8: pch->constrainCorner( 1, 1, 1,dirs,abs(code),basis); break;
        default: return error("vertex index",lndx,true);
        }
      break;

    case 1: // Edge constraints
      if (lndx > 0 && lndx <= 12)
        pch->constrainEdge(lndx,open,dirs,code,basis);
      else
        return error("edge index",lndx,true);
      break;

    case 2: // Face constraints
      switch (lndx)
        {
        case  1: pch->constrainFace(-1,open,dirs,code,basis); break;
        case  2: pch->constrainFace( 1,open,dirs,code,basis); break;
        case  3: pch->constrainFace(-2,open,dirs,code,basis); break;
        case  4: pch->constrainFace( 2,open,dirs,code,basis); break;
        case  5: pch->constrainFace(-3,open,dirs,code,basis); break;
        case  6: pch->constrainFace( 3,open,dirs,code,basis); break;
        case -1:
          ngnod += pch->constrainFaceLocal(-1,open,dirs,code,project,ldim);
          break;
        case -2:
          ngnod += pch->constrainFaceLocal( 1,open,dirs,code,project,ldim);
          break;
        case -3:
          ngnod += pch->constrainFaceLocal(-2,open,dirs,code,project,ldim);
          break;
        case -4:
          ngnod += pch->constrainFaceLocal( 2,open,dirs,code,project,ldim);
          break;
        case -5:
          ngnod += pch->constrainFaceLocal(-3,open,dirs,code,project,ldim);
          break;
        case -6:
          ngnod += pch->constrainFaceLocal( 3,open,dirs,code,project,ldim);
          break;
        default:
          return error("face index",lndx,true);
        }
      break;

    default:
      return error("local dimension switch",ldim,true);
    }

  return true;
}


bool SIM3D::addLineConstraint (int patch, int lndx, int line, double xi,
                               int dirs, char basis)
{
  // Lambda function for error message generation
  auto&& error = [](const char* message, int idx)
  {
    std::cerr <<" *** SIM3D::addLineConstraint: Invalid "
              << message <<" ("<< idx <<")."<< std::endl;
    return false;
  };

  if (patch < 1 || patch > (int)myModel.size())
    return error("patch index",patch);

  IFEM::cout <<"\tConstraining P"<< myModel[patch-1]->idx+1
             <<" F"<< lndx <<" L"<< line <<" at xi="<< xi
             <<" in direction(s) "<< dirs
             <<" basis = " << (int)basis << std::endl;

  // Must dynamic cast here, since ASM3D is not derived from ASMbase
  ASM3D* pch = dynamic_cast<ASM3D*>(myModel[patch-1]);
  if (!pch) return error("3D patch",myModel[patch-1]->idx+1);

  switch (line)
    {
    case 1: // Face line constraints in local I-direction
      switch (lndx)
        {
        case 1: pch->constrainLine(-1,2,xi,dirs,0,basis); break;
        case 2: pch->constrainLine( 1,2,xi,dirs,0,basis); break;
        case 3: pch->constrainLine(-2,3,xi,dirs,0,basis); break;
        case 4: pch->constrainLine( 2,3,xi,dirs,0,basis); break;
        case 5: pch->constrainLine(-3,1,xi,dirs,0,basis); break;
        case 6: pch->constrainLine( 3,1,xi,dirs,0,basis); break;
        default: return error("face index",lndx);
        }
      break;

    case 2: // Face line constraints in local J-direction
      switch (lndx)
        {
        case 1: pch->constrainLine(-1,3,xi,dirs,0,basis); break;
        case 2: pch->constrainLine( 1,3,xi,dirs,0,basis); break;
        case 3: pch->constrainLine(-2,1,xi,dirs,0,basis); break;
        case 4: pch->constrainLine( 2,1,xi,dirs,0,basis); break;
        case 5: pch->constrainLine(-3,2,xi,dirs,0,basis); break;
        case 6: pch->constrainLine( 3,2,xi,dirs,0,basis); break;
        default: return error("face index",lndx);
        }
      break;

    default:
      return error("face line index",line);
    }

  return true;
}


ASMbase* SIM3D::readPatch (std::istream& isp, int pchInd, const CharVec& unf,
                           const char* whiteSpace) const
{
  const CharVec& uunf = unf.empty() ? nf : unf;
  bool isMixed = uunf.size() > 1 && uunf[1] > 0;
  ASMbase* pch = ASM3D::create(opt.discretization,uunf,isMixed);
  if (pch)
  {
    if (!pch->read(isp) || this->getLocalPatchIndex(pchInd+1) < 1)
    {
      delete pch;
      pch = nullptr;
    }
    else
    {
      if (whiteSpace)
        IFEM::cout << whiteSpace <<"Reading patch "<< pchInd+1 << std::endl;
      if (checkRHSys && dynamic_cast<ASM3D*>(pch)->checkRightHandSystem())
        IFEM::cout <<"\tSwapped."<< std::endl;
      pch->idx = myModel.size();
    }
  }

  return pch;
}


void SIM3D::readNodes (std::istream& isn)
{
  while (isn.good())
  {
    int patch = 0;
    isn >> patch;
    int pid = this->getLocalPatchIndex(patch+1);
    if (pid < 0) return;

    if (!this->readNodes(isn,pid-1))
    {
      std::cerr <<" *** SIM3D::readNodes: Failed to assign node numbers"
                <<" for patch "<< patch+1 << std::endl;
      return;
    }
  }
}


bool SIM3D::readNodes (std::istream& isn, int pchInd, int basis, bool oneBased)
{
  int i;
  ASMs3D::BlockNodes n;

  for (i = 0; i <  8 && isn.good(); i++)
    isn >> n.ibnod[i];
  for (i = 0; i < 12 && isn.good(); i++)
    isn >> n.edges[i].icnod >> n.edges[i].incr;
  for (i = 0; i <  6 && isn.good(); i++)
    isn >> n.faces[i].isnod >> n.faces[i].incrI >> n.faces[i].incrJ;
  isn >> n.iinod;

  if (!isn.good() || pchInd < 0) return true;

  if (!oneBased)
  {
    // We always require the node numbers to be 1-based
    for (i = 0; i <  8; i++) ++n.ibnod[i];
    for (i = 0; i < 12; i++) ++n.edges[i].icnod;
    for (i = 0; i <  6; i++) ++n.faces[i].isnod;
    ++n.iinod;
  }

  return static_cast<ASMs3D*>(myModel[pchInd])->assignNodeNumbers(n,basis);
}


void SIM3D::clonePatches (const PatchVec& patches,
                          const std::map<int,int>& glb2locN)
{
  for (ASMbase* patch : patches)
  {
    ASM3D* pch3D = dynamic_cast<ASM3D*>(patch);
    if (pch3D) myModel.push_back(pch3D->clone(nf));
  }

  g2l = &glb2locN;

  if (nGlPatches == 0)
    nGlPatches = myModel.size();
}


ModelGenerator* SIM3D::getModelGenerator (const tinyxml2::XMLElement* geo) const
{
  return new DefaultGeometry3D(geo);
}


RealArray SIM3D::getSolution (const Vector& psol, double u, double v, double w,
                              int deriv, int patch) const
{
  double par[3] = { u, v, w };
  return this->SIMgeneric::getSolution(psol,par,deriv,patch);
}
