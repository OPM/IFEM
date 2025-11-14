// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D %Lagrange FE models.
//!
//==============================================================================

#include "ASMu2DLag.h"
#include "ElementBlock.h"
#include "GlobalIntegral.h"
#include "Integrand.h"
#include "MPC.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include <numeric>
#include <sstream>
#include <fstream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

extern "C" {
  //! \brief Parametric inversion for a linear triangle.
  void t3_inv_(const int& nsd, const double& tol, const double& eps,
               const double* X, const double* Xn, double* xi,
               const int& ipsw, const int& iwr, int& ierr);
  //! \brief Parametric inversion for a bilinear quadrilateral.
  void q4_inv_(const int& nsd, const double& tol, const double& eps,
               const double* X, const double* Xn, double* xi,
               const int& ipsw, const int& iwr, int& ierr);
  //! \brief Evaluate the bilinear shape functions.
  void shape2_(const int& iop, const int& nenod, const double* xi,
               double* shp, const int& ipsw, const int& iwr, int& ierr);
  //! \brief Opens a temporary file for Fortran print.
  int openftnfile_(const char* fname, const ssize_t nchar);
  //! \brief Closes a Fortran file.
  void closeftnfile_(const int& iunit);
}


ASMu2DLag::ASMu2DLag (unsigned char n_s,
                      unsigned char n_f, char fType) : ASMs2DLag(n_s,n_f)
{
  fileType = fType;
  swapNode34 = false;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p, unsigned char n_f) :
  ASMs2DLag(p,n_f), nodeSets(p.nodeSets)
{
  fileType = 0;
  swapNode34 = p.swapNode34;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p) :
  ASMs2DLag(p), nodeSets(p.nodeSets)
{
  fileType = 0;
  swapNode34 = p.swapNode34;
}


bool ASMu2DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    swapNode34 = true;
    return ASM::readMatlab(is,myMNPC,myCoord,nodeSets);
  case 'x':
  case 'X':
    return ASM::readXML(is,myMNPC,myCoord,nodeSets);
  default:
    std::cerr <<" *** ASMu2DLag::read: Undefined file format."<< std::endl;
    return false;
  }
}


bool ASMu2DLag::generateFEMTopology ()
{
  p1 = p2 = 2; // So far only linear elements supported

  firstEl = gEl;

  nnod = myCoord.size();
  nel  = myMNPC.size();

  bool ok = true;
  if (myMLGN.empty())
  {
    myMLGN.resize(nnod);
    std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  }
  else
    ok = myMLGN.size() == nnod;

  if (myMLGE.empty())
  {
    myMLGE.resize(nel);
    std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);
  }
  else if (ok)
    ok = myMLGE.size() == nel;

  if (!ok)
    std::cerr <<" *** ASMu2DLag::generateFEMTopology: Array mismatch, "
              <<" size(coord)="<< myCoord.size() <<" size(MLGN)="<< MLGN.size()
              <<" size(MNPC)="<< MNPC.size() <<" size(MLGE)="<< MLGE.size()
              << std::endl;

  gNod += nnod;
  gEl  += nel;

  myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  return ok;
}


int ASMu2DLag::getNodeSetIdx (const std::string& setName) const
{
  int iset = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMu2DLag::getNodeSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(nodeSets.size()))
    return nodeSets[iset-1].second;

  return this->ASMbase::getNodeSet(iset);
}


bool ASMu2DLag::getNodeSet (int iset, std::string& name) const
{
  if (iset < 0 || iset > static_cast<int>(nodeSets.size()))
    return false;
  else if (nodeSets[iset-1].second.empty())
    return false;

  name = nodeSets[iset-1].first;
  return true;
}


/*!
  If \a inod is negative, the absolute value is taken as the external node ID.
  Otherwise, it is taken as the 1-based internal node index within the patch.
*/

bool ASMu2DLag::isInNodeSet (int iset, int inod) const
{
  if (iset < 1 || iset > static_cast<int>(nodeSets.size()))
    return false;

  if (inod < 0)
    inod = this->getNodeIndex(-inod);

  return utl::findIndex(nodeSets[iset-1].second,inod) >= 0;
}


int ASMu2DLag::parseNodeSet (const std::string& setName, const char* cset)
{
  int iset = this->getNodeSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = nodeSets[iset].second;
  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  // Transform to internal node indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if (mySet[i] > 0)
    {
      if (int inod = this->getNodeIndex(mySet[i]); inod > 0)
        mySet[ifirst++] = inod;
      else
        IFEM::cout <<"  ** Warning: Non-existing node "<< mySet[i]
                   <<" in the set \""<< setName <<"\" (ignored)."<< std::endl;
    }

  if (ifirst < mySet.size())
    mySet.resize(ifirst);

  return 1+iset;
}


int ASMu2DLag::parseNodeBox (const std::string& setName, const char* data)
{
  if (myCoord.empty())
    return 0; // No nodes yet

  Vec3 X0, X1;
  std::istringstream(data) >> X0 >> X1;

  // Lambda function for checking if a point is within the bounding box
  auto&& isInside=[&X0,&X1](const Vec3& X)
  {
    for (int i = 0; i < 3; i++)
      if (X[i] < X0[i] || X[i] > X1[i])
        return false;
    return true;
  };

  IntVec nodes;
  for (size_t inod = 0; inod < myCoord.size(); inod++)
    if (isInside(myCoord[inod]))
      nodes.push_back(1+inod);

  IFEM::cout <<"\tBounding Box: "<< X0 <<" - "<< X1
             <<": "<< nodes.size() <<" nodes"<< std::endl;

  if (nodes.empty())
    return 0; // No nodes are within the given box

  for (size_t iset = 0; iset < nodeSets.size(); iset++)
    if (nodeSets[iset].first == setName)
    {
      if (nodeSets[iset].second.empty())
        nodeSets[iset].second.swap(nodes);
      else for (int inod : nodes)
        if (utl::findIndex(nodeSets[iset].second,inod) < 0)
          nodeSets[iset].second.push_back(inod);
      return iset+1;
    }

  nodeSets.emplace_back(setName,nodes);
  return nodeSets.size();
}


int ASMu2DLag::addToNodeSet (const std::string& setName, int inod, bool extId)
{
  int iset = this->getNodeSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  int nId = inod;
  if (extId) // Transform to internal node index
    inod = this->getNodeIndex(inod);

  if (inod > 0 && inod <= static_cast<int>(nnod))
    nodeSets[iset].second.push_back(inod);
  else
  {
    iset = -1;
    IFEM::cout <<"  ** Warning: Non-existing node "<< nId
               <<" in node set \""<< setName <<"\""<< std::endl;
  }

  return iset;
}


int ASMu2DLag::getElementSetIdx (const std::string& setName) const
{
  int iset = 1;
  for (const ASM::NodeSet& es : elemSets)
    if (es.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMu2DLag::getElementSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(elemSets.size()))
    return elemSets[iset-1].second;

  return this->ASMbase::getElementSet(iset);
}


bool ASMu2DLag::getElementSet (int iset, std::string& name) const
{
  if (iset < 1 || iset > static_cast<int>(elemSets.size()))
    return false;
  else if (elemSets[iset-1].second.empty())
    return false;

  name = elemSets[iset-1].first;
  return true;
}


/*!
  If \a iel is negative, the absolute value is taken as the external element ID.
  Otherwise, it is taken as the 1-based internal element index within the patch.
*/

bool ASMu2DLag::isInElementSet (int iset, int iel) const
{
  if (iset < 1 || iset > static_cast<int>(elemSets.size()))
    return false;

  if (iel < 0)
    iel = this->getElmIndex(-iel);

  return utl::findIndex(elemSets[iset-1].second,iel) >= 0;
}


int ASMu2DLag::parseElemSet (const std::string& setName, const char* cset)
{
  int iset = this->getElementSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = elemSets.size();
    elemSets.emplace_back(setName,IntVec());
  }

  IntVec& mySet = elemSets[iset].second;
  if (!strcmp(cset,"all"))
  {
    mySet.resize(nel);
    std::iota(mySet.begin(),mySet.end(),1);
    return 1+iset;
  }

  size_t ifirst = mySet.size();
  utl::parseIntegers(mySet,cset);

  // Transform to internal element indices
  for (size_t i = ifirst; i < mySet.size(); i++)
    if (mySet[i] > 0)
    {
      if (int iel = this->getElmIndex(mySet[i]); iel > 0)
        mySet[ifirst++] = iel;
      else
        IFEM::cout <<"  ** Warning: Non-existing element "<< mySet[i]
                   <<" in the set \""<< setName <<"\" (ignored)."<< std::endl;
    }

  if (ifirst < mySet.size())
    mySet.resize(ifirst);

  return 1+iset;
}


int ASMu2DLag::parseElemBox (const std::string& setName,
                             const std::string& unionSet, const char* data)
{
  int iuSet = unionSet.empty() ? 0 : this->getElementSetIdx(unionSet);

  Vec3 X0, X1;
  std::istringstream(data) >> X0 >> X1;

  // Lambda function for checking if an element is within the bounding box
  auto&& isInside=[this,iuSet,&X0,&X1](size_t iel)
  {
    if (iuSet > 0 && !this->isInElementSet(iuSet,1+iel))
      return false; // Filter out all elements not in the set unionSet

    Vec3 XC = std::accumulate(MNPC[iel].begin(), MNPC[iel].end(), Vec3(),
                              [&c = coord](const Vec3& X, int inod)
                              { return X + c[inod]; }) / MNPC[iel].size();
    for (size_t j = 0; j < nsd; j++)
      if (XC[j] < X0[j] || XC[j] > X1[j])
        return false;

    return true;
  };

  Matrix Xnod;
  IntVec elems;
  for (size_t iel = 0; iel < nel; iel++)
    if (isInside(iel))
      elems.push_back(1+iel);

  IFEM::cout <<"\tBounding Box: "<< X0 <<" - "<< X1
             <<": "<< elems.size() <<" elements"<< std::endl;

  if (elems.empty())
    return 0; // No elements are within the given box

  for (size_t iset = 0; iset < elemSets.size(); iset++)
    if (elemSets[iset].first == setName)
    {
      if (elemSets[iset].second.empty())
        elemSets[iset].second.swap(elems);
      else for (int iel : elems)
        if (utl::findIndex(elemSets[iset].second,iel) < 0)
          elemSets[iset].second.push_back(iel);
      return iset+1;
    }

  elemSets.emplace_back(setName,elems);
  return elemSets.size();
}


int ASMu2DLag::addToElemSet (const std::string& setName, int iel, bool extId)
{
  int iset = this->getElementSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = elemSets.size();
    elemSets.emplace_back(setName,IntVec());
  }

  int eId = iel;
  if (extId) // Transform to internal element index
    iel = this->getElmIndex(iel);

  if (iel > 0 && iel <= static_cast<int>(nel))
    elemSets[iset].second.push_back(iel);
  else
  {
    iset = -1;
    IFEM::cout <<"  ** Warning: Non-existing element "<< eId
               <<" in element set \""<< setName <<"\""<< std::endl;
  }

  return iset;
}


void ASMu2DLag::getBoundaryNodes (int lIndex, IntVec& nodes,
                                  int, int, int, bool local) const
{
  nodes = this->getNodeSet(lIndex);
  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


void ASMu2DLag::generateThreadGroups (const Integrand&, bool silence,
                                      bool separateGroup1noded)
{
#ifdef USE_OPENMP
  if (omp_get_max_threads() > 1 && threadGroups.stripDir != ThreadGroups::NONE)
    this->generateThreadGroupsMultiColored(silence, separateGroup1noded);
  else
#endif
    threadGroups.oneGroup(nel); // No threading, all elements in one group
}


void ASMu2DLag::generateThreadGroupsMultiColored (bool silence,
                                                  bool separateGroup1noded)
{
  threadGroups[1].clear();
  threadGroups[0].clear();

  // Status vector for the elements:
  // -1 is unusable for current color, 0 is available,
  // any other value is the assigned color
  IntVec status(nel,0);

  using IntSet = std::set<int>;
  std::vector<IntSet> nodeConn(nnod); // node-to-element connectivity
  std::vector<IntSet> nodeNode(nnod); // node-to-node connectivity via MPCs

  for (const MPC* mpc : mpcs)
    if (int slave = this->getNodeIndex(mpc->getSlave().node); slave > 0)
      for (size_t i = 0; i < mpc->getNoMaster(); i++)
        if (int mastr = this->getNodeIndex(mpc->getMaster(i).node); mastr > 0)
        {
          nodeNode[slave-1].insert(mastr-1);
          nodeNode[mastr-1].insert(slave-1);
        }

  size_t fixedElements = 0;
  for (size_t iel = 0; iel < nel; iel++)
    if (separateGroup1noded && MNPC[iel].size() == 1 &&
        nodeNode[MNPC[iel].front()].empty())
    {
      status[iel] = 1; // Separate group for all single-noded elements
      if (++fixedElements == 1)
        threadGroups[0].resize(1, { static_cast<int>(iel) });
      else
        threadGroups[0].front().push_back(iel);
    }
    else
      for (int node : MNPC[iel])
      {
        nodeConn[node].insert(iel);
        for (int master : nodeNode[node])
          nodeConn[master].insert(iel);
      }

  for (size_t nColors = fixedElements > 0; fixedElements < nel; ++nColors)
  {
    // Reset un-assigned element tags
    std::for_each(status.begin(), status.end(),
                  [](int& s) { if (s < 0) s = 0; });

    // Look for available elements
    IntVec thisColor;
    for (size_t i = 0; i < nel; ++i)
      if (status[i] == 0)
      {
        status[i] = nColors + 1;
        thisColor.push_back(i);
        ++fixedElements;

        for (int node : MNPC[i])
          for (int j : nodeConn[node])
            if (status[j] == 0) // if not assigned a color yet
              status[j] = -1;   // set as unavailable (with current color)
      }

    threadGroups[0].emplace_back(thisColor);
  }

  if (!silence)
    threadGroups.analyzeUnstruct(true);
}


bool ASMu2DLag::tesselate (ElementBlock& grid, const int*) const
{
  size_t i, nmnpc = 0, nelms = nel;
  for (i = 0; i < nel; i++)
    if (MNPC[i].size() > 1 && MLGE[i] > 0)
      nmnpc += MNPC[i].size();
    else // ignore 1-noded elements (point masses, etc.) and collapsed elements
      --nelms;

  grid.unStructResize(nelms,nnod,nmnpc);

  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  size_t j, k, e;
  for (i = k = e = 0; i < nel; i++)
    if (MNPC[i].size() > 1 && MLGE[i] > 0)
    {
      for (j = 0; j < MNPC[i].size(); j++)
        if (j > 1 && swapNode34 && MNPC[i].size() == 4)
          grid.setNode(k++,MNPC[i][5-j]);
        else
          grid.setNode(k++,MNPC[i][j]);
      grid.endOfElm(k);
      grid.setElmId(++e,MLGE[i]);
    }

  return true;
}


bool ASMu2DLag::integrate (Integrand& integrand,
                           GlobalIntegral& glInt,
                           const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches
  if (!myElms.empty() && myElms.front() == -1) return true;

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));
  ASMs2D::BasisFunctionCache& cache = *myCache.front();

  cache.setIntegrand(&integrand);
  cache.init(integrand.getIntegrandType() & Integrand::NO_DERIVATIVES ? 0 : 1);

  ThreadGroups oneGroup;
  if (glInt.threadSafe())
    oneGroup.oneStripe(nel, myElms);
  const IntMat& group = glInt.threadSafe() ? oneGroup[0] : threadGroups[0];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < group.size() && ok; t++)
#pragma omp parallel for schedule(static)
    for (int iel : group[t])
      if (ok)
        ok = this->integrateElm(integrand,glInt,iel,cache,time);

  if (ASM::cachePolicy == ASM::PRE_CACHE)
    cache.clear();

  return ok;
}


bool ASMu2DLag::evalSolPt (int iel, double xi, double eta, size_t nCmp,
                           const Vector& pchSol, RealArray& ptSol,
                           RealArray& N) const
{
  ptSol.clear();
  size_t jel = iel > 0 ? iel : -iel;
  if (jel < 1 || jel > MNPC.size())
    return false;
  else if (MNPC[--jel].size() > 3) // 4-noded 2D element
    return this->ASMs2DLag::evalSolPt(iel,xi,eta,nCmp,pchSol,ptSol,N);
  else if (MNPC[jel].size() < 3)
    return true; // no results for 0D and 1D elements

  // This element is a linear triangle
  RealArray N3(3,1.0/3.0);
  if (iel > 0) N3 = { xi, eta, 1.0-xi-eta };
  return this->ASMs2DLag::evalSolPt(-1-jel,xi,eta,nCmp,pchSol,ptSol,N3);
}


bool ASMu2DLag::writeXML (const char* fname) const
{
  std::ofstream os(fname);
  if (!os) return false;

  os <<"<patch>\n  <nodes>";
  for (const Vec3& X : coord)
    os <<"\n    "<< X;
  os <<"\n  </nodes>";
  for (size_t nen = 1; nen <= 4; nen++)
  {
    bool haveElms = false;
    for (const IntVec& mnpc : MNPC)
      if (mnpc.size() == nen)
      {
        if (!haveElms)
          os <<"\n  <elements nenod=\""<< nen <<"\">";
        os <<"\n   ";
        for (int inod : mnpc)
          os <<" "<< inod;
        haveElms = true;
      }
    if (haveElms)
      os <<"\n  </elements>";
  }
  os <<"\n</patch>\n";

  return true;
}


bool ASMu2DLag::findPoints (Vec3Vec& points, ParamsVec& locs) const
{
  // Calculate the bounding box for all (shell) elements
  size_t nonshell = 0;
  std::vector<Vec3Pair> bbox(nel);
  for (size_t iel = 0; iel < nel; iel++)
    if (MNPC[iel].size() == 3 || MNPC[iel].size() == 4)
      this->getBoundingBox(MNPC[iel],bbox[iel]);
    else
    {
      nonshell++; // make sure any other elements are sorted last
      bbox[iel].first = bbox[iel].second = Vec3(1.0e99,1.0e99,1.0e99);
    }

  std::array<std::vector<size_t>,3> eIdx;
  std::vector<size_t>::iterator it;

  // Index sort the elements in each spatial direction w.r.t. the max coordinate
  for (int i = 0; i < 3; i++)
  {
    eIdx[i].resize(nel);
    std::iota(eIdx[i].begin(), eIdx[i].end(), 0);
    std::sort(eIdx[i].begin(), eIdx[i].end(), [i,&bbox](size_t a, size_t b)
              { return bbox[a].second[i] < bbox[b].second[i]; });
    if (nonshell > 0)
      eIdx[i].resize(nel-nonshell); // Remove the non-shell elements
#if SP_DEBUG > 2
    std::cout <<"\nElement order in "<< char('X'+i) <<"-direction:";
    for (size_t iel = 0; iel < nel-nonshell; iel++)
      std::cout << (iel%10 ? ' ' : '\n') << MLGE[eIdx[i][iel]]
                <<" ("<< bbox[eIdx[i][iel]].second[i] <<")";
    std::cout << std::endl;
#endif
  }

  Vec3Pair BBox; // Bounding box for the whole patch
  double maxDim = 0.0;
  for (int i = 0; i < 3; i++)
  {
    BBox.first[i]  = bbox[eIdx[i].front()].first[i];
    BBox.second[i] = bbox[eIdx[i].back()].second[i];
    if (double dim = BBox.second[i] - BBox.first[i]; dim > maxDim)
      maxDim = dim;
  }
#if SP_DEBUG > 1
  std::cout <<"\nDomain bounding box: [ "
            << BBox.first <<" , "<< BBox.second
            <<" ] maxDim = "<< maxDim << std::endl;
#endif
  const double tol = maxDim*1.0e-8;
  BBox.second += Vec3(RealArray(3,-tol));

  size_t nFound = 0;
  locs.clear();
  locs.reserve(points.size());
  for (Vec3& X : points)
  {
    // Do a binary search in each coordinate direction
    // to find the range of elements to check in more detail
    std::set<size_t> nearElms;
    for (int i = 0; i < 3; i++)
      if (BBox.second[i] > BBox.first[i]) // skip the flat dimension
      {
        it = std::upper_bound(eIdx[i].begin(), eIdx[i].end(), X[i],
                              [i,&bbox](double x, size_t a)
                              { return x < bbox[a].second[i]; });
        if (nearElms.empty())
          nearElms.insert(it,eIdx[i].end());
        else
        {
          std::set<size_t> tmp1(it,eIdx[i].end()), tmp2;
          nearElms.swap(tmp2);
          std::set_intersection(tmp1.begin(), tmp1.end(),
                                tmp2.begin(), tmp2.end(),
                                std::inserter(nearElms,nearElms.begin()));
        }
      }

    // Lambda function doing the Y < X comparison for two Vec3 objects
    auto&& below = [&X,tol](const Vec3& Y)
    {
      for (int i = 0; i < 3; i++)
        if (Y[i] > X[i]+tol)
          return false;
      return true;
    };

#if SP_DEBUG > 1
    std::cout <<"\nSearching spatial point "<< X <<" (max "<< nearElms.size()
              <<" elements to check)."<< std::endl;
#endif
    bool found = false;
    for (size_t iel : nearElms)
      if (below(bbox[iel].first))
      {
#if SP_DEBUG > 1
        std::cout <<"Checking element "<< MLGE[iel] <<" in [ "<< bbox[iel].first
                  <<" - "<< bbox[iel].second <<" ]"<< std::endl;
#endif
        // Perform parametric inversion to see if this point
        // really is within current element
        if (PointParams elem(1+iel); this->paramInvert(X,elem))
        {
          locs.emplace_back(elem);
          found = true;
          break;
        }
      }

    if (found)
      nFound++;
    else
    {
      locs.emplace_back(PointParams());
      IFEM::cout <<"  ** Warning: No element matches the point "<< X
                 << std::endl;
    }
  }

  return nFound == points.size();
}


bool ASMu2DLag::paramInvert (Vec3& X, PointParams& elm) const
{
  Matrix Xn;
  if (!this->getElementCoordinates(Xn,elm.iel))
    return false;

  const int nelnod = Xn.cols();
  if (nelnod == 4) // swap element node 3 and 4
    for (size_t r = 1; r <= Xn.rows(); r++)
      std::swap(Xn(r,3),Xn(r,4));

  const double eps = 1.0e-15;
  const double tol = 1.0e-3;
#if SP_DEBUG > 2
  int ipsw = SP_DEBUG;
#else
  int ipsw = 0;
#endif
  // Open a temporary file for the Fortran output
  std::string tmpfile = "/tmp/" + std::string(getenv("USER")) + ".ftn";
  const int iwr = openftnfile_(tmpfile.c_str(),tmpfile.size());
  if (iwr <= 0) ipsw = -1; // suppress all output
  int ierr = 0;
  if (nelnod == 3)
    t3_inv_(nsd,tol,eps,X.ptr(),Xn.ptr(),elm.u,ipsw,iwr,ierr);
  else if (nelnod == 4)
    q4_inv_(nsd,tol,eps,X.ptr(),Xn.ptr(),elm.u,ipsw,iwr,ierr);
  else
    ierr = -99;
  if (iwr > 0 && ipsw >= 0)
  {
    // Copy fortran output to IFEM::cout
    closeftnfile_(iwr);
    std::string cline;
    std::ifstream is(tmpfile);
    while (std::getline(is,cline))
      IFEM::cout << cline << std::endl;
  }
  if (ierr != 0 && ierr != 2)
    return false;

  // Calculate the distance to the projected point
  Vector N(nelnod);
  shape2_(0,nelnod,elm.u,N.ptr(),ipsw,iwr,ierr);
  Vec3 X0(Xn * N);
  elm.dist = (X0-X).length();

  X = X0; // update the spatial coordinates
  return true;
}
