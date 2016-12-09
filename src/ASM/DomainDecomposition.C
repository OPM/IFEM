// $Id$
//==============================================================================
//!
//! \file DomainDecomposition.C
//!
//! \date Feb 23 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Domain decomposition partitioning for structured models.
//!
//==============================================================================

#include "DomainDecomposition.h"
#include "ASMstruct.h"
#include "ASM2D.h"
#include "ASM3D.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAMpatch.h"
#include "SIMbase.h"
#include "Utilities.h"
#include "Vec3.h"
#include "IFEM.h"
#include <numeric>


//! \brief Iterator for matching nodes on edges/faces with a given orientation and index.
class NodeIterator {
public:
  //! \brief Constructor.
  //! \param pch Slave patch
  //! \param orient Orientation of boundary on slave
  //! \param lIdx Index of boundary on slave
  //! \param basis Basis to use for boundary on slave
  //! \param dim Dimension to iterate over
  NodeIterator(const ASMstruct* pch, int orient, int lIdx, int basis, int dim = 2)
  {
    if (dim == 0) {
      nodes.resize(1);
      return;
    }

    int n1, n2, n3;
    pch->getSize(n1,n2,n3,basis);
    int nsd = pch->getNoSpaceDim();
    if (nsd == 3) {
      if (dim == 1) {
        if (lIdx <= 4)
          nodes.resize(n1);
        else if (lIdx >= 5 && lIdx <= 8)
          nodes.resize(n2);
        else
          nodes.resize(n3);

        if (orient == 0)
          std::iota(nodes.begin(), nodes.end(), 0);
        else
          std::iota(nodes.rbegin(), nodes.rend(), 0);

        return;
      }

      int dim1, dim2;
      if (lIdx == 1 || lIdx == 2)
        dim1 = n2, dim2 = n3;
      else if (lIdx == 3 || lIdx == 4)
        dim1 = n1, dim2 = n3;
      else //if (lIdx == 5 || lIdx == 6)
        dim1 = n1, dim2 = n2;

      nodes.resize(dim1*dim2);
      typedef std::function<std::pair<int,int>(int dim1, int dim2, int n)> NodeOrder;
#define PAIR(x, y) [](int dim1, int dim2, int n) { return std::make_pair(x,y); }
      const std::vector<NodeOrder> orders  = {PAIR(n*dim1,              1),  // 0 1 2 3 4 5 6 7 8
                                              PAIR((dim2-n-1)*dim1,     1),  // 6 7 8 3 4 5 0 1 2
                                              PAIR((n+1)*dim1-1,       -1),  // 2 1 0 5 4 3 8 7 6
                                              PAIR(dim1*dim2-n*dim1-1, -1),  // 8 7 6 5 4 3 2 1 0
                                              PAIR(n,                dim1),  // 0 3 6 1 4 7 2 5 8
                                              PAIR(dim2-n-1,         dim1),  // 2 5 8 1 4 7 0 3 6
                                              PAIR((dim2-1)*dim1+n, -dim1),  // 6 3 0 7 4 1 8 5 2
                                              PAIR(dim2*dim1-n-1,   -dim1)}; // 8 5 2 7 4 1 6 3 0
#undef PAIR

      auto it = nodes.begin();
      for (int n = 0; n < dim2; ++n) {
        auto order = orders[orient](dim1, dim2, n);
        for (int i = 0; i < dim1; ++i, ++it, order.first += order.second)
          *it = order.first;
      }
    } else {
      if (lIdx == 1 || lIdx == 2)
        nodes.resize(n2);
      else
        nodes.resize(n1);

      if (orient == 0)
        std::iota(nodes.begin(), nodes.end(), 0);
      else
        std::iota(nodes.rbegin(), nodes.rend(), 0);
    }
  }

  //! \brief Obtain start of node numbers.
  std::vector<int>::const_iterator begin() { return nodes.begin(); }
  //! \brief Obtain end of node numbers.
  std::vector<int>::const_iterator end() { return nodes.end(); }
  //! \brief Obtain number of nodes
  size_t size() const { return nodes.size(); }
protected:
  //! \brief Node number on boundary.
  std::vector<int> nodes;
};


std::vector<std::set<int>> DomainDecomposition::getSubdomains(int nx, int ny, int nz,
                                                              int overlap, size_t block) const
{
  std::vector<IntSet> result(nx*(ny?ny:1)*(nz?nz:1)*getSAM()->getNoPatches());
  size_t d = 0;
  for (const auto& it : *getSAM()) {
    const ASMstruct* pch = dynamic_cast<const ASMstruct*>(it);
    if (!pch)
      break;
    int n1, n2, n3;
    pch->getNoStructElms(n1,n2,n3);
    auto subdomains = calcSubdomains(n1, n2, n3, nx, ny, nz, overlap);
    for (size_t g = 0; g < subdomains.size(); ++g, ++d) {
      for (const int& iEl : subdomains[g]) {
        if (getNoBlocks() == 0) {
          IntVec eqns;
          int el = it->getElmID(iEl+1);
          if (el == 0)
            continue;

          getSAM()->getElmEqns(eqns, el);
          for (auto& it : eqns) {
            if (it > 0)
              result[d].insert(it);
          }
        } else {
          IntVec nodes;
          getSAM()->getElmNodes(nodes, it->getElmID(iEl+1));
          int basis = blocks[block+1].basis;
          for (auto& node : nodes) {
            char type = getSAM()->getNodeType(node);
            if (type == (basis == 1 ? 'D' : 'P'+basis-2) || (type == ' ' && basis < 2)) {
              IntVec eqns;
              getSAM()->getNodeEqns(eqns, node);
              for (auto& it : eqns) {
                if (it > 0) {
                  auto geq = blocks[block+1].G2LEQ.find(it);
                  result[d].insert(geq->second);
                }
              }
            }
          }
        }
      }
    }
  }

  return result;
}


std::vector<std::vector<int>> DomainDecomposition::calcSubdomains(size_t nel1, size_t nel2, size_t nel3,
                                                                  size_t g1, size_t g2, size_t g3, size_t overlap)
{
  if (nel3 > 0)
    return calcSubdomains3D(nel1, nel2, nel3, g1, g2, g3, overlap);
  else if (nel2 > 0)
    return calcSubdomains2D(nel1, nel2, g1, g2, overlap);
  else
    return calcSubdomains1D(nel1, g1, overlap);
}


std::vector<IntVec> DomainDecomposition::calcSubdomains1D(size_t nel1, size_t g1, size_t overlap)
{
  std::vector<IntVec> subdomains;

  if (g1 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    size_t nel1_sub = floor(double(nel1)/g1 + 1);
    subdomains.resize(g1);
    size_t ofs1 = 0;
    for (size_t gu = 0; gu < g1; ++gu) {
      if (gu == g1 - 1)
        nel1_sub = nel1 - ofs1;
      subdomains[gu].reserve(nel1_sub);
      for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
        subdomains[gu].push_back(ofs1 + i);
      ofs1 += nel1_sub - overlap;
    }
  }

  return subdomains;
}


std::vector<IntVec> DomainDecomposition::calcSubdomains2D(size_t nel1, size_t nel2,
                                                          size_t g1, size_t g2, size_t overlap)
{
  std::vector<IntVec> subdomains;

  if (g1*g2 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1*nel2);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    size_t nel2_sub = floor(double(nel2)/g2 + 1);
    subdomains.resize(g1*g2);
    size_t g = 0;
    size_t ofs2 = 0;
    for (size_t gv = 0; gv < g2; ++gv) {
      if (gv == g2 - 1)
        nel2_sub = nel2 - ofs2;
      size_t nel1_sub = floor(double(nel1)/g1 + 1);
      size_t ofs1 = 0;
      for (size_t gu = 0; gu < g1; ++gu, ++g) {
        if (gu == g1 - 1)
          nel1_sub = nel1-ofs1;

        subdomains[g].reserve(nel1_sub*nel2_sub);
        for (size_t j = 0; j < nel2_sub && ofs2 + j < nel2; ++j)
          for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
            subdomains[g].push_back((ofs2+j)*nel1 + ofs1 + i);
        ofs1 += nel1_sub - overlap;
      }
      ofs2 += nel2_sub - overlap;
    }
  }

  return subdomains;
}


std::vector<IntVec> DomainDecomposition::calcSubdomains3D(size_t nel1, size_t nel2, size_t nel3,
                                                          size_t g1, size_t g2, size_t g3, size_t overlap)
{
  std::vector<IntVec> subdomains;

  if (g1*g2*g3 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1*nel2*nel3);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    subdomains.resize(g1*g2*g3);
    size_t g = 0;
    size_t ofs3 = 0;
    size_t nel3_sub = floor(double(nel3)/g3+1);
    for (size_t gw = 0; gw < g3; ++gw) {
      if (gw == g3 - 1)
        nel3_sub = nel3-ofs3;
      size_t ofs2 = 0;
      size_t nel2_sub = floor(double(nel2)/g2+1);
      for (size_t gv = 0; gv < g2; ++gv) {
        if (gv == g2 - 1)
          nel2_sub = nel2-ofs2;
        size_t ofs1 = 0;
        size_t nel1_sub = floor(double(nel1)/g1+1);
        for (size_t gu = 0; gu < g1; ++gu, ++g) {
          if (gu == g1 - 1)
            nel1_sub = nel1-ofs1;
          subdomains[g].reserve(nel1_sub*nel2_sub*nel3_sub);
          for (size_t k = 0; k < nel3_sub && ofs3 + k < nel3; ++k)
            for (size_t j = 0; j < nel2_sub && ofs2 + j < nel2; ++j)
              for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
                subdomains[g].push_back((ofs3 + k)*nel1*nel2 + (ofs2+j)*nel1 + ofs1 + i);
          ofs1 += nel1_sub - overlap;
        }
        ofs2 += nel2_sub - overlap;
      }
      ofs3 += nel3_sub - overlap;
    }
  }

  return subdomains;
}


void DomainDecomposition::setupNodeNumbers(int basis, IntVec& lNodes,
                                           std::set<int>& cbasis,
                                           const ASMbase* pch,
                                           int dim, int lidx, int thick)
{
  if (basis != 0) // specified base
    cbasis = utl::getDigits(basis);
  else if (dim == 0 || (dim == 1 && pch->getNoSpaceDim() == 3))
    // need to expand to all bases for corners and edges
    for (size_t b = 1; b <= pch->getNoBasis(); ++b)
      cbasis.insert(b);
  else {// directly add nodes, cbasis remains empty
    pch->getBoundaryNodes(lidx, lNodes, 0, thick, false);
  }

  const ASM2D* pch2D = dynamic_cast<const ASM2D*>(pch);
  const ASM3D* pch3D = dynamic_cast<const ASM3D*>(pch);
  for (const int& it2 : cbasis)
    if (dim == 0) {
      int node = 0;
      if (pch2D)
        switch (lidx) {
        case 1: node = pch2D->getCorner(-1,-1, it2); break;
        case 2: node = pch2D->getCorner( 1,-1, it2); break;
        case 3: node = pch2D->getCorner(-1, 1, it2); break;
        case 4: node = pch2D->getCorner( 1, 1, it2); break;
        }
      else if (pch3D)
        switch (lidx) {
        case 1: node = pch3D->getCorner(-1,-1,-1, it2); break;
        case 2: node = pch3D->getCorner( 1,-1,-1, it2); break;
        case 3: node = pch3D->getCorner(-1, 1,-1, it2); break;
        case 4: node = pch3D->getCorner( 1, 1,-1, it2); break;
        case 5: node = pch3D->getCorner(-1,-1, 1, it2); break;
        case 6: node = pch3D->getCorner( 1,-1, 1, it2); break;
        case 7: node = pch3D->getCorner(-1, 1, 1, it2); break;
        case 8: node = pch3D->getCorner( 1, 1, 1, it2); break;
        }
      lNodes.push_back(pch->getNodeID(node));
    } else if (dim == 1 && pch3D) {
      std::vector<int> eNodes = pch3D->getEdge(lidx, false, it2);
      for (const int& it : eNodes)
        lNodes.push_back(pch->getNodeID(it));
    } else {
      pch->getBoundaryNodes(lidx, lNodes, it2, thick, false);
    }
}


bool DomainDecomposition::calcGlobalNodeNumbers(const ProcessAdm& adm,
                                                const SIMbase& sim)
{
  minDof = minNode = 1;
  maxDof = sim.getSAM()->getNoDOFs();
  maxNode = sim.getSAM()->getNoNodes();

#ifdef HAVE_MPI
  if (adm.getNoProcs() == 1)
    return true;

  IFEM::cout << "\tEstablishing global node numbers" << std::endl;

  std::vector<int> locLMs, glbLMs;
  for (size_t n = 1; n <= sim.getPatch(1)->getNoNodes(); ++n) {
    if (sim.getPatch(1)->getLMType(n) == 'G')
      locLMs.push_back(n);
  }

  minDof = 1;
  maxDof = sim.getSAM()->getNoDOFs();
  minNode = 1;
  maxNode = sim.getSAM()->getNoNodes();
  if (adm.getProcId() > 0) {
    adm.receive(minNode, adm.getProcId()-1);
    adm.receive(minDof, adm.getProcId()-1);

    int nLMs;
    adm.receive(nLMs, adm.getProcId()-1);
    if (locLMs.size() != (size_t)nLMs) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers():"
                <<" Non-matching number of multipliers "
                << nLMs << ", expected " << locLMs.size() << std::endl;
      return false;
    }
    glbLMs.resize(nLMs);
    adm.receive(glbLMs, adm.getProcId()-1);

    maxDof  = minDof++;
    maxNode = minNode++;
  }

  MLGN.resize(sim.getSAM()->getNoNodes());
  std::iota(MLGN.begin(), MLGN.end(), minNode);

  std::map<int,int> old2new;
  for (const auto& it : ghostConnections) {
    int sidx = sim.getLocalPatchIndex(it.slave);
    if (sidx < 1)
      continue;

    std::set<int> cbasis;
    IntVec lNodes;

    setupNodeNumbers(it.basis, lNodes, cbasis, sim.getPatch(sidx),
                     it.dim, it.sidx, it.thick);

    int nRecv;
    adm.receive(nRecv, getPatchOwner(it.master));
    if (nRecv =! lNodes.size()) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): "
                <<" Topology error, boundary size "
                << nRecv << ", expected " << lNodes.size() << std::endl;
      return false;
    }
    IntVec glbNodes(lNodes.size());
    adm.receive(glbNodes, getPatchOwner(it.master));
    size_t ofs = 0;
    for (size_t b = 1; b <= sim.getPatch(sidx)->getNoBasis(); ++b) {
      if (cbasis.empty() || cbasis.find(b) != cbasis.end()) {
        NodeIterator iter(dynamic_cast<const ASMstruct*>(sim.getPatch(sidx)),
                          it.orient, it.sidx, b, it.dim);
        auto it_n = iter.begin();
        for (size_t i = 0; i < iter.size(); ++i, ++it_n) {
          for (int t = 0; t < it.thick; ++t) {
            int node = MLGN[lNodes[i*it.thick+t+ofs]-1];
            old2new[node] = glbNodes[*it_n*it.thick+t + ofs];
          }
        }
        ofs += iter.size()*it.thick;
      }
    }
  }
  // add multiplier remappings
  for (size_t i = 0; i < locLMs.size() && adm.getProcId() > 0; ++i)
    old2new[MLGN[sim.getPatch(1)->getNodeID(locLMs[i])-1]] = glbLMs[i];

  // remap ghost nodes
  for (auto& it : MLGN)
    utl::renumber(it, old2new, false);

  // remap rest of our nodes
  for (int i = 0; i < sim.getSAM()->getNoNodes() && adm.getProcId() != 0; ++i)
    if (old2new.find(i + minNode) == old2new.end()) {
      std::map<int,int> old2new2;
      old2new2[i + minNode] = ++maxNode;
      auto dof = sim.getSAM()->getNodeDOFs(i);
      maxDof += dof.second-dof.first+1;
      for (auto& it : MLGN)
        utl::renumber(it, old2new2, false);
    }

  if (adm.getProcId() < adm.getNoProcs()-1) {
    adm.send(maxNode, adm.getProcId()+1);
    adm.send(maxDof, adm.getProcId()+1);
    for (auto& it : locLMs)
      it = MLGN[sim.getPatch(1)->getNodeID(it)-1];

    adm.send((int)locLMs.size(), adm.getProcId()+1);
    adm.send(locLMs, adm.getProcId()+1);
  }

  for (const auto& it : ghostConnections) {
    int midx = sim.getLocalPatchIndex(it.master);
    if (midx < 1)
      continue;

    std::set<int> cbasis;
    std::vector<int> glbNodes;
    setupNodeNumbers(it.basis, glbNodes, cbasis, sim.getPatch(midx),
                     it.dim, it.midx, it.thick);

    for (size_t i = 0; i < glbNodes.size(); ++i)
      glbNodes[i] = MLGN[glbNodes[i]-1];

    adm.send(int(glbNodes.size()), getPatchOwner(it.slave));
    adm.send(glbNodes, getPatchOwner(it.slave));
  }
#endif

  return true;
}


std::vector<int> DomainDecomposition::setupEquationNumbers(const SIMbase& sim,
                                                           int pidx, int lidx,
                                                           const std::set<int>& cbasis,
                                                           int dim, int thick)
{
  std::vector<IntVec> lNodes(sim.getPatch(pidx)->getNoBasis());
  std::vector<int> result;

  for (size_t block = 0; block < blocks.size(); ++block) {
    std::set<int> bases;
    if (block == 0) {
      for (size_t i = 1; i <= sim.getPatch(pidx)->getNoBasis(); ++i)
        bases.insert(i);
    } else
      bases = utl::getDigits(sim.getSolParams()->getBlock(block-1).basis);

    for (const int& basis : bases) {
      if (!cbasis.empty() && cbasis.find(basis) == cbasis.end())
        continue;

      if (lNodes[basis-1].empty()) {
        IntSet dummy;
        setupNodeNumbers(basis, lNodes[basis-1], dummy, sim.getPatch(pidx),
                         dim, lidx, thick);
      }

      std::set<int> components;
      if (block == 0 || sim.getSolParams()->getBlock(block-1).comps == 0) {
        for (size_t i = 1; i <= sim.getPatch(pidx)->getNoFields(basis); ++i)
          components.insert(i);
      } else
        components = utl::getDigits(sim.getSolParams()->getBlock(block-1).comps);

      for (size_t i = 0; i < lNodes[basis-1].size(); ++i) {
        int node = lNodes[basis-1][i];
        for (const int& dof : components) {
          int eq = sim.getSAM()->getEquation(node, dof);
          if (eq > 0) {
            if (block == 0)
              result.push_back(blocks[block].MLGEQ[eq-1]);
            else
              result.push_back(blocks[block].MLGEQ[blocks[block].G2LEQ[eq]-1]);
          } else
            result.push_back(0);
        }
      }
    }
  }

  return result;
}


bool DomainDecomposition::calcGlobalEqNumbers(const ProcessAdm& adm,
                                              const SIMbase& sim)
{
  // defaults in serial and non-parallel runs with MPI enabled builds.
  blocks[0].minEq = 1;
  blocks[0].maxEq = sim.getSAM()->getNoEquations();
  for (size_t i = 1; i < blocks.size(); ++i) {
    blocks[i].minEq = 1;
    blocks[i].maxEq = blocks[i].localEqs.size();
    blocks[i].nGlbEqs = blocks[i].localEqs.size();
  }

#ifdef HAVE_MPI
  if (adm.getNoProcs() == 1)
    return true;

  IFEM::cout << "\tEstablishing global equation numbers" << std::endl;

  std::vector<int> locLMs;
  for (size_t n = 1; n <= sim.getPatch(1)->getNoNodes(); ++n) {
    if (sim.getPatch(1)->getLMType(n) == 'G')
      locLMs.push_back(n);
  }
  std::vector<int> glbLMs;
  std::vector<int> blkLMs;
  std::vector<int> nEqs(blocks.size());
  if (adm.getProcId() > 0) {
    adm.receive(nEqs, adm.getProcId()-1);
    int nLMs;
    adm.receive(nLMs, adm.getProcId()-1);

    if (locLMs.size() != (size_t)nLMs) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalEqNumbers():"
                <<" Non-matching number of multipliers "
                << nLMs << ", expected " << locLMs.size() << std::endl;
      return false;
    }

    glbLMs.resize(nLMs);
    adm.receive(glbLMs, adm.getProcId()-1);
    // HACK: multipliers always in second block
    if (blocks.size() > 1) {
      blkLMs.resize(nLMs);
      adm.receive(blkLMs, adm.getProcId()-1);
    }
  }

  for (size_t block = 0; block < blocks.size(); ++block) {
    size_t size = block == 0 ? sim.getSAM()->getNoEquations() :
                               blocks[block].localEqs.size();
    blocks[block].MLGEQ.resize(size);
    std::iota(blocks[block].MLGEQ.begin(),
              blocks[block].MLGEQ.end(), nEqs[block]+1);
    if (adm.getProcId() > 0) {
      blocks[block].minEq = nEqs[block]+1;
      blocks[block].maxEq = nEqs[block];
    } else {
      blocks[block].minEq = 1;
      blocks[block].maxEq = size;
    }
  }

  std::vector<std::map<int,int>> old2new(blocks.size());
  std::vector<std::vector<bool>> o2nu(blocks.size()); // flat vector cache
  for (size_t block = 0; block < o2nu.size(); ++block)
    if (block == 0)
      o2nu[block].resize(sim.getSAM()->getNoEquations(), false);
    else
      o2nu[block].resize(adm.dd.getBlockEqs(block-1).size(), false);


  size_t n = 0;
  auto blockIt = blkLMs.begin();
  for (const auto& it : glbLMs) {
    int seq = sim.getSAM()->getEquation(sim.getPatch(1)->getNodeID(locLMs[n++]), 1);
    int leq = blocks[0].MLGEQ[seq-1];
    old2new[0][leq] = it;
    o2nu[0][leq-nEqs[0]-1] = true;
    if (blocks.size() > 1 && !blkLMs.empty()) {
      int leq = blocks[2].MLGEQ[blocks[2].G2LEQ[seq]-1];
      old2new[2][leq] = *blockIt;
      o2nu[2][leq-nEqs[2]-1] = true;
      ++blockIt;
    }
  }

  for (const auto& it : ghostConnections) {
    int sidx = sim.getLocalPatchIndex(it.slave);
    if (sidx < 1)
      continue;

    std::set<int> cbasis;
    if (it.basis != 0)
      cbasis = utl::getDigits(it.basis);

    IntVec locEqs = setupEquationNumbers(sim, sidx, it.sidx,
                                         cbasis, it.dim, it.thick);

    int nRecv;
    adm.receive(nRecv, getPatchOwner(it.master));
    if (nRecv != (int)locEqs.size()) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalEqNumbers():"
                <<" Topology error, number of equations "
                << nRecv << ", expected " << locEqs.size() << std::endl;
      return false;
    }

    IntVec glbEqs(locEqs.size());
    adm.receive(glbEqs, getPatchOwner(it.master));

    size_t ofs = 0;
    for (size_t block = 0; block < blocks.size(); ++block) {
      std::set<int> bases;
      if (block == 0) {
        for (size_t i = 1; i <= sim.getPatch(sidx)->getNoBasis(); ++i)
          bases.insert(i);
      } else
        bases = utl::getDigits(sim.getSolParams()->getBlock(block-1).basis);

      for (const int& basis : bases) {
        if (!cbasis.empty() && cbasis.find(basis) == cbasis.end())
          continue;

        std::set<int> components;
        if (block == 0 || sim.getSolParams()->getBlock(block-1).comps == 0) {
          for (size_t i = 1; i <= sim.getPatch(sidx)->getNoFields(basis); ++i)
            components.insert(i);
        } else
          components = utl::getDigits(sim.getSolParams()->getBlock(block-1).comps);

        NodeIterator iter(dynamic_cast<const ASMstruct*>(sim.getPatch(sidx)),
                          it.orient, it.sidx, basis, it.dim);
        auto it_n = iter.begin();
        int nodeDofs = components.size();
        for (size_t i = 0; i < iter.size(); ++i, ++it_n) {
          for (int t = 0; t < it.thick; ++t) {
            for (size_t comp = 0; comp < components.size(); ++comp) {
              int leq = locEqs[ofs+(i*it.thick+t)*nodeDofs+comp];
              int geq = glbEqs[ofs+(*it_n*it.thick+t)*nodeDofs+comp];
              if ((leq < 1 && geq > 0) || (leq > 0 && geq < 1)) {
                std::cerr <<"\n *** DomainDecomposition::calcGlobalEqNumbers(): "
                          << "Topology error on process " << adm.getProcId()
                          << ", dof constraint mismatch "
                          << "local: " << leq << " global " << geq << std::endl;
                return false;
              }

              if (leq < 1)
                continue;

              old2new[block][leq] = geq;
              o2nu[block][leq-nEqs[block]-1] = true;
            }
          }
        }
        ofs += iter.size()*nodeDofs*it.thick;
      }
    }
  }

  for (size_t block = 0; block < blocks.size() && adm.getProcId() > 0; ++block) {
    // remap ghost equations
    for (auto& it : blocks[block].MLGEQ)
      utl::renumber(it, old2new[block], false);

    // remap the rest of our equations
    size_t size = block == 0 ? sim.getSAM()->getNoEquations() :
                               blocks[block].localEqs.size();
    std::map<int,int> old2new2;
    for (size_t i = 1; i <= size; ++i) {
      if (!o2nu[block][i-1])
        old2new2[i + nEqs[block]] = ++blocks[block].maxEq;
    }
    for (auto& it : blocks[block].MLGEQ)
      utl::renumber(it, old2new2, false);
  }

  if (adm.getProcId() < adm.getNoProcs()-1) {
    std::vector<int> maxEqs;
    for (auto& it : blocks)
      maxEqs.push_back(it.maxEq);

    adm.send(maxEqs, adm.getProcId()+1);

    std::vector<int> LM;
    size_t nMult = 0;
    for (size_t n = 0; n < locLMs.size(); ++n) {
      // TODO: > 1 dof for multipliers
      int eq = sim.getSAM()->getEquation(sim.getPatch(1)->getNodeID(locLMs[n]), 1);
      if (eq > 0)
        LM.push_back(blocks[0].MLGEQ[eq-1]), ++nMult;
    }

    adm.send((int)LM.size(), adm.getProcId()+1);
    adm.send(LM, adm.getProcId()+1);

    // HACK: multipliers always in second block
    LM.clear();
    if (blocks.size() > 1) {
      for (size_t i = 0; i < nMult; ++i)
        LM.push_back(blocks[2].MLGEQ[blocks[2].MLGEQ.size()-nMult+i]);
      adm.send(LM, adm.getProcId()+1);
    }
  }

  for (const auto& it : ghostConnections) {
    int midx = sim.getLocalPatchIndex(it.master);
    if (midx < 1)
      continue;

    std::set<int> cbasis;
    if (it.basis != 0)
      cbasis = utl::getDigits(it.basis);

    IntVec glbEqs = setupEquationNumbers(sim, midx, it.midx,
                                         cbasis, it.dim, it.thick);
    adm.send(int(glbEqs.size()), getPatchOwner(it.slave));
    adm.send(glbEqs, getPatchOwner(it.slave));
  }
#endif

  return true;
}


int DomainDecomposition::getGlobalEq(int lEq, size_t idx) const
{
  if (lEq < 1 || idx >= blocks.size())
    return 0;

  if (blocks[0].MLGEQ.empty()) { // serial
    if (lEq > blocks[idx].maxEq)
      return 0;

    return lEq;
  }

  if (!blocks[idx].MLGEQ.empty() && lEq > (int)blocks[idx].MLGEQ.size())
    return 0;

  return blocks[idx].MLGEQ[lEq-1];
}


bool DomainDecomposition::sanityCheckCorners(const SIMbase& sim)
{
#ifdef HAVE_MPI
  const ProcessAdm& adm = sim.getProcessAdm();
  if (!adm.isParallel())
    return true;

  std::vector<int> sizes(adm.getNoProcs());
  for (int i = 1; i <= sim.getNoPatches(); ++i)
    sizes[adm.dd.getPatchOwner(i)] += 5*sim.getPatch(1)->getNoBasis() *
                                      pow(2,sim.getNoSpaceDim());

  std::vector<int> displ(adm.getNoProcs());
  for (int i = 1; i < adm.getNoProcs(); ++i)
    displ[i] = displ[i-1] + sizes[i-1];

  std::vector<double> loc_data;
  for (int i = 1; i <= sim.getNoPatches(); ++i) {
    int pIdx;
    if ((pIdx=sim.getLocalPatchIndex(i)) > 0) {
      const ASM2D* pch2D = nullptr;
      if (sim.getNoSpaceDim() == 2)
        pch2D = dynamic_cast<const ASM2D*>(sim.getPatch(pIdx));
      const ASM3D* pch3D = nullptr;
      if (sim.getNoSpaceDim() == 3)
        pch3D = dynamic_cast<const ASM3D*>(sim.getPatch(pIdx));
      if (pch2D || pch3D) {
        for (size_t c = 0; c < pow(2,sim.getPatch(pIdx)->getNoSpaceDim()); ++c) {
          for (size_t b = 1; b <= sim.getPatch(pIdx)->getNoBasis(); ++b) {
            int node = 0;
            if (pch2D)
              node = pch2D->getCorner((c == 1 || c == 3) ? 1 : -1,
                                                 c >= 2  ? 1 : -1, b);
            else if (pch3D)
              node = pch3D->getCorner((c == 1 || c == 3 || c == 5 || c == 7) ? 1 : -1,
                                      (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : -1,
                                                                     c >= 4  ? 1 : -1, b);

            Vec3 pos = sim.getPatch(pIdx)->getCoord(node);
            loc_data.push_back(pos[0]);
            loc_data.push_back(pos[1]);
            loc_data.push_back(pos[2]);
            loc_data.push_back(b);
            loc_data.push_back(adm.dd.getMLGN()[sim.getPatch(pIdx)->getNodeID(node)-1]);
          }
        }
      }
    }
  }

  std::vector<double> glob_data(displ.back()+sizes.back());
  MPI_Allgatherv(loc_data.data(), loc_data.size(), MPI_DOUBLE,
                 glob_data.data(), sizes.data(), displ.data(), MPI_DOUBLE,
                 *adm.getCommunicator());

  // unpack data
  std::vector<std::array<double,5>> corners;
  for (size_t i = 0; i < glob_data.size(); i += 5)
    corners.push_back({glob_data[i], glob_data[i+1],
                        glob_data[i+2], glob_data[i+3], glob_data[i+4]});

  for (const auto& c : corners) {
    auto fail = std::find_if(corners.begin(), corners.end(),
                            [c](const std::array<double,5>& C)
                            {
                              return std::fabs(C[0]-c[0]) < 1e-6 &&
                                     std::fabs(C[1]-c[1]) < 1e-6 &&
                                     std::fabs(C[2]-c[2]) < 1e-6 &&
                                    (int)c[3] == (int)C[3] &&
                                    (int)c[4] != (int)C[4];
                            });
    if (fail != corners.end()) {
      std::cerr << "** DomainDecomposition::setup ** Corner node " << (int)c[4]
                << " with coordinates (" << c[0] << "," << c[1] << "," << c[2] << ")"
                << " found with different global ID " << (int)(*fail)[4]
                << " for basis " << (int)c[3] << "."
                << " You have to add vertex / edge connections." << std::endl;
      return false;
    }
  }
#endif

  return true;
}


bool DomainDecomposition::setup(const ProcessAdm& adm, const SIMbase& sim)
{
#ifdef HAVE_MPI
  if (!adm.isParallel())
    IFEM::cout << "Establishing domain decomposition" << std::endl;

#if SP_DEBUG > 1
  IFEM::cout << "  Ghost connections:\n";
  for (const auto& it : ghostConnections) {
    IFEM::cout << "    Interface: master/idx=" << it.master << "/" << it.midx <<
                  ", slave/idx=" << it.slave << "/" << it.sidx <<
                  ", orient=" << it.orient << ", dim=" << it.dim <<
                  ", basis=" << it.basis << std::endl;
  }
#endif
#endif

  sam = dynamic_cast<const SAMpatch*>(sim.getSAM());

  int ok = 1;

  // Establish global node numbers
  if (!calcGlobalNodeNumbers(adm, sim))
    ok = 0;

  // sanity check the established domain decomposition
  if (getMinNode() > getMaxNode()) {
    std::cerr << "**DomainDecomposition::setup ** Process "
              << adm.getProcId() << " owns no nodes." << std::endl;
    ok = 0;
  }

#ifdef HAVE_MPI
  int lok = ok;
  MPI_Allreduce(&lok, &ok, 1, MPI_INT, MPI_SUM, *adm.getCommunicator());
#endif

  if (ok < adm.getNoProcs())
    return false;

  // sanity check all corners of the patches
  if (!sanityCheckCorners(sim))
    return false;

  ok = 1;

  // Establish local equation mappings for each block.
  if (sim.getSolParams() && sim.getSolParams()->getNoBlocks() > 1) {
    IFEM::cout << "\tEstablishing local block equation numbers" << std::endl;
    const LinSolParams& solParams = *sim.getSolParams();
    blocks.clear();
    blocks.resize(solParams.getNoBlocks()+1);

    // Find local equations for each block
    for (size_t i = 0; i < solParams.getNoBlocks(); ++i) {
      // grab DOFs of given type(s)
      blocks[i+1].basis = solParams.getBlock(i).basis;
      blocks[i+1].components = solParams.getBlock(i).comps;
      char dofType = blocks[i+1].basis == 1 ? 'D' : 'P'+blocks[i+1].basis-2;
      if (solParams.getBlock(i).comps != 0) {
        std::set<int> comps = utl::getDigits(solParams.getBlock(i).comps);
        for (auto& c : comps) {
          std::set<int> tmp = sam->getEquations(dofType, c);
          blocks[i+1].localEqs.insert(tmp.begin(), tmp.end());
        }
      } else {
        std::set<int> bases = utl::getDigits(blocks[i+1].basis);
        for (auto& b : bases) {
          int cb = b;
          dofType = cb == 1 ? 'D' : 'P'+cb-2;
          std::set<int> tmp = adm.dd.getSAM()->getEquations(dofType);
          blocks[i+1].localEqs.insert(tmp.begin(), tmp.end());
          // HACK: multipliers always in second block
          // Correct thing to do for average pressure constraint in Stokes.
          if (i == 1) {
            for (size_t n = 1; n <= sim.getPatch(1)->getNoNodes(); ++n) {
              if (sim.getPatch(1)->isLMn(n) && sim.getPatch(1)->getLMType(n) == 'G') {
                int lEq = sam->getEquation(sim.getPatch(1)->getNodeID(n), 1);
                if (lEq > 0)
                  blocks[i+1].localEqs.insert(lEq);
              }
            }
          }
        }
      }

      size_t idx = 1;
      for (auto& it : blocks[i+1].localEqs)
        blocks[i+1].G2LEQ[it] = idx++;
    }
  }

  // Establish global equation numbers for all blocks.
  if (!calcGlobalEqNumbers(adm, sim))
    ok = 0;

#ifdef HAVE_MPI
  if (!adm.isParallel())
    return true;

  std::vector<int> nEqs(blocks.size());
  if (adm.getProcId() == adm.getNoProcs()-1)
    for (size_t i = 0; i < blocks.size(); ++i)
      nEqs[i] = getMaxEq(i);

  MPI_Bcast(&nEqs[0], nEqs.size(), MPI_INT, adm.getNoProcs()-1,
            *adm.getCommunicator());

  for (size_t i = 0; i < blocks.size(); ++i)
    blocks[i].nGlbEqs = nEqs[i];

  IFEM::cout << "\n >>> Domain decomposition summary <<<"
             << "\nNumber of domains     " << adm.getNoProcs();
  IFEM::cout << "\nNumber of equations   " << nEqs[0] << " (" << getMaxEq()-getMinEq()+1 << " on process)";
  for (size_t i = 1; i < blocks.size(); ++i)
    IFEM::cout << "\n  Block " << i << "             " << nEqs[i]  << " (" << getMaxEq(i)-getMinEq(i)+1 << " on process)";
  IFEM::cout << std::endl;

  // sanity check the established domain decomposition
  for (size_t i = 0; i < blocks.size(); ++i) {
    if (getMinEq(i) > getMaxEq(i)) {
      std::cerr << "Process " << adm.getProcId() << " owns no equations";
      if (i > 0)
        std::cerr << " in block " << i;
      std::cerr << "." << std::endl;

      ok = 0;
    }
  }

  lok = ok;
  MPI_Allreduce(&lok, &ok, 1, MPI_INT, MPI_SUM, *adm.getCommunicator());

  if (ok < adm.getNoProcs())
    return false;

#endif

  return true;
}


int DomainDecomposition::getPatchOwner(size_t p) const
{
  auto it = patchOwner.find(p);
  if (it == patchOwner.end())
    return -1;

  return it->second;
}
