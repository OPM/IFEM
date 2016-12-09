// $Id$
//==============================================================================
//!
//! \file DomainDecomposition.h
//!
//! \date Feb 23 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Domain decomposition related partitioning for structured models.
//!
//==============================================================================

#ifndef _DOMAIN_DECOMPOSITION_H
#define _DOMAIN_DECOMPOSITION_H

#include <map>
#include <set>
#include <vector>
#include <cstddef>

class ASMbase;
class ProcessAdm;
class SAMpatch;
class SIMbase;


/*!
  \brief Class containing domain decomposition related partitioning.
*/

class DomainDecomposition
{
public:
  //! \brief Struct defining a domain interface.
  struct Interface {
    int master; //!< Master patch (global number).
    int slave;  //!< Slave patch (global number).
    int midx;   //!< Index of boundary on master.
    int sidx;   //!< Index of boundary on slave.
    int orient; //!< Orientation.
    int dim;    //!< Dimension of boundary.
    int basis;  //!< Basis of boundary.
    int thick;  //!< Thickness of connection.
  };

  //! \brief Functor to order ghost connections.
  class SlaveOrder {
    public:
      //! \brief The constructor initializes the DomainDecomposition reference.
      SlaveOrder(const DomainDecomposition& dd_) : dd(dd_) {}
      //! \brief Hide ill-formed default assignment operator.
      SlaveOrder& operator=(const SlaveOrder&) { return *this; }
      //! \brief Compare interfaces.
      bool operator()(const Interface& A, const Interface& B) const
      {
        // order by master owner id first
        if (dd.getPatchOwner(A.master) != dd.getPatchOwner(B.master))
          return dd.getPatchOwner(A.master) < dd.getPatchOwner(B.master);

        // then by lower slave owner id
        if (dd.getPatchOwner(A.slave) != dd.getPatchOwner(B.slave))
          return dd.getPatchOwner(A.slave) < dd.getPatchOwner(B.slave);

        // then by slave id
        if (A.slave != B.slave)
          return A.slave < B.slave;

        // then by dim
        if (A.dim != B.dim)
          return A.dim < B.dim;

        // finally by index on master
        return A.midx < B.midx;
      }
    protected:
      const DomainDecomposition& dd;
  };

  std::set<Interface, SlaveOrder> ghostConnections; //!< Connections to other processes.

  //! \brief Default constructor.
  DomainDecomposition() : ghostConnections(SlaveOrder(*this)), blocks(1) {}

  //! \brief Setup domain decomposition.
  bool setup(const ProcessAdm& adm, const SIMbase& sim);

  //! \brief Obtain local subdomains for an equation block.
  //! \param nx Number of domains in x
  //! \param ny Number of domains in y
  //! \param nz Number of domains in z
  //! \param overlap Overlap
  //! \param block Block to obtain equations for
  //! \return Vector with equations on each subdomain
  std::vector<std::set<int>> getSubdomains(int nx, int ny, int nz, int overlap, size_t block) const;

  //! \brief Calculates subdomains with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] nel3 Number of knot-spans in third parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] g3 Number of subdomains in third parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  //! \details nel values determine the dimensionality.
  static std::vector<std::vector<int>> calcSubdomains(size_t nel1, size_t nel2, size_t nel3,
                                                      size_t g1, size_t g2, size_t g3, size_t overlap);

  //! \brief Get first equation owned by this process.
  int getMinEq(size_t idx = 0) const { return blocks[idx].minEq; }
  //! \brief Get last equation owned by this process.
  int getMaxEq(size_t idx = 0) const { return blocks[idx].maxEq; }
  //! \brief Get total number of equations in a block.
  int getNoGlbEqs(size_t idx = 0) const { return blocks[idx].nGlbEqs; }
  //! \brief Get first node owned by this process.
  int getMinNode() const { return minNode; }
  //! \brief Get last node owned by this process.
  int getMaxNode() const { return maxNode; }
  //! \brief Get first DOF owned by this process.
  int getMinDOF() const { return minDof; }
  //! \brief Get last DOF owned by this process.
  int getMaxDOF() const { return maxDof; }

  //! \brief Set owner for a patch.
  void setPatchOwner(size_t p, size_t owner) { patchOwner[p] = owner; }

  //! \brief Get process owning patch p.
  int getPatchOwner(size_t p) const;

  //! \brief Get global equation number
  //! \param lEq Local equation number
  //! \param idx Block to get index for
  int getGlobalEq(int lEq, size_t idx=0) const;

  //! \brief Obtain local-to-global equation mapping.
  const std::vector<int>& getMLGEQ(size_t idx = 0) const { return blocks[idx].MLGEQ; }

  //! \brief Obtain local-to-global equation mapping.
  //! \param idx Block equation index (global block not included).
  const std::set<int>& getBlockEqs(size_t idx) const { return blocks[idx+1].localEqs; }

  //! \brief Obtain global-to-local equation mapping.
  //! \param idx Block equation index (global block not included).
  const std::map<int,int>& getG2LEQ(size_t idx) const { return blocks[idx+1].G2LEQ; }

  //! \brief Obtain local-to-global node mapping.
  const std::vector<int>& getMLGN() const { return MLGN; }

  //! \brief Returns associated SAM
  const SAMpatch* getSAM() const { return sam; }

  //! \brief Returns number of matrix blocks.
  size_t getNoBlocks() const { return blocks.size()-1; }

private:
  //! \brief Calculates a 1D partitioning with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  static std::vector<std::vector<int>> calcSubdomains1D(size_t nel1, size_t g1, size_t overlap);

  //! \brief Calculates 2D subdomains with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  static std::vector<std::vector<int>> calcSubdomains2D(size_t nel1, size_t nel2,
                                                        size_t g1, size_t g2, size_t overlap);

  //! \brief Calculates 3D subdomains with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] nel2 Number of knot-spans in third parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] g3 Number of subdomains in third parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  static std::vector<std::vector<int>> calcSubdomains3D(size_t nel1, size_t nel2, size_t nel3,
                                                        size_t g1, size_t g2, size_t g3, size_t overlap);


  //! \brief Setup equation numbers for all blocks on a boundary.
  //! \param sim Simulator with patches and linear solver block information
  //! \param pidx Patch index
  //! \param lidx Boundary index on patch
  //! \param cbasis If non-empty, bases to connect
  //! \param thick Thickness of connection (subdivisions)
  std::vector<int> setupEquationNumbers(const SIMbase& sim,
                                        int pidx, int lidx,
                                        const std::set<int>& cbasis,
                                        int dim, int thick);

  //! \brief Setup node numbers for all bases on a boundary.
  //! \param basis Bases to grab nodes for
  //! \param lNodes Resulting nodes
  //! \param cbasis Bases to connect
  //! \param pch Patch to obtain nodes for
  //! \param dim Dimension of boundary to obtain nodes for
  //! \param lidx Local index of boundary to obtain nodes for
  //! \param thick Thickness of connection (subdivisions)
  void setupNodeNumbers(int basis, std::vector<int>& lNodes,
                        std::set<int>& cbasis,
                        const ASMbase* pch,
                        int dim, int lidx, int thick);

  //! \brief Calculate the global node numbers for given finite element model.
  bool calcGlobalNodeNumbers(const ProcessAdm& adm, const SIMbase& sim);

  //! \brief Calculate the global equation numbers for given finite element model.
  bool calcGlobalEqNumbers(const ProcessAdm& adm, const SIMbase& sim);

  //! \brief Sanity check model.
  //! \details Collects the corners of all patches in the model and make sure nodes
  //!          with matching coordinates have the same global node ID.
  bool sanityCheckCorners(const SIMbase& sim);

  std::map<int,int> patchOwner; //!< Process that owns a particular patch

  //! \brief Struct with information per matrix block.
  struct BlockInfo {
    int basis;      //!< Bases for block
    int components; //!< Components in block
    std::vector<int> MLGEQ; //!< Process-local-to-global equation numbers for block.
    int minEq; //!< First equation we own in block.
    int maxEq; //!< Last equation we own in block.
    int nGlbEqs; //!< Total matrix size
    std::set<int> localEqs; //!< Local equations belonging to the block.
    std::map<int,int> G2LEQ; //!< Maps from local total matrix index to local block index
  };

  std::vector<int> MLGN; //!< Process-local-to-global node numbers
  std::vector<BlockInfo> blocks; //!< Equation mappings for all matrix blocks.
  int minDof; //!< First DOF we own
  int maxDof; //!< Last DOF we own
  int minNode; //!< First node we own
  int maxNode; //!< Last node we own

  const SAMpatch* sam; //!< The assembly handler the DD is constructed for.
};

#endif
