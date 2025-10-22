// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.h
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D %Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_LAG_H
#define _ASM_U2D_LAG_H

#include "ASMs2DLag.h"
#include "ASMutils.h"


/*!
  \brief Driver for assembly of unstructured 2D %Lagrange FE models.
  \details This class overrides the methods of its parent class such that
  it does not depend on a surface spline object for geometry discretization.
  It can therefore be used for any unstructured grid read from mesh files.
*/

class ASMu2DLag : public ASMs2DLag
{
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ASMs2DLag::BasisFunctionCache
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    //! \param pch Patch the cache is for
    BasisFunctionCache(const ASMu2DLag& pch)
      : ASMs2DLag::BasisFunctionCache(pch) {}

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

    //! \brief No integration point parameters for unstructured patches.
    double getParam(int, size_t, size_t, bool) const override { return 0.0; }

  protected:
    //! \brief No integration point parameters for unstructured patches.
    void setupParameters() override {}
  };

public:
  //! \brief Default constructor.
  ASMu2DLag(unsigned char n = 2, unsigned char n_f = 2, char fType = 'm');
  //! \brief Special copy constructor for sharing of FE data.
  ASMu2DLag(const ASMu2DLag& pch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMu2DLag(const ASMu2DLag& pch);
  //! \brief Empty destructor.
  virtual ~ASMu2DLag() {}

  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);
  //! \brief Generates the finite element topology data for the patch.
  virtual bool generateFEMTopology();
  //! \brief Checks if this patch is empty.
  virtual bool empty() const { return nel == 0; }

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed predefined node set.
  virtual const IntVec& getNodeSet(int iset) const;
  //! \brief Checks if node \a inod is within predefined node set \a iset.
  virtual bool isInNodeSet(int iset, int inod) const;
  //! \brief Defines a node set by parsing a list of node numbers.
  virtual int parseNodeSet(const std::string& setName, const char* cset);
  //! \brief Defines a node set by parsing a 3D bounding box.
  virtual int parseNodeBox(const std::string& setName, const char* bbox);
  //! \brief Adds a node to a named node set.
  //! \param[in] setName Name of the node set to receive the node
  //! \param[in] inod Identifier of the node to add to the set
  //! \param[in] extId If \e true, \a inod is the external node Id,
  //! otherwise it is a one-based internal node index
  //! \return Node set index associated with \a setName
  int addToNodeSet(const std::string& setName, int inod, bool extId = false);
  //! \brief Returns the name of an indexed predefined node set.
  bool getNodeSet(int iset, std::string& name) const;

  //! \brief Returns (1-based) index of a predefined element set in the patch.
  virtual int getElementSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed predefined element set.
  virtual const IntVec& getElementSet(int iset) const;
  //! \brief Returns the name of an indexed predefined element set.
  virtual bool getElementSet(int iset, std::string& name) const;
  //! \brief Checks if element \a iel is within predefined element set \a iset.
  virtual bool isInElementSet(int iset, int iel) const;
  //! \brief Defines an element set by parsing a list of element numbers.
  virtual int parseElemSet(const std::string& setName, const char* cset);
  //! \brief Defines an element set by parsing a 3D bounding box.
  virtual int parseElemBox(const std::string& setName,
                           const std::string& unionSet, const char* bbox);
  //! \brief Adds an element to a named element set.
  //! \param[in] setName Name of the element set to receive the element
  //! \param[in] iel Identifier of the element to add to the set
  //! \param[in] extId If \e true, \a iel is the external element Id,
  //! otherwise it is a one-based internal element index
  //! \return Element set index associated with \a setName
  int addToElemSet(const std::string& setName, int iel, bool extId = false);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary node set
  //! \param nodes Array of node numbers
  //! \param[in] local If \e true, return patch-local numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int, int, int, bool local) const;

  using ASMs2DLag::integrate;
  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  using ASMs2DLag::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool silence,
                                    bool separateGroup1noded);

  // Post-processing methods
  // =======================

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int*) const;

protected:
  //! \brief Generate thread groups using multi-coloring.
  void generateThreadGroupsMultiColored(bool silence, bool separateGroup1Noded);

  bool swapNode34; //!< If \e true, element nodes 3 and 4 should be swapped

  std::vector<ASM::NodeSet> nodeSets; //!< Node sets for Dirichlet BCs
  std::vector<ASM::NodeSet> elemSets; //!< Element sets for properties

private:
  char fileType; //!< Mesh file format
};

#endif
