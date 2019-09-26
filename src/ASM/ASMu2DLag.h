// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.h
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_MATLAB_H
#define _ASM_S2D_MATLAB_H

#include "ASMs2DLag.h"


/*!
  \brief Driver for assembly of unstructured 2D Lagrange FE models.
  \details This class overrides the methods of its parent class such that
  it does not depend on a surface spline object for geometry discretization.
  It can therefore be used for any unstructured grid read from mesh files.
*/

class ASMu2DLag : public ASMs2DLag
{
public:
  //! \brief Default constructor.
  ASMu2DLag(unsigned char n = 2, unsigned char n_f = 2, char fType = 'm');
  //! \brief Special copy constructor for sharing of FE data.
  ASMu2DLag(const ASMu2DLag& pch, unsigned char nf);
  //! \brief Default copy constructor copying everything.
  ASMu2DLag(const ASMu2DLag& pch);
  //! \brief Empty destructor.
  virtual ~ASMu2DLag() {}

  // Methods for model generation
  // ============================

private:
  //! \brief Creates an instance by reading Matlab commands from a stream.
  bool readMatlab(std::istream& is);

public:
  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);
  //! \brief Generates the finite element topology data for the patch.
  virtual bool generateFEMTopology();
  //! \brief Checks if this patch is empty.
  virtual bool empty() const { return nel == 0; }

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed pre-defined node set.
  virtual const IntVec& getNodeSet(int idx) const;
  //! \brief Returns a named node set for update.
  virtual IntVec& getNodeSet(const std::string& setName);

  using ASMs2DLag::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool, bool);

  // Post-processing methods
  // =======================

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int*) const;

private:
  //! \brief No assignment.
  ASMu2DLag& operator=(const ASMu2DLag&) = delete;

  char fileType; //!< Mesh file format
  typedef std::pair<std::string,IntVec> NodeSet; //!< Named node set container
  std::vector<NodeSet> nodeSets; //!< Pre-defined node sets for Dirichlet BCs
};

#endif
