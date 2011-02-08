// $Id: LinearEl.h,v 1.19 2009-10-30 12:12:40 kmo Exp $
//==============================================================================
//!
//! \file LinearEl.h
//!
//! \date Jan 27 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solver for linear elasticity problems using NURBS-based FEM.
//!
//==============================================================================

#ifndef _LINEAR_EL_H
#define _LINEAR_EL_H

#include "Vec3.h"
#include "Function.h"
#include "PressureLoad.h"
#include "LinEqSystem.h"
#include <map>


class VolumePatch;
class LocalSystem;
class SAMSpline;
class VTF;


/*!
  \brief Struct for storage of data associated with one mode shape.
*/

struct Mode
{
  double eigVal; //!< Eigenvalue associated with this mode
  Vector eigVec; //!< Eigenvector associated with this mode
  // \brief Default constructor setting \a eigVal to zero.
  Mode() { eigVal = 0.0; }
};


/*!
  \brief Result container for passing results to GLview export module.
*/

struct Result
{
  Vector            load;  //!< The external load vector
  Vector            displ; //!< The displacement vector due to the external load
  Matrix            norms; //!< Element norms
  std::vector<Mode> modes; //!< Mode shapes
  bool freq; //!< Indicates whether eigenvalues are frequencies or not
  // \brief Default constructor setting \a freq to \e false.
  Result() { freq = false; }
};


/*!
  \brief Driver class for the NURBS-based linear elasticity solver.
  \details The class incapsulates data and methods for solving linear elasticity
  problems using NURBS-based finite elements.
*/

class LinearEl
{
public:
  //! \brief The constructor generates the model from the given input file.
  //! \param[in] fileName Name of input file with model description
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  //! \param[in] free If \e true, any specified boundary conditions are ignored
  LinearEl(const char* fileName = 0, bool checkRHS = false, bool free = false);
  //! \brief The destructor frees the dynamically allocated objects.
  ~LinearEl();

  //! \brief Reads model data from the specified input file \a *fileName.
  bool read(const char* fileName, bool checkRHS = false, bool free = false);

  //! \brief Performs some preprocessing tasks on the finite element model.
  //! \details This method should be invoked after \a readModel and before
  //! any of the assembleK.... methods.
  //! The main purpose of this method is to fill the SAMSpline object with data.
  bool preprocess(const std::vector<int>& ignoredPatches = std::vector<int>(),
		  bool fixDup = false);

  //! \brief Administers assembly of system stiffness matrix and load vector.
  //! \param[in] solver Which equation solver (matrix format) to use
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[out] R The external load vector in DOF-order (for visualization)
  bool assembleKandR(SystemMatrix::Type solver, int nGauss, Vector* R);
  //! \brief Administers assembly of the system stiffness- and mass matrices.
  //! \param[in] solver Which equation solver (matrix format) to use
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  bool assembleKandM(SystemMatrix::Type solver, int nGauss);
  //! \brief Administers assembly of the system stiffness matrices.
  //! \param[in] solver Which equation solver (matrix format) to use
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[in] dis Solution vector in DOF-order (for geometric stiffness)
  bool assembleKandKg(SystemMatrix::Type solver, int nGauss, const Vector& dis);
  //! \brief Administers assembly of the system stiffness matrix.
  //! \param[in] solver Which equation solver (matrix format) to use
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  bool assembleKonly(SystemMatrix::Type solver, int nGauss);

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Solution vector, displacements at nodal points
  bool solve(Vector& solution);

  //! \brief Integrates and prints out some solution norm quantities.
  //! \details If an analytical solution is provided, norms of the exact
  //! error in the solution is computed as well.
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[in] dis Solution vector in DOF-order, displacements at nodal points
  bool solutionNorms(Matrix& eNorm, int nGauss, const Vector& dis);

  //! \brief Performs a generalized eigenvalue analysis of the assembled system.
  //! \param[in] iop Which eigensolver method to use
  //! \param[in] nev Number of eigenvalues/vector (see ARPack documentation)
  //! \param[in] ncv Number of Arnoldi vectors (see ARPack documentation)
  //! \param[in] shift Eigenvalue shift
  //! \param[out] solution Computed eigenvalues and associated eigenvectors
  bool modes(int iop, int nev, int ncv, double shift,
	     std::vector<Mode>& solution);

  //! \brief Writes the global grid and boundary conditions to a file.
  //! \param[in] inputFile File name used to construct the grid file name from
  //! \param[in] n         Number of FE nodes over each knot-span
  //! \param[in] nenod     Number of nodes in each element to generate
  bool writeGlobalGrid(const char* inputFile, const int* n,
		       int nenod = 8) const;

  //! \brief Writes a VTF-file with the geometry and solution fields.
  //! \param[in] inputFile File name used to construct the VTF-file name from
  //! \param[in] solution  The solution vector(s) to output
  //! \param[in] nViz      Number of visualization points over a knot-span
  //! \param[in] format    Format of VTF-file (0=ASCII, 1=BINARY)
  //! \param[in] debug     If \e true, output some additional debug quantities
  //!
  //! \details The NURBS patches are tesselated into linear hexahedrons with
  //! a fixed number of HEX8-elements within each knot-span of non-zero length.
  //! The solution fields (displacements and stresses) are then evaluated at the
  //! nodal points of the generated HEX8-mesh and written to the VTF-file as
  //! vector and scalar result fields.
  bool writeGlv(const char* inputFile, const Result& solution,
		const int* nViz, int format = 0, bool debug = false) const;

  //! \brief Writes the geometry to g2-file for external viewing.
  //! \param[in] g2file Filename for GoTools output file
  void dumpGeometry(const char* g2file) const;

  //! \brief Writes the solution vector to files for external viewing.
  //! \param[in] solfile Filename prefix
  //! \param[in] solution The global solution vector in DOF-order
  //!
  //! \details The solution vector is written to the ASCII file \a solfile.vec
  //! in the same order as the control points are ordered in the NURBS patches.
  //! In addition, each component is written to files \a solfile.u, \a solfile.v
  //! and \a solfile.w respectively.
  void dumpSolution(const char* solfile, const Vector& solution) const;

  //! \brief Writes the system matrices to files for debugging.
  //! \param[in] Kfile Filename for stiffness matrix
  //! \param[in] Mfile Filename for mass matrix
  //! \param[in] Rfile Filename for load vector
  void dumpMat(const char* Kfile,
	       const char* Mfile = 0,
	       const char* Rfile = 0) const;

protected:
  //! \brief Creates the computational FEM model.
  bool createFEMmodel();

  //! \brief Writes the grid geometry to the VTF-file.
  //! \param vtf The VTF-file object to receive the geometry data
  //! \param[in] nViz Number of visualization points over each knot-span
  bool writeGlvG(VTF& vtf, const int* nViz) const;

  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvT(VTF& vtf, int iStep, int& nBlock) const;

  //! \brief Writes boundary condition codes as scalar fields to the VTF-file.
  //! \param vtf The VTF-file object to receive the boundary condition data
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param nBlock Running result block counter
  //! \param[in] debug If \e true, output some additional debug quantities
  bool writeGlvBC(VTF& vtf, const int* nViz, int& nBlock, bool debug) const;

  //! \brief Writes the load vector for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the load vector
  //! \param[in] load The load vector to output
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvR(VTF& vtf, const Vector& load,
		 const int* nViz, int iStep, int& nBlock) const;

  //! \brief Writes solution fields for a given load/time step to the VTF-file.
  //! \details If an analytical solution is provided, the exact stress fields
  //! are written to the VTF-file as well.
  //! \param vtf The VTF-file object to receive the solution fields
  //! \param[in] solution The solution vector to derive the result fields from
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvS(VTF& vtf, const Vector& solution,
		 const int* nViz, int iStep, int& nBlock) const;

  //! \brief Writes an eigenvector and associated eigenvalue to the VTF-file.
  //! \details The eigenvalue is used as a label on the step state info
  //! that is associated with the eigenvector.
  //! \param vtf The VTF-file object to receive the eigenvector
  //! \param[in] mode The eigenvector to output
  //! \param[in] freq \e true if the eigenvalue is a frequency
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] iMode Mode shape identifier
  //! \param nBlock Running result block counter
  bool writeGlvM(VTF& vtf, const Mode& mode, bool freq,
		 const int* nViz, int iMode, int& nBlock) const;

  //! \brief Writes element norms for a given load/time step to the VTF-file.
  //! \details This method can be used only when the number of visualization
  //! points over each knot-span equals 2 (that is, no additonal points).
  //! \param vtf The VTF-file object to receive the element norms
  //! \param[in] norms The element norms to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvN(VTF& vtf, const Matrix& norms, int iStep, int& nBlock) const;

private:
  std::vector<VolumePatch*> model;   //!< The actual NURBS/spline model
  std::vector<PLoad>        load;    //!< Surface pressure loads
  Vec3                      gravity; //!< Gravitation vector
  LocalSystem*              rCS;     //!< Local coordinate system for results
  std::map<Vec3,Vec3>       trac;    //!< Evaluated surface tractions
  TensorFunc*               asol;    //!< Analytical stress field
  SAMSpline*                sam;     //!< Data for FE assembly management
  LinEqSystem               sys;     //!< The linear equation system
};

#endif
