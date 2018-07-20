// $Id$
//==============================================================================
//!
//! \file SIMoutput.h
//!
//! \date Sep 05 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Sub-class with functionality for result output to VTF and terminal.
//!
//==============================================================================

#ifndef _SIM_OUTPUT_H
#define _SIM_OUTPUT_H

#include "SIMinput.h"
#include "Vec3.h"

class VTF;

typedef std::pair<Vec3,double> PointValue; //!< Convenience type


/*!
  \brief Sub-class with additional functionality for result output.
  \details This class extends the SIMbase class with some added functionalities
  for dumping simulation results to VTF-file, and terminal printout. These items
  are put in a separate sub-class here to hide them away from the SIMbase class,
  which contains the main simulation driver.
*/

class SIMoutput : public SIMinput
{
protected:
  //! \brief The constructor just forwards to the base class constructor.
  explicit SIMoutput(IntegrandBase* itg);

public:
  //! \brief The destructor frees the dynamically allocated VTF object.
  virtual ~SIMoutput();

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

  using SIMinput::parse;
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

protected:
  //! \brief Parses a subelement of the \a resultoutput XML-tag.
  virtual bool parseOutputTag(const TiXmlElement* elem);

public:
  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //! \param[in] doClear If \e true, clear geometry block if \a inpFile is null
  //!
  //! \details The spline patches are tesselated into linear finite elements
  //! with a fixed number of elements within each knot-span of non-zero length.
  //! The solution fields are then evaluated at the nodal points of the
  //! generated FE mesh and written to the VTF-file as vector and scalar fields
  //! by the other \a writeGlv* methods.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear = true);

  //! \brief Writes additional, problem-specific, results to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  virtual bool writeGlvA(int& nBlock, int iStep = 1) const { return true; }

  //! \brief Writes boundary conditions as scalar fields to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  bool writeGlvBC(int& nBlock, int iStep = 1) const;

  //! \brief Writes boundary tractions for a given time step to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  bool writeGlvT(int iStep, int& geoBlk, int& nBlock) const;

  //! \brief Writes a vector field for a given load/time step to the VTF-file.
  //! \param[in] vec The vector field to output (nodal values in DOF-order)
  //! \param[in] fieldName Name identifying the vector field
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Result block ID number
  //! \param[in] ncmp Number of components in vector field
  bool writeGlvV(const Vector& vec, const char* fieldName,
                 int iStep, int& nBlock, int idBlock = 2, int ncmp = 0) const;

  //! \brief Writes a scalar field for a given load/time step to the VTF-file.
  //! \param[in] scl The scalar field to output (nodal values)
  //! \param[in] fieldName Name identifying the scalar field
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Result block ID number
  bool writeGlvS(const Vector& scl, const char* fieldName,
                 int iStep, int& nBlock, int idBlock = 2) const;

  //! \brief Writes solution fields for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  bool writeGlvS(const Vector& psol, int iStep, int& nBlock, double time = 0.0,
                 const char* pvecName = nullptr, int idBlock = 10,
                 int psolComps = 0);

  //! \brief Writes primary solution for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  //! \param[in] scalarOnly If \e true, write vector as scalar components only
  int writeGlvS1(const Vector& psol, int iStep, int& nBlock, double time = 0.0,
                 const char* pvecName = nullptr, int idBlock = 10,
                 int psolComps = 0, bool scalarOnly = false);

  //! \brief Writes secondary solution for a load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  bool writeGlvS2(const Vector& psol, int iStep, int& nBlock, double time = 0.0,
                  int idBlock = 20, int psolComps = 0);

  //! \brief Evaluates the secondary solution for a given load/time step.
  //! \param[in] psol Primary solution vector
  //! \param[in] time Load/time step parameter
  //! \param[in] psolComps Optional number of primary solution components
  //!
  //! \details This method only evaluates the solutions fields, and does not
  //! return any data. It is used only for load/time steps that are not saved
  //! when the solution has to be evaluated at every increment in any case to
  //! ensure consistency (e.g., when material models with history variables
  //! are in use).
  bool eval2ndSolution(const Vector& psol, double time, int psolComps = 0);

  //! \brief Writes projected solutions for a given time step to the VTF-file.
  //! \param[in] ssol Secondary solution vector (control point values)
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] prefix Common prefix for the field components
  //! \param[in] maxVal Optional array of maximum values
  bool writeGlvP(const Vector& ssol, int iStep, int& nBlock,
                 int idBlock = 100, const char* prefix = "Global projected",
                 std::vector<PointValue>* maxVal = nullptr);

  //! \brief Evaluates the projected solution for a given load/time step.
  //! \param[in] ssol Secondary solution vector (control point values)
  //! \param[in] maxVal Array of maximum values
  //!
  //! \details This method only evaluates/updates the maximum values of the
  //! secondary solution fields (i.e. same as writeGlvP but without VTF output).
  bool evalProjSolution(const Vector& ssol, std::vector<PointValue>& maxVal);

  //! \brief Writes a mode shape to the VTF-file.
  //! \param[in] mode The mode shape eigenvector and associated eigenvalue
  //! \param[in] freq \e true if the eigenvalue is a frequency
  //! \param nBlock Running result block counter
  //!
  //! \details The eigenvalue is used as a label on the step state info.
  bool writeGlvM(const Mode& mode, bool freq, int& nBlock);

  //! \brief Writes element field for a given load/time step to the VTF-file.
  //! \param[in] field The element field to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] name Name of field
  bool writeGlvE(const Vector& field, int iStep, int& nBlock, const char* name);

  //! \brief Writes element norms for a given load/time step to the VTF-file.
  //! \param[in] norms The element norms to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] prefix Prefices for projected solutions
  //! \param[in] idBlock Starting value of result block numbering
  bool writeGlvN(const Matrix& norms, int iStep, int& nBlock,
                 const char** prefix = nullptr, int idBlock = 200);

  //! \brief Writes a scalar function to the VTF-file.
  //! \param[in] f The function to output
  //! \param[in] fname Name of the function
  //! \param[in] iStep Load/time step identifier
  //! \param[in] idBlock Starting value of result block numbering
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] time Load/time step parameter
  bool writeGlvF(const RealFunc& f, const char* fname,
                 int iStep, int& nBlock, int idBlock = 50, double time = 0.0);

  //! \brief Writes load/time step info to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] value Load parameter or time of the step
  //! \param[in] itype Type identifier of the step
  bool writeGlvStep(int iStep, double value = 0.0, int itype = 0);

  //! \brief Closes the current VTF-file.
  void closeGlv();

  //! \brief Returns the current VTF-file object.
  VTF* getVTF() const { return myVtf; }
  //! \brief Defines the VTF-file for subsequent results output.
  void setVTF(VTF* vtf) { myVtf = vtf; }

  //! \brief Initializes the geometry block counter.
  void setStartGeo(int gID) { myGeomID = gID; }

  //! \brief Dumps the FE model to Matlab format.
  //! \param os Output stream to write the grid data to
  //! \param[in] name Name of the Matlab function returning the grid
  //! \param[in] sets Names of nodal sets to be printed out in addition
  //! \param[in] scale Scaling factor
  bool dumpMatlabGrid(std::ostream& os, const std::string& name = "IFEM_Mesh",
                      const std::vector<std::string>& sets = {},
                      double scale = 1.0) const;

  //! \brief Dumps the (possibly refined) spline geometry in g2-format.
  //! \param os Output stream to write the geometry data to
  bool dumpGeometry(std::ostream& os) const;

  //! \brief Dumps the primary solution in ASCII format for inspection.
  //! \param[in] psol Primary solution vector
  //! \param os Output stream to write the solution data to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpPrimSol(const Vector& psol, utl::LogStream& os,
                   bool withID = true) const;
  //! \brief Dumps the entire solution in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param os Output stream to write the solution data to
  bool dumpSolution(const Vector& psol, utl::LogStream& os) const;
  //! \brief Dumps solution results at specified points in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(const Vector& psol, double time, utl::LogStream& os,
                   bool formatted = false, std::streamsize precision = 3) const;
  //! \brief Dumps vector solution at specified points in ASCII format.
  //! \param[in] vsol Solution vector
  //! \param[in] fname Name of vector field
  //! \param os Output stream to write the solution data to
  //! \param[in] precision Number of digits after the decimal point
  bool dumpVector(const Vector& vsol, const char* fname,
                  utl::LogStream& os, std::streamsize precision = 3) const;
  //! \brief Dumps additional problem-specific results in ASCII format.
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] precision Number of digits after the decimal point
  virtual void dumpMoreResults(double time, utl::LogStream& os,
                               std::streamsize precision = 3) const {}

  //! \brief Saves point results to output file for a given time step.
  //! \param[in] psol Primary solution vector
  //! \param[in] time Load/time step parameter
  //! \param[in] step Load/time step counter
  bool savePoints(const Vector& psol, double time, int step) const;

  //! \brief Sets the file name for result point output.
  //! \param[in] filename The file name prefix (optionally with extension)
  //! \param[in] dumpCoord If \e true, write point coordinates to separate file
  void setPointResultFile(const std::string& filename, bool dumpCoord = false);
  //! \brief Checks whether point result files have been defined or not.
  bool hasPointResultFile() const;
  //! \brief Checks whether result points have been defined or not.
  bool hasResultPoints() const { return !myPoints.empty(); }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] norms The norm values
  //! \param[in] w Total number of characters in the norm labels
  virtual void printNorms(const Vectors& norms, size_t w = 36) const;
  //! \brief Prints a norm group to the log stream (app-specific).
  virtual void printNormGroup(const Vector&, const Vector&,
                              const std::string&) const {}

  //! \brief Returns the reference norm to base mesh adaptation upon.
  virtual double getReferenceNorm(const Vectors& gNorm, size_t adaptor) const;

  //! \brief Returns the global effectivity index.
  virtual double getEffectivityIndex(const Vectors& gNorm,
                                     size_t idx, size_t inorm) const;

  //! \brief Serialization support.
  virtual bool serialize(std::map<std::string,std::string>&) { return false; }

protected:
  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints();

private:
  //! \brief Private helper to write out scalar fields to VTF-file.
  bool writeScalarFields(const Matrix& field, int geomID,
                         int& nBlock, std::vector< std::vector<int> >& sID,
                         size_t* i = nullptr);

  //! \brief Struct defining a result sampling point.
  struct ResultPoint
  {
    unsigned char npar;  //!< Number of parameters
    size_t        patch; //!< Patch index [1,nPatch]
    int           inod;  //!< Local node number of the closest node
    double        u[3];  //!< Parameters of the point (u,v,w)
    Vec3          X;     //!< Spatial coordinates of the point
    // \brief Default constructor.
    ResultPoint() : npar(0), patch(1), inod(0) { u[0] = u[1] = u[2] = 0.0; }
  };

  typedef std::vector<ResultPoint> ResPointVec; //!< Result point container

  //! \brief Preprocesses a result sampling point group.
  //! \param ptFile Name of file that these result points are dumped to
  //! \param points Group of result points that are dumped to the given file
  void preprocessResPtGroup(std::string& ptFile, ResPointVec& points);

  //! \brief Dumps solution results at the given points in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] gPoints Group of result points to write solution data for
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(const Vector& psol, double time,
                   utl::LogStream& os, const ResPointVec& gPoints,
                   bool formatted, std::streamsize precision) const;

  typedef std::pair<std::string,ResPointVec> ResPtPair; //!< Result point group

  std::vector<ResPtPair> myPoints; //!< User-defined result sampling points
  int                    myPrec;   //!< Output precision for result sampling

  int  myGeomID; //!< VTF geometry block ID for the first patch
  VTF* myVtf;    //!< VTF-file for result visualization
};

#endif
