// $Id$
//==============================================================================
//!
//! \file SIMinput.h
//!
//! \date Sep 24 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Sub-class with functionality for model input and setup.
//!
//==============================================================================

#ifndef _SIM_INPUT_H
#define _SIM_INPUT_H

#include "SIMbase.h"
#include "TopologySet.h"
#include "Interface.h"
#include "Vec3.h"

class ModelGenerator;

namespace LR { struct RefineData; }


/*!
  \brief Sub-class with functionality for model input and setup.
  \details This class extends the SIMbase class with some added functionalities
  for reading model definitions from input file. This also includes reading
  initial conditions from HDF5-files and adaptive refinement. These items
  are put in a separate sub-class to hide them away from the SIMbase class,
  which contains the actual simulation driver (i.e., the computational methods).
*/

class SIMinput : public SIMbase
{
public:
  //! \brief Struct holding information about an initial condition.
  struct ICInfo
  {
    int file_level; //!< The time level for the field in the file
    int geo_level;  //!< The time level for the (adapted) geometry in the file
    char basis;     //!< The basis to inject field into (for mixed)
    char component; //!< Component for field (for functions)
    std::string sim_field;  //!< The name of the field in the %SIM class
    std::string file_basis; //!< Basis for field on file
    std::string file_field; //!< The field name in the file or type of function
    std::string function;   //!< Function expression if given in function form

    //! \brief Default constructor.
    ICInfo() : file_level(-1), geo_level(0), basis(1), component(0) {}
    //! \brief Constructor providing the field name.
    explicit ICInfo(const std::string& f) : file_level(-1), geo_level(0),
                                            basis(1), component(0),
                                            sim_field(f), file_field(f) {}
  };
  typedef std::vector<ICInfo> InitialCondVec; //!< Convenience declaration
  typedef std::vector<unsigned char> CharVec; //!< Convenience declaration
  typedef std::pair<int,Vec3>        IdxVec3; //!< Convenience declaration

protected:
  //! \brief The constructor just forwards to the base class constructor.
  explicit SIMinput(IntegrandBase* itg);

public:
  //! \brief Empty destructor.
  virtual ~SIMinput() {}

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

  //! \brief Reads patch topology from the specified input file.
  //! \param[in] fileName Path to XML file to read patch topology from
  bool readTopologyOnly(const std::string& fileName);

  //! \brief Returns a unique integer code for a Property set.
  //! \param[in] setName Name of the topology set the property is defined on
  //! \param[in] comp The solution components on which the property is applied
  //!
  //! \details The actual Property objects are also created (one for each entity
  //! in the topology set) and their type is set to UNDEFINED. The method
  //! setPropertyType must be used to assign the actual Property type.
  int getUniquePropertyCode(const std::string& setName, int comp = 0);

  //! \brief Defines a vector field property.
  //! \param[in] code The property code to be associated with the property
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] field The vector field representing the physical property
  //! \param[in] pflag Flag for local axis directions (see setPropertyType)
  size_t setVecProperty(int code, Property::Type ptype,
                        VecFunc* field = nullptr, int pflag = -1);

  //! \brief Defines a traction field property.
  //! \param[in] code The property code to be associated with the property
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] field The traction field representing the physical property
  bool setTracProperty(int code, Property::Type ptype,
                       TractionFunc* field = nullptr);

private:
  //! \brief Parses a subelement of the \a geometry XML-tag.
  bool parseGeometryTag(const tinyxml2::XMLElement* elem);
  //! \brief Parses a subelement of the \a boundaryconditions XML-tag.
  bool parseBCTag(const tinyxml2::XMLElement* elem);
  //! \brief Parses the \a initialconditions XML-tag.
  bool parseICTag(const tinyxml2::XMLElement* elem);
  //! \brief Parses a subelement of the \a linearsolver XML-tag.
  bool parseLinSolTag(const tinyxml2::XMLElement* elem);

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Parses a dimension-specific subelement of the \a geometry XML-tag.
  virtual bool parseGeometryDimTag(const tinyxml2::XMLElement* elem) = 0;
  //! \brief Parses the \a periodic XML-tag.
  bool parsePeriodic(const tinyxml2::XMLElement* elem);
  //! \brief Parses a subelement of the \a resultoutput XML-tag.
  virtual bool parseOutputTag(const tinyxml2::XMLElement* elem);
  //! \brief Parses the \a dualfield tag.
  FunctionBase* parseDualTag(const tinyxml2::XMLElement* elem, int ftype = 1);
  //! \brief Parses the "set" attribute of a material XML-tag.
  //! \param[in] elem The XML element extract the set name from
  //! \param[in] mindex Index into problem-dependent material property container
  //! \return The property code to be associated with the material
  int parseMaterialSet(const tinyxml2::XMLElement* elem, int mindex);
  //! \brief Parses the "set" attribute of a refine/raiseorder XML-tag.
  //! \param[in] elem The XML element extract the set name from
  //! \param[in] patches List of patch indices of the specified set
  bool parseTopologySet(const tinyxml2::XMLElement* elem,
                        std::vector<int>& patches) const;
  //! \brief Parses a list of patch indices from an XML-tag.
  //! \param[in] elem The XML element extract the patch indices from
  //! \param[in] patches List of patch indices
  bool parsePatchList(const tinyxml2::XMLElement* elem,
                      std::vector<int>& patches) const;

  //! \brief Creates a set of Property objects.
  //! \param[in] setName Name of the topology set the property is defined on
  //! \param[in] pc The property code to be associated with this set
  bool createPropertySet(const std::string& setName, int pc);

public:
  //! \brief Defines the type of a property set.
  //! \param[in] code The property code to be associated with the property type
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] pindex 0-based index into problem-dependent property container
  //! \param[in] basis 1-based index into basis associated with property
  size_t setPropertyType(int code, Property::Type ptype, int pindex = -1,
                         char basis = 1);

  //! \brief Defines a Neumann boundary condition property by parsing a string.
  //! \param[in] prop The string to parse the property definition from
  //! \param[in] type Additional option defining the type of property definition
  //! \param[in] ndir Direction of the surface traction on the Neumann boundary
  //! \param[in] code The property code to be associated with this property
  bool setNeumann(const std::string& prop, const std::string& type,
                  int ndir, int code);

  //! \brief Finds the set of basis functions with support on a set of elements.
  //! \param[in] elements 0-based element indices
  //! \return 0-based node indices with support on the given elements
  std::vector<int> getFunctionsForElements(const std::vector<int>& elements);

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  bool refine(const LR::RefineData& prm);
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  //! \param[in] sol Vector to interpolate onto refined mesh
  bool refine(const LR::RefineData& prm, Vector& sol);
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  //! \param[in] sol Vectors to interpolate onto refined mesh
  bool refine(const LR::RefineData& prm, Vectors& sol);

  //! \brief Refines the mesh based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  //! \return Positive value on refinement, negative value on error,
  //! zero if no elements are refined
  int refine(const RealFunc& refC, double refTol);

  //! \brief Reads patches from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] whiteSpace For message formatting
  bool readPatches(std::istream& isp, const char* whiteSpace = "");

  //! \brief Connects two patches.
  //! \param[in] mst Master patch
  //! \param[in] slv Slave patch
  //! \param[in] mIdx Boundary index on master patch
  //! \param[in] sIdx Boundary index on slave patch
  //! \param[in] coordCheck If \e true, check for matching coordinates
  //! \param[in] dim Dimensionality of connection
  //! \param[in] thick Thickness of connection
  bool addConnection(int mst, int slv, int mIdx, int sIdx,
                     bool coordCheck = false, int dim = 1, int thick = 1)
  {
    return this->connectPatches(ASM::Interface{mst,slv,mIdx,sIdx,0,dim,0,thick},
                                coordCheck);
  }

  //! \brief Obtain a const reference to model topology.
  const std::vector<ASM::Interface>& getInterfaces() const
  { return myInterfaces; }

protected:
  //! \brief Helper method returning a stream for patch geometry input.
  //! \param[in] tag The XML-tag containing the patch geometry definition
  //! \param[in] patch The value of the \a tag, either a file name or g2 string
  static std::istream* getPatchStream(const char* tag, const char* patch);

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] unf Number of unknowns per basis function for each field
  //! \param[in] whiteSpace For message formatting
  virtual ASMbase* readPatch(std::istream& isp, int pchInd,
                             const CharVec& unf = CharVec(),
                             const char* whiteSpace = "") const = 0;

  //! \brief Reads global node data for a patch from given input stream.
  //! \param[in] isn The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read node data for
  //! \param[in] basis The basis to read node data for (when mixed FEM, 0 = all)
  //! \param[in] oneBased If \e true the read node numbers are assumed
  //! one-based. If \e false they are assumed to be zero-based.
  virtual bool readNodes(std::istream& isn, int pchInd, int basis = 0,
                         bool oneBased = false) { return false; }
  //! \brief Reads node numbers from given input stream.
  //! \param[in] isn The input stream to read from
  virtual void readNodes(std::istream& isn) {}

  //! \brief Instantiates a FEM model generator.
  //! \param[in] geo XML element containing geometry definition
  virtual ModelGenerator* getModelGenerator(const tinyxml2::XMLElement* geo) const = 0;

  //! \brief Connects two patches.
  //! \param[in] ifc Patch interface definition
  //! \param[in] coordCheck If \e false, do not check for matching coordinates
  virtual bool connectPatches(const ASM::Interface& ifc, bool coordCheck) = 0;

  typedef std::map<std::string,std::string> SerializeMap; //!< Convenience type

  //! \brief Writes current basis to a serialization container.
  //! \param data Container for serialized data
  bool saveBasis(SerializeMap& data) const;

  //! \brief Restores the basis from a serialization container.
  //! \param[in] data Container for serialized data
  bool restoreBasis(const SerializeMap& data);

public:
  //! \brief Handles application restarts by reading a serialized basis.
  //! \param[in] restartFile Path to file to read restart basis from
  //! \param[in] restartStep Index of the time step to read restart basis for
  //! \return One-based index of the first time step to solve after restart.
  //! If zero, no restart specified. If one, no serialized basis stored.
  //! If negative, read failure.
  int restartBasis(const std::string& restartFile, int restartStep);

  //! \brief Creates the computational FEM model from the spline patches.
  //! \param[in] resetNumb If \e 'y', start element and node numbers from zero
  virtual bool createFEMmodel(char resetNumb = 'y');

  //! \brief Creates the computational FEM model by copying the given patches.
  virtual void clonePatches(const PatchVec&, const std::map<int,int>&) {}

  //! \brief Checks whether a named initial condition is present.
  virtual bool hasIC(const std::string& name) const;

  //! \brief Sets the initial conditions.
  //! \param fieldHolder The SIM-object to inject the initial conditions into,
  //! nullptr means this.
  bool setInitialConditions(SIMdependency* fieldHolder = nullptr);

  //! \brief Deserialization support (for simulation restart).
  virtual bool deSerialize(const SerializeMap&);

  //! \brief Returns reference to a named topology entity.
  const TopEntity& getEntity(const std::string& name) const;
  //! \brief Returns access to a named topology entity (for model generators).
  TopEntity& topology(const std::string& name) { return myEntitys[name]; }
  //! \brief Returns the whole topology set container (for testing only).
  const TopologySet& getTopology() const { return myEntitys; }

  //! \brief Returns the node numbers associated with a topological item.
  bool getTopItemNodes(const TopItem& titem, std::vector<int>& glbNodes) const;

  //! \brief Returns pointer to an indexed discrete point.
  IdxVec3* getDiscretePoint(int idx);

  //! \brief Returns the element-to-element connectivities.
  virtual std::vector<std::vector<int>> getElmConnectivities() const;

private:
  //! \brief Sets initial conditions from a file.
  //! \param fieldHolder The SIM-object to inject the initial conditions into
  //! \param[in] fileName Path to file to read the initial conditions from
  //! \param[in] info Initial condition information
  bool setInitialCondition(SIMdependency* fieldHolder,
                           const std::string& fileName,
                           const InitialCondVec& info);

protected:
  ModelGenerator* myGen; //!< Model generator

  TopologySet myEntitys; //!< Set of named topological entities

  std::vector<IdxVec3> myTopPts; //!< Discrete points not belonging to any patch

  std::vector<ASM::Interface> myInterfaces; //!< Topology interface descriptions

  std::map<std::string,InitialCondVec> myICs; //!< Initial condition definitions

private:
  bool isReading; //!< If \e true, we are reading multiple patch tags
  bool isSaved;   //!< If \e true, the basis has been serialized
};

#endif
