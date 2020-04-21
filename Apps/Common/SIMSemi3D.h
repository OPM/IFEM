// $Id$
//==============================================================================
//!
//! \file SIMSemi3D.h
//!
//! \date Jun 5 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for plane-decoupled solvers.
//!
//==============================================================================

#ifndef _SIM_SEMI3D_H
#define _SIM_SEMI3D_H

#include "SIMadmin.h"
#include "SIMdependency.h"
#include "SIMconfigure.h"
#include "IFEM.h"
#include "Utilities.h"
#include "MatVec.h"
#include "Vec3.h"
#include "DataExporter.h"
#include "HDF5Restart.h"
#include "HDF5Writer.h"
#include "tinyxml.h"
#include <fstream>
#include <memory>

class TimeStep;
class VTF;

//! Data exporters for the planes
extern std::vector<DataExporter*> plane_exporters;
//! Log files for the planes
extern std::vector< std::shared_ptr<std::ostream> > plane_log_files;


/*!
  \brief Driver class for plane-decoupled 3D problems.
*/

template<class PlaneSolver>
class SIMSemi3D : public SIMadmin, public SIMdependency
{
public:
  typedef typename PlaneSolver::SetupProps SetupProps; //!< Convenience type

  //! \brief Enum announcing the dimensionality.
  enum { dimension = 2 };

  //! \brief The constructor initializes the setup properties.
  explicit SIMSemi3D(const SetupProps& props_) :
    startCtx(0), planes(1), procs_per_plane(1), output_plane(-1),
    direction('Z'), props(props_)
  {
    SIMadmin::myHeading = "Plane-decoupled 3D simulation driver";
  }

  //! \brief The destructor deletes the plane-wise sub-step solvers.
  virtual ~SIMSemi3D()
  {
    if (!m_planes.empty())
      delete m_planes.front();

    for (size_t i = 1; i < m_planes.size(); i++)
    {
      // The other planes do not own the integrand, so clear the pointer to it
      // to avoid that the SIMbase destructor tries to delete it again
      m_planes[i]->clearProblem();
      delete m_planes[i];
    }
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp)
  {
    bool ok = true;
    for (size_t i = 0; i < m_planes.size() && ok; i++)
      ok = m_planes[i]->advanceStep(tp);

    return ok;
  }

  //! \brief Returns the order of the BDF scheme.
  int getBDForder() const
  {
    return m_planes.empty() ? 0 : m_planes.front()->getBDForder();
  }

  //! \brief Returns the visualization dump interval.
  int getDumpInterval() const
  {
    return m_planes.empty() ? 1 : m_planes.front()->getDumpInterval();
  }

  //! \brief Dummy method.
  void printFinalNorms(const TimeStep&) {}

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocess()
  {
    for (size_t i=0;i<m_planes.size();++i)
      if (!m_planes[i]->preprocess())
        return false;

    this->grabPlaneNodes();
    return true;
  }

  //! \brief Get FSI nodes for all planes.
  void grabPlaneNodes()
  {
    planeNodes.resize(this->getNoPlanes());
    for (size_t i=0;i<m_planes.size();++i)
      planeNodes[startCtx+i] = m_planes[i]->getNoNodes();

#ifdef HAVE_MPI
    std::vector<int> send(planeNodes);
    MPI_Allreduce(&send[0], &planeNodes[0], planeNodes.size(),
                  MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
#endif
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "Semi3D"; }

  //! \brief Adds fields to a data exporter.
  void registerFields(DataExporter& exporter)
  {
    std::string name = exporter.getName();
    int plane = 1 + startCtx;
    for (size_t i = 0; i < m_planes.size(); i++, plane++)
      if (plane_exporters.size() <= i) {
        DataExporter* exp = new DataExporter(true,exporter.getStride());
        std::stringstream str;
        str << "_plane" << plane;
        HDF5Writer* hdf = new HDF5Writer(name+str.str(),
                                         m_planes[i]->getProcessAdm(),false);
        exp->registerWriter(hdf);
        plane_exporters.push_back(exp);
        m_planes[i]->registerFields(*exp);
      }
      else
        m_planes[i]->registerFields(*plane_exporters[i]);
  }

  //! \brief Initializes for time-dependent simulation.
  bool init(const TimeStep& tp)
  {
    bool ok = true;
    for (size_t i = 0; i < m_planes.size() && ok; i++)
      ok = m_planes[i]->init(tp);

    return ok;
  }

  //! \brief Dummy method (VTF export is not supported).
  bool saveModel(char*,int&,int&) { return true; }

  //! \brief Dumps all registered fields for each plane.
  bool saveStep(const TimeStep& tp, int&)
  {
    for (size_t i = 0; i < plane_exporters.size(); i++)
      plane_exporters[i]->dumpTimeLevel(&tp);

    // Note: VTF export is not supported for the Semi3D simulators.
    // Instead this method is used for the HDF5 export.
    return true;
  }

  //! \brief Dummy method, no serialization support.
  bool serialize(HDF5Restart::SerializeData&) { return false; }
  //! \brief Dummy method, no deserialization support.
  bool deSerialize(const HDF5Restart::SerializeData&) { return false; }

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  bool solveStep(TimeStep& tp)
  {
    bool ok = true;
    for (size_t i = 0; i < m_planes.size() && ok; i++) {
      m_planes[i]->getProcessAdm().cout <<"\n  Plane = "<< startCtx+i+1 <<":";
      ok = m_planes[i]->solveStep(tp);
    }

    return ok;
  }

  //! \brief Sets the initial conditions.
  bool setInitialConditions()
  {
    bool ok = true;
    for (size_t i=0;i<m_planes.size();++i)
      ok &= m_planes[i]->setInitialConditions();

    return ok;
  }

  //! \brief Initialize the FEM system.
  void initSystem()
  {
    for (size_t i=0;i<m_planes.size();++i)
      m_planes[i]->initSystem();
  }

  //! \brief Reads model data from the specified input file.
  //! \param[in] fileName Name of input file to read data from
  virtual bool read(const char* fileName)
  {
    if (!this->SIMadmin::read(fileName))
      return false;

    // Setup our communicator
#ifdef HAVE_MPI
    size_t loc_planes = planes/(nProc/procs_per_plane);
    MPI_Comm comm;
    MPI_Comm_split(PETSC_COMM_WORLD,
                   myPid/procs_per_plane,
                   myPid%procs_per_plane, &comm);
    startCtx = loc_planes*(myPid/procs_per_plane);
#else
    size_t loc_planes = planes;
#endif

    for (size_t i=0;i<loc_planes;++i) {
      m_planes.push_back(new PlaneSolver(props));
#ifdef HAVE_MPI
      m_planes.back()->setCommunicator(&comm);
#endif
      if (!log_files.empty()) {
        int pid = 0;
#ifdef HAVE_MPI
        MPI_Comm_rank(comm, &pid);
#endif
        if (plane_log_files.size() < loc_planes) {
          std::stringstream str;
          str << log_files <<"_plane"<< startCtx+i+1 <<".log";
          std::shared_ptr<std::ostream> file(new std::ofstream(str.str()));
          plane_log_files.push_back(file);
        }
        m_planes[i]->getProcessAdm().cout.addExtraLog(plane_log_files[i],true);
        m_planes[i]->getProcessAdm().cout.setPIDs(0, pid);
      }
      if (output_plane != -1 && output_plane != (int)(i+startCtx+1))
        m_planes[i]->getProcessAdm().cout.setNull();
      else
        m_planes[i]->getProcessAdm().cout.setStream(std::cout);
    }

    return true;
  }

  using SIMadmin::parse;
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"postprocessing")) {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        opt.parseOutputTag(child);
    }
    if (strcasecmp(elem->Value(),"semi3d"))
      return true;

    utl::getAttribute(elem, "nplanes", planes);
#ifdef HAVE_MPI
    utl::getAttribute(elem, "procs_per_plane", procs_per_plane);
    if (procs_per_plane > (size_t)nProc) procs_per_plane = nProc;
#endif
    std::string dir("Z");
    if (utl::getAttribute(elem, "direction", dir))
      direction = toupper(dir[0]);

    utl::getAttribute(elem,"output_prefix", log_files);
    utl::getAttribute(elem,"output_plane", output_plane);

    IFEM::cout <<"\tSemi3D: "<< direction
               <<" "<< planes <<" planes, "<< procs_per_plane
               <<" proces"<< (procs_per_plane > 1 ? "ses":"s") <<" per plane."
               <<"\n\tSemi3D: Printing output from ";
    if (output_plane == -1)
      IFEM::cout <<"all planes to screen."<< std::endl;
    else
      IFEM::cout <<"plane "<< output_plane <<" to screen."<< std::endl;
    if (!log_files.empty())
      IFEM::cout <<"\tSemi3D: Logging to files with prefix "
                 << log_files <<"."<< std::endl;

    return true;
  }

  using SIMdependency::registerDependency;
  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] patches The geometry the field is defined over
  //! \param[in] diffBasis Different basis for the SIM class and the field
  //! \param[in] component Component to use from field
  template<class T>
  void registerDependency(SIMSemi3D<T>* sim, const std::string& name,
                          short int nvc, const PatchVec& patches,
                          char diffBasis = 0, int component = 0)
  {
    for (size_t i=0;i<m_planes.size(); ++i)
      if (diffBasis)
        m_planes[i]->registerDependency(sim->getPlane(i), name, nvc,
                                        sim->getPlane(i)->getFEModel(),
                                        diffBasis, component);
      else
        m_planes[i]->registerDependency(sim->getPlane(i), name, nvc);
  }

  //! \brief Returns the spatial dimension of plane solvers.
  virtual size_t getNoSpaceDim() const { return 2; }
  //! \brief Returns the number of plane solvers.
  size_t getNoPlanes() const { return planes; }

  //! \brief Returns a pointer to a given plane solver.
  PlaneSolver* getPlane(size_t i) { return m_planes[i]; }
  //! \brief Returns a const reference to the plane solvers.
  const std::vector<PlaneSolver*>& getPlanes() const { return m_planes; }

  //! \brief Returns the context of the first plane on this process.
  size_t getStartContext() const { return startCtx; }

  //! \brief Returns the number or processes participating in each planar solve.
  size_t getProcsPerPlane() const { return procs_per_plane; }

  //! \brief Returns the number of FSI nodes for a plane.
  int getPlaneNodes(int plane) const { return planeNodes[plane]; }

  //! \brief Returns the (unoriented) normal direction of the planes.
  char getDirection() const { return direction; }

  //! \brief Dummy method.
  VTF* getVTF() { return NULL; }
  //! \brief Dummy method.
  void setVTF(VTF*) {}

  //! \brief Dummy method.
  Vec3 getForce(int) const { return Vec3(); }
  //! \brief Dummy method.
  const PatchVec& getFEModel() { static PatchVec vec; return vec; }

  //! \brief Setup inter-SIM dependencies.
  void setupDependencies()
  {
    for (size_t i=0;i<m_planes.size();++i)
      m_planes[i]->setupDependencies();
  }

  //! \brief Updating the grid in an ALE solver.
  bool updateALE()
  {
    bool ok = true;
    for (size_t i = 0; i < m_planes.size() && ok; i++)
      ok = m_planes[i]->updateALE();

    return ok;
  }

  //! \brief Dummy method.
  int getLocalNode(int) const { return -1; }
  //! \brief Dummy method.
  int getGlobalNode(int) const { return -1; }

  //! \brief Returns the number of bases in the model.
  int getNoBasis() const { return m_planes.front()->getNoBasis(); }

  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const { return m_planes.front()->getNoSolutions(); }

protected:
  std::vector<PlaneSolver*> m_planes; //!< Planar solvers

private:
  size_t startCtx;             //!< Context for first plane on this process
  size_t planes;               //!< Total number of planes
  size_t procs_per_plane;      //!< Number of processes per plane
  int    output_plane;         //!< Plane to print to screen for (-1 for all)
  char   direction;            //!< (Unoriented) normal direction of plane
  std::string log_files;       //!< Log file prefix for planes
  std::vector<int> planeNodes; //!< FSI nodes for all planes
  SetupProps props;            //!< Setup properties to configure planar solvers
};


//! \brief Configuration template for a SIMSemi3D
template<class PlaneSolver>
struct SolverConfigurator< SIMSemi3D<PlaneSolver> > {
  //! \brief Configure a SIMSemi3D
  //! \param simulator The simulator to configure
  //! \param props The setup properties to use
  //! \param infile The input file to parse
  int setup(SIMSemi3D<PlaneSolver>& simulator,
            const typename PlaneSolver::SetupProps& props,
            char* infile)
  {
    int retval = simulator.read(infile) ? 0 : 1;

    for (size_t i = 0; i < simulator.getPlanes().size() && retval == 0; i++) {
      simulator.getPlane(i)->setContext(i+simulator.getStartContext()+1);
      retval = ConfigureSIM(*simulator.getPlane(i), infile, props);
    }
    simulator.opt = simulator.getPlane(0)->opt;

    return retval;
  }
};

#endif
