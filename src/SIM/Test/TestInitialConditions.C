//==============================================================================
//!
//! \file TestInitialConditions.C
//!
//! \date Apr 6 2016
//!
//! \author Timo van Opstal / NTNU
//!
//! \brief Tests initial condition handling.
//!
//==============================================================================

#include "SIM2D.h"
#include "Function.h"

#include "Catch2Support.h"


class TestSIM : public SIM2D
{
public:
  TestSIM() : SIM2D(4) {}
  const InitialCondVec& getIC() const { return myICs.begin()->second; }
};


TEST_CASE("TestInitialConditions.Parse")
{
  const char* geo = "<geometry><topologysets>"
    "<set name='left' type='edge'><item patch='1'>1</item></set>"
    "<set name='right' type='edge'><item patch='1'>2</item></set>"
    "<set name='bottom' type='edge'><item patch='1'>3</item></set>"
    "</topologysets></geometry>";

  const char* dbc = "<boundaryconditions>"
    "<dirichlet set='left' comp='1'>1</dirichlet>"
    "<dirichlet set='right' component='2'>2</dirichlet>"
    "<dirichlet set='bottom' comp='3' component='-3'>3</dirichlet>"
    "<dirichlet set='bottom' component='-4' comp='4'>4</dirichlet>"
    "</boundaryconditions>";

  TestSIM sim;
  REQUIRE(sim.loadXML(geo));
  REQUIRE(sim.loadXML(dbc));
  REQUIRE(sim.loadXML("<initialcondition field='solution' type='expression' comp='1'>1</initialcondition>"));
  REQUIRE(sim.loadXML("<initialcondition field='solution' type='expression' component='2'>2</initialcondition>"));
  REQUIRE(sim.loadXML("<initialcondition field='solution' type='expression' comp='3' component='-3'>3</initialcondition>"));
  REQUIRE(sim.loadXML("<initialcondition field='solution' type='expression' component='-4' comp='4'>4</initialcondition>"));

  // Recognize both comp and component attributes and correct priority
  // Boundary conditions
  for (int i = 1; i <= 4; i++)
    CHECK(static_cast<double>(i) == (*sim.getSclFunc(i))(Vec3()));

  // Initial conditions
  for (const SIMinput::ICInfo& info : sim.getIC())
    switch (info.component)
      {
      case 1:
        CHECK(info.function == "1");
        break;
      case 2:
        CHECK(info.function == "2");
        break;
      case 3:
        CHECK(info.function == "3");
        break;
      case 4:
        CHECK(info.function == "4");
        break;
      default:
        std::cerr <<" *** ICInfo: Invalid component, "<< info.component << std::endl;
        CHECK(false);
      }
}


#ifdef HAS_HDF5

#include "ASMbase.h"
#include "ASMmxBase.h"
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "IntegrandBase.h"
#include "ProcessAdm.h"

#include <filesystem>

namespace {

//! \brief Carries nothing but the field count the writer asks it for.
class FieldIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor records how many fields a node holds.
  FieldIntegrand(int dim, unsigned char nf) : IntegrandBase(dim), myNF(nf) {}

  //! \brief Returns the number of primary fields, on the one basis there is.
  size_t getNoFields(int fld) const override { return fld < 2 ? myNF : 0; }

  //! \brief Names the field, which is what the condition asks for it by.
  std::string getField1Name(size_t, const char* prefix) const override
  {
    return prefix ? std::string(prefix) + " u" : "u";
  }

private:
  unsigned char myNF; //!< Number of fields on the basis
};


//! \brief A simulator over the unit square carrying one basis of \a nf fields.
class FieldSIM : public SIM2D
{
public:
  //! \brief Builds the model and registers the solution under "solution".
  //! \param[in] nf Number of fields on the basis
  //! \param[in] name What the writer names the basis after
  //! \param[in] nLag Number of global Lagrange multipliers to carry
  FieldSIM(unsigned char nf, const char* name, int nLag = 0)
    : SIM2D(nf), myName(name), myNLag(nLag)
  {
    myProblem = new FieldIntegrand(2,nf);
    REQUIRE(this->createDefaultModel());
    REQUIRE(this->preprocess());
    sol.resize(this->getNoDOFs());
    this->registerField("solution",sol);
  }

  //! \brief Returns the name the writer puts on the basis.
  std::string getName() const override { return myName; }

  Vector sol; //!< The field an initial condition is read onto

protected:
  //! \brief Hangs the multipliers off every patch, as a constraint does.
  bool preprocessBeforeAsmInit(int& nnod) override
  {
    for (int i = 0; i < myNLag; i++)
      for (int p = 1; p <= this->getNoPatches(); p++)
        if (int lp = this->getLocalPatchIndex(p); lp > 0)
          this->getPatch(lp)->addGlobalLagrangeMultipliers({++nnod},1);

    return true;
  }

private:
  std::string myName; //!< What the writer names the basis after
  int         myNLag; //!< Number of global Lagrange multipliers
};


//! \brief A simulator over the unit square carrying two bases of one field.
class MixedFieldSIM : public SIM2D
{
public:
  //! \brief Builds the mixed model and registers the solution.
  MixedFieldSIM() : SIM2D(CharVec{1,1})
  {
    myProblem = new FieldIntegrand(2,1);
    REQUIRE(this->createDefaultModel());
    REQUIRE(this->preprocess());
    sol.resize(this->getNoDOFs());
    this->registerField("solution",sol);
  }

  //! \brief Returns the name the writer puts on the bases.
  std::string getName() const override { return "Mx"; }

  Vector sol; //!< The solution of the mixed model
};


//! \brief Removes the file it names when it goes out of scope.
class ScratchFile
{
public:
  //! \brief The constructor puts the file in the directory for temporaries.
  //! \param[in] name What to call it, without a suffix
  explicit ScratchFile(const char* name)
    : myStem((std::filesystem::temp_directory_path() / name).string()) {}
  //! \brief The destructor removes the file, whether it was written or not.
  ~ScratchFile()
  {
    std::error_code ignored;
    std::filesystem::remove(myStem + ".hdf5", ignored);
  }

  //! \brief No copying, since the copy would remove the original's file.
  ScratchFile(const ScratchFile&) = delete;
  //! \brief No assigning, for the reason there is no copying.
  ScratchFile& operator=(const ScratchFile&) = delete;

  //! \brief Returns the name without its suffix, which is what both ends take.
  const std::string& stem() const { return myStem; }

private:
  std::string myStem; //!< Path of the file without its suffix
};

} // anonymous namespace


TEST_CASE("TestInitialConditions.PlaceFewerComponentsThanTheBasis")
{
  ScratchFile out("TestInitialConditions_place");

  // A stored field of two components, each node holding a pair which says
  // where it came from and which of the two it is
  {
    FieldSIM src(2,"Src");
    for (size_t n = 0; 2*n < src.sol.size(); n++) {
      src.sol[2*n]   = 1.0 + n;
      src.sol[2*n+1] = 100.0 + n;
    }

    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY));
    REQUIRE(exporter.setFieldValue("u",&src,&src.sol));
    REQUIRE(exporter.dumpTimeLevel());
  }

  // Said nothing about where they go, so they go first and the rest of the
  // node keeps what it had
  const std::string head = "<initialcondition field='solution' file='" +
                           out.stem() + "' file_basis='Src-1' file_field='u'"
                           " file_level='0'";
  {
    FieldSIM dst(3,"Dst");
    std::fill(dst.sol.begin(),dst.sol.end(),-1.0);
    REQUIRE(dst.loadXML((head + "/>").c_str()));
    REQUIRE(dst.setInitialConditions());

    REQUIRE(dst.sol.size() == 3*dst.getNoNodes());
    for (size_t n = 0; n < dst.getNoNodes(); n++) {
      REQUIRE_THAT(dst.sol[3*n],   WithinAbs(1.0 + n,   1.0e-13));
      REQUIRE_THAT(dst.sol[3*n+1], WithinAbs(100.0 + n, 1.0e-13));
      REQUIRE_THAT(dst.sol[3*n+2], WithinAbs(-1.0,      1.0e-13));
    }
  }

  // The component the condition names is where the first of them goes
  {
    FieldSIM dst(3,"Dst");
    std::fill(dst.sol.begin(),dst.sol.end(),-1.0);
    REQUIRE(dst.loadXML((head + " comp='2'/>").c_str()));
    REQUIRE(dst.setInitialConditions());

    for (size_t n = 0; n < dst.getNoNodes(); n++) {
      REQUIRE_THAT(dst.sol[3*n],   WithinAbs(-1.0,      1.0e-13));
      REQUIRE_THAT(dst.sol[3*n+1], WithinAbs(1.0 + n,   1.0e-13));
      REQUIRE_THAT(dst.sol[3*n+2], WithinAbs(100.0 + n, 1.0e-13));
    }
  }

  // Two components from the third of three do not fit, and that is an error
  // rather than a write past the end of the node
  {
    FieldSIM dst(3,"Dst");
    REQUIRE(dst.loadXML((head + " comp='3'/>").c_str()));
    REQUIRE(!dst.setInitialConditions());
  }

  // A field which fills the basis exactly leaves no room to be offset, and
  // an offset asked for anyway is refused rather than quietly dropped
  {
    FieldSIM dst(2,"Dst");
    REQUIRE(dst.loadXML((head + " comp='2'/>").c_str()));
    REQUIRE(!dst.setInitialConditions());
  }

  // The same field without an offset is what that basis is for
  {
    FieldSIM dst(2,"Dst");
    std::fill(dst.sol.begin(),dst.sol.end(),-1.0);
    REQUIRE(dst.loadXML((head + "/>").c_str()));
    REQUIRE(dst.setInitialConditions());
    for (size_t n = 0; n < dst.getNoNodes(); n++) {
      REQUIRE_THAT(dst.sol[2*n],   WithinAbs(1.0 + n,   1.0e-13));
      REQUIRE_THAT(dst.sol[2*n+1], WithinAbs(100.0 + n, 1.0e-13));
    }
  }

  // A model of one basis answers for basis zero too, and the stored basis
  // is a basis of its own whatever the condition names, so asking for zero
  // is the same read
  {
    FieldSIM dst(2,"Dst");
    std::fill(dst.sol.begin(),dst.sol.end(),-1.0);
    REQUIRE(dst.loadXML((head + " basis='0'/>").c_str()));
    REQUIRE(dst.setInitialConditions());
    for (size_t n = 0; n < dst.getNoNodes(); n++) {
      REQUIRE_THAT(dst.sol[2*n],   WithinAbs(1.0 + n,   1.0e-13));
      REQUIRE_THAT(dst.sol[2*n+1], WithinAbs(100.0 + n, 1.0e-13));
    }
  }
}


TEST_CASE("TestInitialConditions.ReadOntoABasisCarryingMultipliers")
{
  ScratchFile out("TestInitialConditions_lagrange");

  // Three components to the node, which is what the reader's basis holds
  {
    FieldSIM src(3,"Src");
    for (size_t n = 0; 3*n < src.sol.size(); n++) {
      src.sol[3*n]   = 1.0 + n;
      src.sol[3*n+1] = 100.0 + n;
      src.sol[3*n+2] = 200.0 + n;
    }

    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY));
    REQUIRE(exporter.setFieldValue("u",&src,&src.sol));
    REQUIRE(exporter.dumpTimeLevel());
  }

  // The multipliers are nodes of the patch but carry no field, so counting
  // them in would make the three components read as two
  FieldSIM dst(3,"Dst",2);
  std::fill(dst.sol.begin(),dst.sol.end(),-1.0);
  const size_t nnod = dst.getNoNodes(1);
  REQUIRE(dst.sol.size() > 3*nnod); // the multipliers are in there too

  REQUIRE(dst.loadXML(("<initialcondition field='solution' file='" +
                       out.stem() + "' file_basis='Src-1' file_field='u'"
                       " file_level='0' basis='0'/>").c_str()));
  REQUIRE(dst.setInitialConditions());

  for (size_t n = 0; n < nnod; n++) {
    REQUIRE_THAT(dst.sol[3*n],   WithinAbs(1.0 + n,   1.0e-13));
    REQUIRE_THAT(dst.sol[3*n+1], WithinAbs(100.0 + n, 1.0e-13));
    REQUIRE_THAT(dst.sol[3*n+2], WithinAbs(200.0 + n, 1.0e-13));
  }
}


TEST_CASE("TestInitialConditions.ReadOntoAFieldOfItsOwn")
{
  ScratchFile out("TestInitialConditions_own");

  // One value to the node, which is what a field a solver carries alongside
  // its solution rather than within it holds
  {
    FieldSIM src(1,"Src");
    for (size_t n = 0; n < src.sol.size(); n++)
      src.sol[n] = 1.0 + n;

    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY));
    REQUIRE(exporter.setFieldValue("u",&src,&src.sol));
    REQUIRE(exporter.dumpTimeLevel());
  }

  // The basis of the reader holds two fields, but the vector read onto is one
  // of its own with a single value to the node, and that is its layout
  FieldSIM dst(2,"Dst");
  Vector scalar(dst.getNoNodes(1));
  std::fill(scalar.begin(),scalar.end(),-1.0);
  dst.registerField("scalar",scalar);

  REQUIRE(dst.loadXML(("<initialcondition field='scalar' file='" + out.stem() +
                       "' file_basis='Src-1' file_field='u'"
                       " file_level='0'/>").c_str()));
  REQUIRE(dst.setInitialConditions());

  REQUIRE(scalar.size() == dst.getNoNodes(1));
  for (size_t n = 0; n < scalar.size(); n++)
    REQUIRE_THAT(scalar[n], WithinAbs(1.0 + n, 1.0e-13));

  // A mixed patch takes a vector apart over the degrees of freedom of the
  // whole solution, so a vector of the field's own is refused there rather
  // than read past both ends
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  MixedFieldSIM mx;
  Vector own(mx.getNoNodes(1));
  mx.registerField("own",own);

  REQUIRE(mx.loadXML(("<initialcondition field='own' file='" + out.stem() +
                      "' file_basis='Src-1' file_field='u'"
                      " file_level='0'/>").c_str()));
  REQUIRE(!mx.setInitialConditions());
}


TEST_CASE("TestInitialConditions.ReadOntoAMixedBasisNamingNone")
{
  ScratchFile out("TestInitialConditions_mixed");

  {
    FieldSIM src(1,"Src");
    for (size_t n = 0; n < src.sol.size(); n++)
      src.sol[n] = 1.0 + n;

    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY));
    REQUIRE(exporter.setFieldValue("u",&src,&src.sol));
    REQUIRE(exporter.dumpTimeLevel());
  }

  // Each basis of a mixed patch is packed with its own field count, so a
  // condition naming no basis has no one layout to be read over, and is
  // refused rather than laid out over a count which fits none of them
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  MixedFieldSIM mx;

  REQUIRE(mx.loadXML(("<initialcondition field='solution' file='" +
                      out.stem() + "' file_basis='Src-1' file_field='u'"
                      " file_level='0' basis='0'/>").c_str()));
  REQUIRE(!mx.setInitialConditions());
}


#endif
