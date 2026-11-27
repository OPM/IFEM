// $Id$
//==============================================================================
//!
//! \file TestHDF5Writer.C
//!
//! \date Oct 15 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for the HDF5 result output.
//!
//==============================================================================

#ifdef HAS_HDF5

#include "DataExporter.h"
#include "HDF5Writer.h"
#include "IntegrandBase.h"
#include "ProcessAdm.h"
#include "ASMbase.h"
#include "SIM2D.h"

#include "Catch2Support.h"

#include <hdf5.h>
#include <algorithm>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>


namespace {

class DummyIntegrand : public IntegrandBase
{
public:
  explicit DummyIntegrand(int dim) : IntegrandBase(dim) {}

  using IntegrandBase::evalSol;
  bool evalSol(Vector& s, const FiniteElement&, const Vec3& X,
               const std::vector<int>&) const override
  {
    s.resize(2);
    s(1) = X.sum();
    s(2) = 2.0*X.sum();
    return true;
  }

  size_t getNoFields(int fld) const override { return fld < 2 ? 1 : 2; }
};


class DummySIM : public SIM2D
{
public:
  explicit DummySIM(bool multiPatch, const CharVec& unf = {1,1}) : SIM2D(unf)
  {
    myProblem = new DummyIntegrand(2);
    if (multiPatch) {
      REQUIRE(this->loadXML("<geometry>"
        "<patchfile>src/ASM/Test/refdata/square-4-orient0.g2</patchfile>"
        "<topology>"
        "  <connection master='1' medge='4' slave='2' sedge='3'/>"
        "  <connection master='1' medge='2' slave='3' sedge='1'/>"
        "  <connection master='2' medge='2' slave='4' sedge='1'/>"
        "  <connection master='3' medge='4' slave='4' sedge='3'/>"
        "</topology></geometry>"));
      REQUIRE(this->createFEMmodel());
    }
    else
      REQUIRE(this->createDefaultModel());
    REQUIRE(this->preprocess());
  }

  std::string getName() const override { return "Dummy"; }

  //! \brief Keeps only the given patches, to simulate patch DD.
  void ownOnly(const IntVec& mine, int nproc)
  {
    myPatches = mine;
    nProc = nproc;
  }
};


class ScratchFile
{
public:
  explicit ScratchFile(const char* name)
    : myStem((std::filesystem::temp_directory_path() / name).string()) {}
  ~ScratchFile()
  {
    std::error_code ignored;
    std::filesystem::remove(myStem + ".hdf5", ignored);
  }

  ScratchFile(const ScratchFile&) = delete;
  ScratchFile& operator=(const ScratchFile&) = delete;

  //! \brief Returns the name the writer is to be opened on.
  const std::string& stem() const { return myStem; }

private:
  std::string myStem; //!< Path of the file without its suffix
};


//! \brief What a walk of the file collects.
struct Walk
{
  std::string prefix; //!< Path of the group being walked
  std::vector<std::pair<std::string,hsize_t>>* found; //!< Datasets so far
};


//! \brief Records a dataset, or walks into a group.
herr_t collect (hid_t loc, const char* name, void* data)
{
  Walk* walk = static_cast<Walk*>(data);
  const std::string path = walk->prefix.empty() ? name
                                                : walk->prefix + "/" + name;

  if (hid_t set = H5Dopen2(loc,name,H5P_DEFAULT); set >= 0) {
    hid_t space = H5Dget_space(set);
    walk->found->emplace_back(path, H5Sget_simple_extent_npoints(space));
    H5Sclose(space);
    H5Dclose(set);
  }
  else if (hid_t grp = H5Gopen2(loc,name,H5P_DEFAULT); grp >= 0) {
    Walk below{path, walk->found};
    int idx = 0;
    H5Giterate(grp, ".", &idx, collect, &below);
    H5Gclose(grp);
  }

  return 0;
}


//! \brief The datasets of an HDF5 file, by path, with their lengths.
std::vector<std::pair<std::string,hsize_t>> contents (const std::string& file)
{
  std::vector<std::pair<std::string,hsize_t>> found;
  hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  REQUIRE(f >= 0);

  H5E_auto2_t func = nullptr;
  void* client = nullptr;
  H5Eget_auto2(H5E_DEFAULT,&func,&client);
  H5Eset_auto2(H5E_DEFAULT,nullptr,nullptr);

  Walk walk{"",&found};
  int idx = 0;
  H5Giterate(f, "/", &idx, collect, &walk);

  H5Eset_auto2(H5E_DEFAULT,func,client);
  H5Fclose(f);
  return found;
}


//! \brief Writes one time level of a simulator and returns what it holds.
std::vector<std::pair<std::string,hsize_t>> dump (SIMbase& sim,
                                                  const Vector& sol,
                                                  const std::string& stem)
{
  ProcessAdm adm;
  {
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(stem,adm));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY |
                                   DataExporter::SECONDARY));
    REQUIRE(exporter.setFieldValue("u",&sim,&sol));
    REQUIRE(exporter.dumpTimeLevel());
  }
  return contents(stem + ".hdf5");
}

//! \brief Returns the contents of one dataset.
std::vector<double> values (const std::string& file, const std::string& path)
{
  hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  REQUIRE(f >= 0);
  hid_t set = H5Dopen2(f, path.c_str(), H5P_DEFAULT);
  REQUIRE(set >= 0);
  hid_t space = H5Dget_space(set);
  hsize_t len = 0;
  H5Sget_simple_extent_dims(space,&len,nullptr);
  std::vector<double> read(len);
  if (len > 0)
    H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,read.data());
  H5Sclose(space);
  H5Dclose(set);
  H5Fclose(f);
  return read;
}


//! \brief Writes a field of the projection basis and returns what it holds.
std::vector<std::pair<std::string,hsize_t>> dumpProjected (SIMbase& sim,
                                                           const Vector& sol,
                                                           const std::string& stem,
                                                           int ncmp)
{
  ProcessAdm adm;
  {
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(stem,adm));
    REQUIRE(exporter.registerField("rs","a&&b&&c&&d&&e&&f",DataExporter::SIM,
                                   -(DataExporter::PRIMARY |
                                     DataExporter::PROJECTED),"",ncmp));
    REQUIRE(exporter.setFieldValue("rs",&sim,&sol));
    REQUIRE(exporter.dumpTimeLevel());
  }
  return contents(stem + ".hdf5");
}


//! \brief The projection nodes of every patch, in order.
std::vector<size_t> projNodes (const SIMbase& sim)
{
  std::vector<size_t> n;
  for (int i = 1; i <= sim.getNoPatches(); ++i)
    n.push_back(sim.getPatch(i)->getNoProjectionNodes());
  return n;
}

} // anonymous namespace


TEST_CASE("TestHDF5Writer.OnePatch")
{
  DummySIM sim(false);
  Vector sol(sim.getNoDOFs());
  std::iota(sol.begin(), sol.end(), 1.0);

  ScratchFile out("TestHDF5Writer_one");
  const auto held = dump(sim, sol, out.stem());

  auto has = [&held](const std::string& path)
  {
    return std::find_if(held.begin(), held.end(),
                        [&path](const auto& e) { return e.first == path; })
           != held.end();
  };

  REQUIRE(has("0/Dummy-1/basis/1"));
  REQUIRE(has("0/Dummy-2/basis/1"));

  size_t fields = 0;
  for (const auto& [path,len] : held)
    if (path.find("0/Dummy-") == 0 && path.find("/fields/") != std::string::npos) {
      INFO(path << " is empty");
      REQUIRE(len > 0);
      ++fields;
    }
  REQUIRE(fields >= 3);
}


TEST_CASE("TestHDF5Writer.EveryPatchGetsARecord")
{
  DummySIM sim(true);
  Vector sol(sim.getNoDOFs());
  std::iota(sol.begin(), sol.end(), 1.0);
  sim.ownOnly({1,3}, 2);   // the other two are another process's

  ScratchFile out("TestHDF5Writer_two");
  const auto held = dump(sim, sol, out.stem());

  auto recordsOf = [&held](const std::string& patch)
  {
    std::vector<std::string> names;
    for (const auto& [path,len] : held) {
      const std::string tail = "/" + patch;
      if (path.size() > tail.size() &&
          path.compare(path.size()-tail.size(), tail.size(), tail) == 0)
        names.push_back(path.substr(0,path.size()-tail.size()));
    }
    std::sort(names.begin(), names.end());
    return names;
  };

  const std::vector<std::string> mine = recordsOf("1");
  REQUIRE(!mine.empty());
  for (const char* other : {"2","3","4"}) {
    INFO("patch " << other << " does not carry the same records as patch 1");
    REQUIRE(recordsOf(other) == mine);
  }

  // and what is written for a patch of another process holds nothing
  for (const auto& [path,len] : held)
    if (path.find("/fields/") != std::string::npos &&
        (path.back() == '2' || path.back() == '4')) {
      INFO(path << " should be empty");
      REQUIRE(len == 0);
    }
}


/*!
  A projected field is of no basis of the model's own numbering: it is laid
  out over the projection nodes of each patch in turn. So it is sliced by
  those rather than extracted, it is written under the basis it was fitted
  onto, and six components of it are six components and not a translation
  and a rotation.
*/

TEST_CASE("TestHDF5Writer.AProjectedFieldIsWrittenOnItsOwnBasis")
{
  const int ncmp = 6;

  DummySIM sim(true);
  const std::vector<size_t> nodes = projNodes(sim);
  size_t total = 0;
  for (size_t n : nodes)
    total += ncmp*n;

  Vector sol(total);
  std::iota(sol.begin(), sol.end(), 1.0);

  ScratchFile out("TestHDF5Writer_proj");
  const auto held = dumpProjected(sim, sol, out.stem(), ncmp);

  // written under the basis it was fitted onto, which this model has one of
  REQUIRE(sim.fieldProjections());
  const std::string group = "0/Dummy-proj/fields/";

  // and that basis is in the file, this being the only field of it
  REQUIRE(std::find_if(held.begin(), held.end(), [](const auto& e)
                       { return e.first == "0/Dummy-proj/basis/1"; })
          != held.end());
  size_t ofs = 0;
  for (size_t p = 0; p < nodes.size(); ++p)
  {
    const std::string path = group + "a&&b&&c&&d&&e&&f/" + std::to_string(p+1);
    const auto it = std::find_if(held.begin(), held.end(),
                                 [&path](const auto& e)
                                 { return e.first == path; });
    INFO(path << " is missing");
    REQUIRE(it != held.end());

    // sliced by the projection nodes, and not halved as rotations would be
    INFO(path << " is not as long as the patch is");
    REQUIRE(it->second == ncmp*nodes[p]);

    const std::vector<double> got = values(out.stem() + ".hdf5", path);
    REQUIRE(got.size() == ncmp*nodes[p]);
    REQUIRE_THAT(got.front(), WithinRel(sol[ofs]));
    REQUIRE_THAT(got.back(), WithinRel(sol[ofs + ncmp*nodes[p] - 1]));
    ofs += ncmp*nodes[p];
  }
}


/*!
  Where the model projects onto a basis of its own numbering there is no
  separate one written, and the field belongs under that basis like any
  other rather than under a group with no geometry.
*/

TEST_CASE("TestHDF5Writer.AProjectedFieldOfOneBasisGoesUnderIt")
{
  const int ncmp = 6;

  DummySIM sim(true,{1});
  REQUIRE(!sim.fieldProjections());

  const std::vector<size_t> nodes = projNodes(sim);
  size_t total = 0;
  for (size_t n : nodes)
    total += ncmp*n;

  Vector sol(total);
  std::iota(sol.begin(), sol.end(), 1.0);

  ScratchFile out("TestHDF5Writer_proj1");
  const auto held = dumpProjected(sim, sol, out.stem(), ncmp);

  // the basis it falls back to is in the file, this being its only field
  REQUIRE(std::find_if(held.begin(), held.end(), [](const auto& e)
                       { return e.first == "0/Dummy-1/basis/1"; })
          != held.end());

  size_t ofs = 0;
  for (size_t p = 0; p < nodes.size(); ++p)
  {
    const std::string path = "0/Dummy-1/fields/a&&b&&c&&d&&e&&f/"
                           + std::to_string(p+1);
    const auto it = std::find_if(held.begin(), held.end(),
                                 [&path](const auto& e)
                                 { return e.first == path; });
    INFO(path << " is missing");
    REQUIRE(it != held.end());
    REQUIRE(it->second == ncmp*nodes[p]);

    // laid out over the patches in turn here as well, not over the model
    const std::vector<double> got = values(out.stem() + ".hdf5", path);
    REQUIRE(got.size() == ncmp*nodes[p]);
    REQUIRE_THAT(got.front(), WithinRel(sol[ofs]));
    REQUIRE_THAT(got.back(), WithinRel(sol[ofs + ncmp*nodes[p] - 1]));
    ofs += ncmp*nodes[p];
  }

  for (const auto& [path,len] : held)
  {
    INFO(path << " should not be there, the model projects onto basis one");
    REQUIRE(path.find("Dummy-proj") == std::string::npos);
  }
}


/*!
  A projected field of basis one reaches that basis before an ordinary field
  of it does, the fields being taken in the order of their names. What the
  file holds is the same either way; the test is here for the arrangement,
  which is the one in which the basis was written twice.
*/

TEST_CASE("TestHDF5Writer.AProjectedFieldSharesABasisWithAnOrdinaryOne")
{
  const int ncmp = 6;

  DummySIM sim(true,{1});
  REQUIRE(!sim.fieldProjections());

  const std::vector<size_t> nodes = projNodes(sim);
  size_t total = 0;
  for (size_t n : nodes)
    total += ncmp*n;

  Vector rs(total);
  std::iota(rs.begin(), rs.end(), 1.0);
  Vector sol(sim.getNoDOFs());
  std::iota(sol.begin(), sol.end(), 1.0);

  ScratchFile out("TestHDF5Writer_once");
  {
    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));

    // "rs" is reached before "u" is, the fields being taken in order
    REQUIRE(exporter.registerField("rs","a&&b&&c&&d&&e&&f",DataExporter::SIM,
                                   -(DataExporter::PRIMARY |
                                     DataExporter::PROJECTED),"",ncmp));
    REQUIRE(exporter.setFieldValue("rs",&sim,&rs));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY |
                                   DataExporter::SECONDARY));
    REQUIRE(exporter.setFieldValue("u",&sim,&sol));

    REQUIRE(exporter.dumpTimeLevel());
  }

  const auto held = contents(out.stem() + ".hdf5");
  const auto it = std::find_if(held.begin(), held.end(), [](const auto& e)
                               { return e.first == "0/Dummy-1/basis/1"; });
  REQUIRE(it != held.end());
  REQUIRE(it->second > 0);
}


/*!
  One basis, two fields of it, and the node numbers asked for by the second
  of them. The first writes the geometry; what the second asks for is still
  owed, the geometry and the numbering being two datasets and not one.
*/

TEST_CASE("TestHDF5Writer.ASharedBasisStillOwesTheNodeNumbersAskedForLater")
{
  const int ncmp = 6;

  DummySIM sim(true,{1});
  const std::vector<size_t> nodes = projNodes(sim);
  size_t total = 0;
  for (size_t n : nodes)
    total += ncmp*n;

  Vector rs(total);
  std::iota(rs.begin(), rs.end(), 1.0);
  Vector sol(sim.getNoDOFs());
  std::iota(sol.begin(), sol.end(), 1.0);

  ScratchFile out("TestHDF5Writer_l2g");
  {
    ProcessAdm adm;
    DataExporter exporter(true,1);
    exporter.registerWriter(new HDF5Writer(out.stem(),adm));

    // "rs" is reached first and asks for no numbering; "u" asks for it
    REQUIRE(exporter.registerField("rs","a&&b&&c&&d&&e&&f",DataExporter::SIM,
                                   -(DataExporter::PRIMARY |
                                     DataExporter::PROJECTED),"",ncmp));
    REQUIRE(exporter.setFieldValue("rs",&sim,&rs));
    REQUIRE(exporter.registerField("u","primary",DataExporter::SIM,
                                   DataExporter::PRIMARY |
                                   DataExporter::SECONDARY |
                                   DataExporter::L2G_NODE));
    REQUIRE(exporter.setFieldValue("u",&sim,&sol));
    REQUIRE(exporter.dumpTimeLevel());
  }

  const auto held = contents(out.stem() + ".hdf5");
  auto has = [&held](const std::string& path)
  {
    return std::find_if(held.begin(), held.end(), [&path](const auto& e)
                        { return e.first == path; }) != held.end();
  };

  for (size_t p = 1; p <= nodes.size(); ++p)
  {
    const std::string n = std::to_string(p);
    INFO("the geometry of patch " << n << " is missing");
    REQUIRE(has("0/Dummy-1/basis/" + n));
    INFO("the node numbers of patch " << n << " were never written");
    REQUIRE(has("0/Dummy-1/l2g-node/" + n));
  }
}

#endif // HAS_HDF5
