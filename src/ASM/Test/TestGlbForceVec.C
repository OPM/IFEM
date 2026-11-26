//==============================================================================
//!
//! \file TestGlbForceVec.C
//!
//! \date Nov 26 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for the global nodal force vector.
//!
//==============================================================================

#include "GlbForceVec.h"

#include "ASMbase.h"
#include "ASMmxBase.h"
#include "ElmMats.h"
#include "SAM.h"
#include "SIM2D.h"
#include "Vec3.h"

#include "Catch2Support.h"

#include <numeric>
#include <vector>


namespace {

using CharVec = std::vector<unsigned char>;

//! \brief A simulator read from a file and preprocessed.
class MySIM : public SIM2D
{
public:
  template<class T> MySIM (const T& nf) : SIM2D(nf)
  {
    REQUIRE(this->read("src/ASM/Test/refdata/forcevec_2D_1P.xinp"));
    REQUIRE(this->preprocess());
  }
};

//! \brief Element forces holding one value per degree of freedom.
class MyForces : public ElmMats
{
public:
  //! \brief The constructor takes the right-hand side to hand out.
  explicit MyForces (const Vector& v) : ElmMats(false)
  {
    this->resize(0,1);
    b.front() = v;
  }
};

//! \brief Numbers the degrees of freedom of an element from one.
//! \param[in] sam Assembly management data
//! \param[in] iel 1-based element number
Vector countingForces (const SAM& sam, int iel)
{
  std::vector<int> mnpc;
  REQUIRE(sam.getElmNodes(mnpc,iel));

  size_t ndof = 0;
  for (int inod : mnpc) {
    auto dofs = sam.getNodeDOFs(inod);
    ndof += dofs.second-dofs.first+1;
  }

  Vector v(ndof);
  std::iota(v.begin(),v.end(),1.0);
  return v;
}

}


TEST_CASE("TestGlbForceVec.EveryComponentOnANode")
{
  // Where a node carries every component, its force is its own degrees of
  // freedom in the order they are numbered in.
  ASMmxBase::Type = ASMmxBase::NONE;
  MySIM sim(2);
  const SAM& sam = *sim.getSAM();

  std::vector<int> nodes(sam.getNoNodes());
  std::iota(nodes.begin(),nodes.end(),1);

  GlbForceVec force(sam);
  REQUIRE(force.initNodeMap(nodes,2));
  force.initialize();

  const Vector ES = countingForces(sam,1);
  MyForces elm(ES);
  REQUIRE(force.assemble(&elm,1));

  std::vector<int> mnpc;
  REQUIRE(sam.getElmNodes(mnpc,1));
  for (size_t i = 0; i < mnpc.size(); ++i) {
    const Vec3 f = force.getForce(mnpc[i]);
    CHECK_THAT(f.x, WithinAbs(-ES[2*i],1e-14));
    CHECK_THAT(f.y, WithinAbs(-ES[2*i+1],1e-14));
  }
}


TEST_CASE("TestGlbForceVec.OneComponentToABasis")
{
  // A div-compatible pair keeps one component of the velocity to a basis of
  // its own, so a node gives the one component its single degree of freedom
  // is, and the pressure nodes give none. Reading two values from a node
  // which has one would take the next node's force for this one's.
  ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  ASMmxBase::itgBasis = 3;
  MySIM sim(CharVec{1,1,1});
  const SAM& sam = *sim.getSAM();
  const ASMbase* patch = sim.getPatch(1);

  const size_t nb1 = patch->getNoNodes(1);
  const size_t nb2 = patch->getNoNodes(2);

  std::vector<int> nodes(sam.getNoNodes());
  std::iota(nodes.begin(),nodes.end(),1);

  // The first basis is the x-component of the velocity, the second the y,
  // and the third is the pressure, which is no part of a force
  std::vector<char> firstComp(nodes.size(),3);
  for (size_t i = 0; i < nodes.size(); ++i)
    firstComp[i] = i < nb1 ? 1 : (i < nb1+nb2 ? 2 : 0);

  GlbForceVec force(sam);
  REQUIRE(force.initNodeMap(nodes,2,firstComp));
  force.initialize();

  const Vector ES = countingForces(sam,1);
  MyForces elm(ES);
  REQUIRE(force.assemble(&elm,1));

  std::vector<int> mnpc;
  REQUIRE(sam.getElmNodes(mnpc,1));
  for (size_t i = 0; i < mnpc.size(); ++i) {
    const size_t inod = mnpc[i] - 1;
    const Vec3 f = force.getForce(mnpc[i]);
    CHECK_THAT(f.x, WithinAbs(inod < nb1 ? -ES[i] : 0.0, 1e-14));
    CHECK_THAT(f.y, WithinAbs(inod >= nb1 && inod < nb1+nb2 ? -ES[i] : 0.0,
                              1e-14));
    CHECK_THAT(f.z, WithinAbs(0.0,1e-14));
  }
}
