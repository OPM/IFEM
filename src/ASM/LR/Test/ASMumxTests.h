//==============================================================================
//!
//! \file ASMumxTests.h
//!
//! \date Nov 10 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for mixed LR spline FE models.
//!
//==============================================================================

#include "ASMbase.h"
#include "ASMmxBase.h"
#include "ASMunstruct.h"
#include "SAMpatch.h"

#include "Catch2Support.h"


namespace {

//! \brief Returns the value field \a c of basis \a b is filled with.
double value(double base, int b, size_t c)
{
  return base + 10.0*b + c;
}

//! \brief Returns the length of a solution vector over all bases.
size_t length(const ASMbase& pch)
{
  size_t len = 0;
  for (int b = 1; b <= static_cast<int>(pch.getNoBasis()); b++)
    len += pch.getNoFields(b)*pch.getNoNodes(b);
  return len;
}

//! \brief Fills a solution vector with one constant per field.
void fill(const ASMbase& pch, Vector& v, double base)
{
  v.resize(length(pch));
  size_t ofs = 0;
  for (int b = 1; b <= static_cast<int>(pch.getNoBasis()); b++) {
    const size_t nf = pch.getNoFields(b);
    const size_t nn = pch.getNoNodes(b);
    for (size_t n = 0; n < nn; n++)
      for (size_t c = 0; c < nf; c++)
        v[ofs+n*nf+c] = value(base,b,c);
    ofs += nf*nn;
  }
}

}


template<class Patch>
class ASMumxTests
{
  //! \brief Checks that a solution vector still holds the constants given.
  static void check(const ASMbase& pch, const Vector& v, double base)
  {
    REQUIRE(v.size() == length(pch));
    size_t ofs = 0;
    for (int b = 1; b <= static_cast<int>(pch.getNoBasis()); b++) {
      const size_t nf = pch.getNoFields(b);
      const size_t nn = pch.getNoNodes(b);
      for (size_t n = 0; n < nn; n++)
        for (size_t c = 0; c < nf; c++)
          REQUIRE_THAT(v[ofs+n*nf+c], WithinRel(value(base,b,c)));
      ofs += nf*nn;
    }
  }

public:
  //! \brief Tests that a refinement carries every solution vector along.
  static void SolutionTransfer()
  {
    ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
    ASMbase::resetNumbering();
    Patch pch({2,1}); // Taylor-Hood, so the two bases differ in field count
    for (unsigned char d = 0; d < pch.getNoParamDim(); d++)
      REQUIRE(pch.uniformRefine(d,1));
    REQUIRE(pch.generateFEMTopology());

    // The solution vectors are taken apart over MADOF, which the assembly
    // sets up, so the patch needs a SAM before it can carry one along
    SAMpatch sam;
    REQUIRE(sam.init({&pch},0,{}));

    // Two vectors, since the offset of the second has to start afresh rather
    // than carry on from where the first one ended
    Vectors sol(2);
    for (size_t i = 0; i < sol.size(); i++)
      fill(pch,sol[i],100.0*(i+1));

    LR::RefineData prm;
    prm.elements = { 0 };
    REQUIRE(pch.refine(prm,sol));

    // A field that is constant over the coarse basis is the same constant over
    // the refined one, so every value carried over is known
    for (size_t i = 0; i < sol.size(); i++)
      check(pch,sol[i],100.0*(i+1));
  }
};
