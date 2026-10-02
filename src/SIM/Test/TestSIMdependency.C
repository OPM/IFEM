//==============================================================================
//!
//! \file TestSIMdependency.C
//!
//! \date Oct 2 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for dependent fields between simulators.
//!
//==============================================================================

#include "SIM2D.h"
#include "Field.h"
#include "IntegrandBase.h"

#include "Catch2Support.h"

#include <memory>


namespace {

//! \brief An integrand receiving two scalar dependent fields.
class DependentIntegrand : public IntegrandBase
{
public:
  DependentIntegrand() : IntegrandBase(2)
  {
    registerVector("nodal",&nodal);
    registerVector("solution",&solution);
  }

  void setNamedField(const std::string& name, Field* field) override
  {
    if (name == "nodal")
      nodalField.reset(field);
    else if (name == "solution")
      solutionField.reset(field);
    else
      delete field;
  }

  Vector nodal;    //!< Patch-level nodal values
  Vector solution; //!< Patch-level solution values
  std::unique_ptr<Field> nodalField;    //!< Field of the nodal values
  std::unique_ptr<Field> solutionField; //!< Field of a solution component
};


//! \brief A simulator exposing the extraction of its dependent fields.
class DependentSIM : public SIM2D
{
public:
  explicit DependentSIM(unsigned char nf) : SIM2D(nf) {}
  using SIMdependency::extractPatchDependencies;
};


//! \brief Loads a refined unit square into a simulator.
void loadSquare(SIM2D& sim)
{
  REQUIRE(sim.loadXML("<geometry><refine patch='1' u='2' v='3'/></geometry>"));
  REQUIRE(sim.preprocess());
}

}


TEST_CASE("TestSIMdependency.ScalarFieldOnMultiFieldPatch")
{
  // The fields come from a simulator with two fields per node
  SIM2D source(2);
  loadSquare(source);
  const size_t nNodes = source.getNoNodes();

  // A nodal field with one value per node, and a solution vector
  Vector nodal(nNodes), solution(2*nNodes);
  for (size_t i = 1; i <= nNodes; i++) {
    nodal(i) = i;
    solution(2*i-1) = 10.0*i;
    solution(2*i) = 100.0*i;
  }
  source.registerField("nodal",nodal);
  source.registerField("solution",solution);

  DependentSIM target(1);
  loadSquare(target);
  DependentIntegrand problem;
  target.registerDependency(&source,"nodal",1,source.getFEModel(),1);
  target.registerDependency(&source,"solution",1,source.getFEModel(),1,2);
  REQUIRE(target.extractPatchDependencies(&problem,target.getFEModel(),0));

  // The nodal field is extracted with one value per node,
  // and the second component of the solution
  REQUIRE(problem.nodal.size() == nNodes);
  REQUIRE(problem.nodalField);
  REQUIRE(problem.solutionField);
  for (size_t i = 1; i <= nNodes; i++) {
    REQUIRE_THAT(problem.nodalField->valueNode(i), WithinRel(double(i)));
    REQUIRE_THAT(problem.solutionField->valueNode(i), WithinRel(100.0*i));
  }
}
