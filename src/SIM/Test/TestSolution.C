// $Id$
//==============================================================================
//!
//! \file TestSolution.C
//!
//! \date Feb 14 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for the SIMsolution container.
//!
//==============================================================================

#include "SIMsolution.h"

#include <catch2/catch_test_macros.hpp>


class DummySIM : public SIMsolution
{
  int nVecState;

public:
  DummySIM(int nState, int nvs) : nVecState(nvs)
  {
    solution.resize(nState*nVecState,Vector(5));
    for (size_t i = 1; i <= solution.size(); i++)
      solution[i-1].fill(i);
  }

  void push()
  {
    this->pushSolution(nVecState);
    for (int i = 0; i < nVecState; i++)
      solution[i].fill(0.0);
  }

  void printMe()
  {
    int i = 0;
    for (const Vector& v : solution) std::cout <<"s"<< ++i <<":"<< v;
  }
};


TEST_CASE("TestSIMsolution.Push")
{
  DummySIM s1(2,1);
  s1.printMe();
  s1.push();
  s1.printMe();
  REQUIRE(s1.getSolution(0)[0] == 0);
  REQUIRE(s1.getSolution(1)[0] == 1);

  DummySIM s2(4,3);
  s2.printMe();
  s2.push();
  s2.printMe();
  REQUIRE(s2.getSolution(2)[0] == 0);
  REQUIRE(s2.getSolution(5)[0] == 3);
  REQUIRE(s2.getSolution(8)[0] == 6);
  REQUIRE(s2.getSolution(9)[0] == 7);
  s2.push();
  s2.printMe();
  REQUIRE(s2.getSolution(2)[0] == 0);
  REQUIRE(s2.getSolution(5)[0] == 0);
  REQUIRE(s2.getSolution(8)[0] == 3);
  REQUIRE(s2.getSolution(9)[0] == 4);
}
