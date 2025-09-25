//==============================================================================
//!
//! \file TestControlFIFO.C
//!
//! \date Oct 10 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for application control over a FIFO.
//!
//==============================================================================

#include "ControlFIFO.h"

#include <catch2/catch_test_macros.hpp>

#include "tinyxml2.h"
#include <fcntl.h>
#include <unistd.h>


class MockCallback : public ControlCallback
{
  public:
  MockCallback() : callback1(false), callback2(false) {}

  void OnControl(const tinyxml2::XMLElement* context)
  {
    if (!context)
      return;

    const tinyxml2::XMLElement* elem2 = context->FirstChildElement();
    if (elem2 && strcmp(elem2->Value(), "callback1") == 0)
      callback1 = true;
    if (elem2 && strcmp(elem2->Value(), "callback2") == 0)
      callback2 = true;
  }
  std::string GetContext() const { return "test"; }

  bool callback1;
  bool callback2;
};

TEST_CASE("TestControlFIFO.General")
{
  ControlFIFO fifo;
  MockCallback callback;
  fifo.open("/tmp/ifem-control");
  fifo.registerCallback(callback);

  int f = ::open("/tmp/ifem-control", O_WRONLY);
  std::string data;
  data = "<callback><test><callback1/></test></callback>";
  if (write(f, data.c_str(), data.size() + 1) != static_cast<int>(data.size() + 1))
    REQUIRE(false);

  fifo.poll();
  data = "<callback><test><callback2/></test></callback>";
  if (write(f, data.c_str(), data.size() + 1) != static_cast<int>(data.size() + 1))
    REQUIRE(false);
  fifo.poll();
  close(f);

  REQUIRE(callback.callback1);
  REQUIRE(callback.callback2);
}
