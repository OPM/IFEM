//==============================================================================
//!
//! \file TestPythonFunctions.C
//!
//! \date Sep 7 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for python functions.
//!
//==============================================================================

#include "PythonFunctions.h"
#include <pybind11/embed.h>

#include "gtest/gtest.h"


TEST(TestPythonFunctions, ScalarFunc)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();
  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");
  double x = 2.0;
  double scale = 1.0;

  PythonFunc f("scalar", R"({"scale" : 1.0})");
  EXPECT_DOUBLE_EQ(f(x), scale*x);

  PythonFunc f2("scalar", R"({"scale" : 10.0})");
  x = 3.0;
  scale = 10.0;
  EXPECT_DOUBLE_EQ(f2(x), scale*x);
}


TEST(TestPythonFunctions, RealFunc)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonFunction f("real", R"({"scale" : 1.0})");
  Vec3 X;
  const double& x = X[0] = 1.0;
  const double& y = X[1] = 2.0;
  double scale = 1.0;

  EXPECT_DOUBLE_EQ(f(X), scale*(x + 100.0*y));
  scale = 10.0;

  PythonFunction f2("real", R"({"scale" : 10.0})");
  EXPECT_DOUBLE_EQ(f2(X), scale*(x + 100.0*y));
}


TEST(TestPythonFunctions, VecFunc)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonVecFunc f("vec", R"({"scale" : 1.0})");
  Vec4 X;
  const double& x = X[0] = 1.0;
  const double& y = X[1] = 2.0;
  const double& z = X[2] = 3.0;
  const double& t = X.t = 4.0;
  double scale = 1.0;

  Vec3 res = f(X);

  EXPECT_DOUBLE_EQ(res[0], scale*t*x);
  EXPECT_DOUBLE_EQ(res[1], scale*t*y);
  EXPECT_DOUBLE_EQ(res[2], scale*t*z);

  PythonVecFunc f2("vec", R"({"scale" : 10.0})");

  res = f2(X);
  scale = 10.0;

  EXPECT_DOUBLE_EQ(res[0], scale*t*x);
  EXPECT_DOUBLE_EQ(res[1], scale*t*y);
  EXPECT_DOUBLE_EQ(res[2], scale*t*z);
}


TEST(TestPythonFunctions, TensorFunc2D)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonTensorFunc f("tensor", R"({"scale" : 1.0})");
  Vec4 X;
  const double& x = X[0] = 1.0;
  const double& y = X[1] = 2.0;
  const double& t = X.t = 4.0;
  double scale = 1.0;

  Tensor res = f(X);

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);

  PythonTensorFunc f2("tensor", R"({"scale" : 10.0})");

  res = f2(X);
  scale = 10.0;

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);
}


TEST(TestPythonFunctions, TensorFunc3D)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonTensorFunc f("tensor", R"({"scale" : 1.0})");
  Vec4 X;
  const double& x = X[0] = 1.0;
  const double& y = X[1] = 2.0;
  const double& z = X[2] = 3.0;
  const double& t = X.t = 4.0;
  double scale = 1.0;

  Tensor res = f(X);

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(1,3), scale*t*x*z);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);
  EXPECT_DOUBLE_EQ(res(2,3), scale*t*y*z);
  EXPECT_DOUBLE_EQ(res(3,1), scale*t*z*x);
  EXPECT_DOUBLE_EQ(res(3,2), scale*t*z*y);
  EXPECT_DOUBLE_EQ(res(3,3), scale*t*z*z);

  PythonTensorFunc f2("tensor", R"({"scale" : 10.0})");

  res = f2(X);
  scale = 10.0;

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(1,3), scale*t*x*z);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);
  EXPECT_DOUBLE_EQ(res(2,3), scale*t*y*z);
  EXPECT_DOUBLE_EQ(res(3,1), scale*t*z*x);
  EXPECT_DOUBLE_EQ(res(3,2), scale*t*z*y);
  EXPECT_DOUBLE_EQ(res(3,3), scale*t*z*z);
}


TEST(TestPythonFunctions, STensorFunc2D)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonSTensorFunc f(2, "symmtensor", R"({"scale" : 1.0})");
  Vec4 X;
  double& x = X[0] = 1.0;
  double& y = X[1] = 2.0;
  double& t = X.t = 4.0;
  double scale = 1.0;

  SymmTensor res = f(X);

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);

  PythonSTensorFunc f2(2, "symmtensor", R"({"scale" : 10.0, "withtt" : true})", true);

  auto res2 = f2(X);
  scale = 10.0;

  EXPECT_DOUBLE_EQ(res2(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res2(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res2(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res2(2,2), scale*t*y*y);
  EXPECT_DOUBLE_EQ(res2(3,3), scale*t*3.0);
}


TEST(TestPythonFunctions, STensorFunc3D)
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  pybind11::module sys = pybind11::module::import("sys");
  pybind11::list path = sys.attr("path");
  path.append("./src/Utility/Test/refdata");

  PythonSTensorFunc f(3, "symmtensor", R"({"scale" : 1.0})");
  Vec4 X;
  double& x = X[0] = 1.0;
  double& y = X[1] = 2.0;
  double& z = X[2] = 3.0;
  double& t = X.t = 4.0;
  double scale = 1.0;

  SymmTensor res = f(X);

  EXPECT_DOUBLE_EQ(res(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res(1,3), scale*t*x*z);
  EXPECT_DOUBLE_EQ(res(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res(2,2), scale*t*y*y);
  EXPECT_DOUBLE_EQ(res(2,3), scale*t*y*z);
  EXPECT_DOUBLE_EQ(res(3,1), scale*t*z*x);
  EXPECT_DOUBLE_EQ(res(3,2), scale*t*z*y);
  EXPECT_DOUBLE_EQ(res(3,3), scale*t*z*z);

  PythonSTensorFunc f2(3, "symmtensor", R"({"scale" : 10.0})");

  auto res2 = f2(X);
  scale = 10.0;

  EXPECT_DOUBLE_EQ(res2(1,1), scale*t*x*x);
  EXPECT_DOUBLE_EQ(res2(1,2), scale*t*x*y);
  EXPECT_DOUBLE_EQ(res2(1,3), scale*t*x*z);
  EXPECT_DOUBLE_EQ(res2(2,1), scale*t*y*x);
  EXPECT_DOUBLE_EQ(res2(2,2), scale*t*y*y);
  EXPECT_DOUBLE_EQ(res2(2,3), scale*t*y*z);
  EXPECT_DOUBLE_EQ(res2(3,1), scale*t*z*x);
  EXPECT_DOUBLE_EQ(res2(3,2), scale*t*z*y);
  EXPECT_DOUBLE_EQ(res2(3,3), scale*t*z*z);
}
