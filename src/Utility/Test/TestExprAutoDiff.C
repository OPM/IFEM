//==============================================================================
//!
//! \file TestExprAutoDiff.C
//!
//! \date Jun 6 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for autodifferentiation of expression functions.
//!
//==============================================================================

#include "ExprFunctions.h"

#include <autodiff/reverse/var.hpp>
#include <expreval.h>

#include <iostream>

#include "gtest/gtest.h"

struct ExprAutoDiffImpl {
  using Type = autodiff::var;
  ExprEval::Expression<Type> expr;
  ExprEval::FunctionList<Type> flist;
  ExprEval::ValueList<Type> vlist;
  Type* xval;
  Type* yval;

  ExprAutoDiffImpl(const std::string& func,
                   Type xVal, Type yVal)
  {
    flist.AddDefaultFunctions();

    vlist.Add("x", 0.0, false);
    vlist.Add("y", 0.0, false);

    expr.SetFunctionList(&flist);
    expr.SetValueList(&vlist);
    expr.Parse(func);

    xval = vlist.GetAddress("x");
    yval = vlist.GetAddress("y");
    *xval = xVal;
    *yval = yVal;
  }
};

struct TestDef {
  std::string name;
  std::string func;
  double x;
  double y;
  std::array<double, 6> vals;
};

class TestExprAutoDiff : public testing::Test,
                         public testing::WithParamInterface<TestDef>
{
};


TEST_P(TestExprAutoDiff, Derivatives)
{
  ExprAutoDiffImpl e(GetParam().func, GetParam().x, GetParam().y);

  auto val = e.expr.Evaluate();
  auto [ux, uy] = derivativesx(val, wrt(*e.xval, *e.yval));
  auto [uxx, uxy] = derivativesx(ux, wrt(*e.xval, *e.yval));
  auto [uyx, uyy] = derivativesx(uy, wrt(*e.xval, *e.yval));
  EXPECT_DOUBLE_EQ(val.expr->val, GetParam().vals[0]);
  EXPECT_DOUBLE_EQ(ux.expr->val, GetParam().vals[1]);
  EXPECT_DOUBLE_EQ(uy.expr->val, GetParam().vals[2]);
  EXPECT_DOUBLE_EQ(uxx.expr->val, GetParam().vals[3]);
  EXPECT_DOUBLE_EQ(uxy.expr->val, GetParam().vals[4]);
  EXPECT_DOUBLE_EQ(uyx.expr->val, GetParam().vals[4]);
  EXPECT_DOUBLE_EQ(uyy.expr->val, GetParam().vals[5]);
}

INSTANTIATE_TEST_SUITE_P(TestExprAutoDiff, TestExprAutoDiff,
                         testing::Values(
                             TestDef{"Abs", "abs(x) + abs(y)", 1.0, -2.0,
                                     {1.0 + 2.0,
                                      1.0,
                                      -2.0 / 2.0,
                                      0.0,
                                      0.0,
                                      0.0}},
                             TestDef{"ACos", "acos(rad(x)) + acos(y)", 45.0, M_PI / 4.0,
                                     {acos(45.0 * M_PI / 180.0) + acos(M_PI / 4.0),
                                      -M_PI / sqrt(32400 - M_PI * M_PI * 45 * 45),
                                      -1.0 / sqrt(1.0 - M_PI / 4.0 * M_PI / 4.0),
                                      -pow(M_PI, 3.0) * 45.0 / pow(32400 - M_PI*M_PI*45*45, 3.0 / 2.0),
                                      0.0,
                                      -M_PI / (4.0 * pow(1.0 - M_PI / 4.0 * M_PI / 4.0, 3.0 / 2.0))}},
                             TestDef{"ASin", "asin(rad(x)) + asin(y)", 45.0, M_PI / 4.0,
                                     {asin(45.0 * M_PI / 180.0) + asin(M_PI / 4.0),
                                      M_PI / sqrt(32400 - M_PI * M_PI * 45 * 45),
                                      1.0 / sqrt(1.0 - M_PI / 4.0 * M_PI / 4.0),
                                      pow(M_PI, 3.0) * 45.0 / pow(32400 - M_PI * M_PI * 45.0 * 45.0, 3.0 / 2.0),
                                      0.0,
                                      M_PI / (4.0 * pow(1.0 - M_PI * M_PI / 4.0 / 4.0, 3.0 / 2.0))}},
                             TestDef{"ATan", "atan(2*x) * atan(y)", 1.0, 2.0,
                                     {atan(2.0 * 1.0) * atan(2.0),
                                      2.0 / (4.0 * 1.0 * 1.0 + 1.0) * atan(2.0),
                                      atan(2.0 * 1.0) * 1.0 / (2.0 * 2.0 + 1.0),
                                      -16.0 * 1.0 / pow(4.0 * 1.0 * 1.0 + 1.0, 2.0) * atan(2.0),
                                      2.0 / (4.0 * 1.0 * 1.0 + 1.0) * 1.0 / (2.0 * 2.0 + 1.0),
                                      atan(2.0 * 1.0) * -2.0 * 2.0 / pow(2.0 * 2.0 + 1.0, 2.0)}},
                             TestDef{"ATan2", "atan2(x,y)", 1.0, 2.0,
                                     {atan2(1.0, 2.0),
                                      2.0 / (1.0 * 1.0 + 2.0 * 2.0),
                                      -1.0 / (1.0 * 1.0 + 2.0 * 2.0),
                                      -2.0 * 1.0 * 2.0 / pow(1.0 * 1.0 + 2.0 * 2.0, 2.0),
                                      (1.0 * 1.0 - 2.0 * 2.0) / pow(1.0 * 1.0 + 2.0 * 2.0, 2.0),
                                      2.0 * 1.0 * 2.0 / pow(1.0 * 1.0 + 2.0 * 2.0, 2.0)}},
                             TestDef{"Cos", "cos(rad(x)) + cos(y)", 45.0, M_PI / 4.0,
                                     {cos(45 * M_PI / 180.0) + cos(M_PI / 4.0),
                                      -M_PI / 180.0 * sin(45 * M_PI / 180.0),
                                      -sin(M_PI / 4.0),
                                      -M_PI * M_PI / (180.0 * 180.0) * cos(45 * M_PI / 180.0),
                                      0.0,
                                      -cos(M_PI / 4.0)}},
                             TestDef{"Cosh", "cosh(2*x) * cosh(y)", 1.0, 2.0,
                                     {cosh(2.0 * 1.0) * cosh(2.0),
                                      2.0 * sinh(2.0 * 1.0) * cosh(2.0),
                                      cosh(2.0 * 1.0) * sinh(2.0),
                                      2.0 * 2.0 * cosh(2.0 * 1.0) * cosh(2.0),
                                      2.0 * sinh(2.0 * 1.0) * sinh(2.0),
                                      cosh(2.0 * 1.0) * cosh(2.0)}},
                             TestDef{"Deg", "deg(x*x)", 1.0, 2.0,
                                     {1.0 * 1.0 / M_PI * 180.0,
                                      360.0 * 1.0 / M_PI,
                                      0.0,
                                      360.0 / M_PI,
                                      0.0,
                                      0.0}},
                             TestDef{"Exp", "exp(2*x) * exp(y)", 1.0, 2.0,
                                     {exp(2.0 * 1.0) * exp(2.0),
                                      2.0 * exp(2.0 * 1.0) * exp(2.0),
                                      exp(2.0 * 1.0) * exp(2.0),
                                      4.0 * exp(2.0 * 1.0) * exp(2.0),
                                      2.0 * exp(2.0 * 1.0) * exp(2.0),
                                      exp(2.0 * 1.0) * exp(2.0)}},
                             TestDef{"Ln", "ln(2*x) * ln(y)", 1.0, 2.0,
                                     {log(2.0 * 1.0) * log(2.0),
                                      1.0 / 1.0 * log(2.0),
                                      log(2.0 * 1.0) / 2.0,
                                      -1.0 / (1.0 * 1.0) * log(2.0),
                                      1.0 / (1.0) * 1.0 / (2.0),
                                      log(2.0 * 1.0) * -1.0 / (2.0 * 2.0)}},
                             TestDef{"Log", "log(2*x) * log(y)", 1.0, 2.0,
                                     {log10(2.0 * 1.0) * log10(2.0),
                                      1.0 / (1.0 * log(10)) * log10(2.0),
                                      log10(2.0 * 1.0) / (2.0 * log(10)),
                                      -1.0 / (1.0 * 1.0 * log(10)) * log10(2.0),
                                      1.0 / (1.0 * log(10)) * 1.0 / (2.0 * log(10)),
                                      log10(2.0 * 1.0) * -1.0 / (2.0 * 2.0 * log(10))}},
                             TestDef{"Logn", "logn(2*x, 5)", 1.0, 2.0,
                                     {log(2.0 * 1.0) / log(5.0),
                                      1.0 / (1.0 * log(5.0)),
                                      0.0,
                                      -1.0 / (1.0 * 1.0 * log(5.0)),
                                      0.0,
                                      0.0}},
                             TestDef{"Pow", "pow(x,2) * pow(y,2)", 1.0, 2.0,
                                     {pow(1.0, 2.0) * pow(2.0, 2.0),
                                      2.0 * 1.0 * pow(2.0, 2.0),
                                      pow(1.0, 2.0) * 2.0 * 2.0,
                                      2.0 * pow(2.0, 2.0),
                                      2.0 * 1.0 * 2.0 * 2.0,
                                      pow(1.0, 2.0) * 2.0}},
                             TestDef{"Rad", "rad(x*x)", 1.0, 2.0,
                                     {1.0 * 1.0 * M_PI / 180.0,
                                      M_PI * 1.0 / 90.0,
                                      0.0,
                                      M_PI / 90.0,
                                      0.0,
                                      0.0}},
                             TestDef{"Sin", "sin(rad(x)) + sin(y)", 45.0, M_PI / 4.0,
                                     {sin(45 * M_PI / 180.0) + sin(M_PI / 4.0),
                                      M_PI / 180.0 * cos(45 * M_PI / 180.0),
                                      cos(M_PI / 4.0),
                                      -M_PI * M_PI / (180.0 * 180.0) * sin(45 * M_PI / 180.0),
                                      0.0,
                                      -sin(M_PI / 4.0)}},
                             TestDef{"Sqrt", "sqrt(x) * sqrt(y)", 1.0, 2.0,
                                     {sqrt(1.0) * sqrt(2.0),
                                      1.0 / (2.0 * sqrt(1.0)) * sqrt(2.0),
                                      sqrt(1.0) / (2.0 * sqrt(2.0)),
                                      -1.0 / 4.0 * sqrt(2.0),
                                      1.0 / (2.0 * sqrt(1.0)) * 1.0 / (2.0*sqrt(2.0)),
                                      sqrt(1.0) * -1.0 / (4.0*pow(2.0, 3.0 / 2.0))}},
                             TestDef{"Sinh", "sinh(2*x) * sinh(y)", 1.0, 2.0,
                                     {sinh(2.0 * 1.0) * sinh(2.0),
                                      2.0*cosh(2.0 * 1.0) * sinh(2.0),
                                      sinh(2.0 * 1.0) * cosh(2.0),
                                      2.0 * 2.0 * sinh(2.0 * 1.0) * sinh(2.0),
                                      2.0 * cosh(2.0 * 1.0) * cosh(2.0),
                                      sinh(2.0 * 1.0) * sinh(2.0)}},
                             TestDef{"Tan", "tan(rad(x)) + tan(y)", 45.0, M_PI / 4.0,
                                     {tan(45 * M_PI / 180.0) + tan(M_PI / 4.0),
                                      M_PI / (180.0 * pow(cos(45 * M_PI / 180.0), 2.0)),
                                      1.0 / pow (cos(M_PI / 4.0), 2.0),
                                      M_PI * M_PI / (180.0 * 180.0) * 2.0 * sin(45 * M_PI / 180.0) / pow(cos(45 * M_PI / 180.0), 3.0),
                                      0.0,
                                      2.0 * sin(M_PI / 4.0) / pow(cos(M_PI / 4.0), 3.0)}},
                             TestDef{"Tanh", "tanh(2.0*x) * tanh(y)", 1.0, 2.0,
                                     {tanh(2.0 * 1.0) * tanh(2.0),
                                      2.0 / pow(cosh(2.0 * 1.0), 2.0) * tanh(2.0),
                                      tanh(2.0 * 1.0) / pow(cosh(2.0), 2.0),
                                      -2.0 * 2.0 * 2.0 * sinh(2.0 * 1.0) / pow(cosh(2.0 * 1.0), 3.0) * tanh(2.0),
                                      2.0 / pow(cosh(2.0 * 1.0), 2.0) * 1.0 / pow(cosh(2.0), 2.0),
                                      tanh(2.0 * 1.0) * -2.0 * sinh(2.0) / pow(cosh(2.0), 3.0)}}
                         ), [](const testing::TestParamInfo<TestDef>& info)
                            {
                                return info.param.name;
                            });
