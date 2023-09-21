//==============================================================================
//!
//! \file TestFunctions.C
//!
//! \date Apr 28 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for parsing of functions.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "Functions.h"
#include <cstdlib>
#include <cmath>

#include "gtest/gtest.h"


TEST(TestScalarFunc, ParseDerivative)
{
  const char* func1 = "sin(1.5*t)*t";
  const char* func2 = "sin(1.5*t)*t:1.5*cos(1.5*t)*t+sin(1.5*t)";

  ScalarFunc* f1 = utl::parseTimeFunc(func1,"expression");
  ScalarFunc* f2 = utl::parseTimeFunc(func2,"expression");

  ASSERT_TRUE(f1 != nullptr);
  ASSERT_TRUE(f2 != nullptr);

  EXPECT_FALSE(f1->isConstant());
  EXPECT_FALSE(f2->isConstant());

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    std::cout <<"f("<< t <<") = "<< (*f1)(t)
	      <<"  f'("<< t <<") = "<< f1->deriv(t) << std::endl;
    EXPECT_FLOAT_EQ((*f1)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ((*f2)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ(f1->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
    EXPECT_FLOAT_EQ(f2->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
  }
}


TEST(TestEvalFunction, ExtraParam)
{
  EvalFunction f("x*foo");
  f.setParam("foo", 2.0);
  Vec3 X(1.0,0.0,0.0);
  EXPECT_FLOAT_EQ(f(X), 2.0);
  X.x = 0.5;
  f.setParam("foo", 4.0);
  EXPECT_FLOAT_EQ(f(X), 2.0);
}


TEST(TestEvalFunction, isConstant)
{
  EXPECT_TRUE (EvalFunction("2.0*x*y").isConstant());
  EXPECT_FALSE(EvalFunction("2.0*x*t").isConstant());
  EXPECT_TRUE (EvalFunction("1.8*tan(x)*x").isConstant());
  EXPECT_FALSE(EvalFunction("2.0*x*tan(t*3)+y").isConstant());
}


TEST(TestEvalFunction, Derivatives)
{
  const char* g    = "sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_x  = "cos(x)*sin(y)*sin(z)*sin(t)";
  const char* g_y  = "sin(x)*cos(y)*sin(z)*sin(t)";
  const char* g_z  = "sin(x)*sin(y)*cos(z)*sin(t)";
  const char* g_t  = "sin(x)*sin(y)*sin(z)*cos(t)";
  const char* g_xx = "-sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_xy = "cos(x)*cos(y)*sin(z)*sin(t)";
  const char* g_xz = "cos(x)*sin(y)*cos(z)*sin(t)";
  const char* g_yy = "-sin(x)*sin(y)*sin(z)*sin(t)";
  const char* g_yz = "sin(x)*cos(y)*cos(z)*sin(t)";
  const char* g_zz = "-sin(x)*sin(y)*sin(z)*sin(t)";

  EvalFunction f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);
  f.addDerivative(g_z,"",3);
  f.addDerivative(g_t,"",4);
  f.addDerivative(g_xx,"",1,1);
  f.addDerivative(g_xy,"",1,2);
  f.addDerivative(g_xz,"",1,3);
  f.addDerivative(g_yy,"",2,2);
  f.addDerivative(g_yz,"",2,3);
  f.addDerivative(g_zz,"",3,3);

  const Vec4 X(1.0,2.0,3.0,4.0);
  for (const double t : {0.1, 0.2, 0.3})
    for (const double x : {0.4, 0.5, 0.6})
        for (const double y : {0.7, 0.8, 0.9})
          for (const double z : {1.0, 1.1, 1.2}) {
            const Vec4 X(x,y,z,t);
            EXPECT_DOUBLE_EQ(f(X), sin(x)*sin(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.deriv(X,1), cos(x)*sin(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.deriv(X,2), sin(x)*cos(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.deriv(X,3), sin(x)*sin(y)*cos(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.deriv(X,4), sin(x)*sin(y)*sin(z)*cos(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,1,1), -sin(x)*sin(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,1,2), cos(x)*cos(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,1,3), cos(x)*sin(y)*cos(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,2,1), cos(x)*cos(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,2,2), -sin(x)*sin(y)*sin(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,2,3), sin(x)*cos(y)*cos(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,3,1), cos(x)*sin(y)*cos(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,3,2), sin(x)*cos(y)*cos(z)*sin(t));
            EXPECT_DOUBLE_EQ(f.dderiv(X,3,3), -sin(x)*sin(y)*sin(z)*sin(t));
        }
}
