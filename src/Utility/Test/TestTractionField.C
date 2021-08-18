//==============================================================================
//!
//! \file TestTractionField.C
//!
//! \date Aug 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for evaluation of traction field functions.
//!
//==============================================================================

#include "TractionField.h"
#include "Functions.h"
#include "Vec3Oper.h"

#include "gtest/gtest.h"


/*!
  \brief A simple linearly varying rotation state.
*/

class TestAngleFunc : public VecTimeFunc
{
  Real Roll;  //!< X-rotation gradient
  Real Pitch; //!< Y-rotation gradient
  Real Yaw;   //!< Z-rotation gradient

public:
  //! \brief The constructor initializes the rotation gradients.
  TestAngleFunc(Real rx, Real ry, Real rz) : Roll(rx), Pitch(ry), Yaw(rz) {}
  //! \brief Evaluates the function at time \a t.
  virtual Vec3 evaluate(const Real& t) const
  {
    return Vec3(Roll*t, Pitch*t, Yaw*t);
  }
};


TEST(TestForceDirField, Evaluate)
{
  const Real A = 10, freq = 0.5;        // load magnitude and frequency
  const Real a = 0.2, b = 0.4, c = 0.1; // load orientation angles
  const Real L = 30;                    // half of load domain length

  SineFunc       force(A,freq);
  QuadraticXFunc shape(1,-L,L);
  VecTimeFunc*   norot(new TestAngleFunc(0,0,0));
  VecTimeFunc*   rot_X(new TestAngleFunc(a,0,0));
  VecTimeFunc*   rot_Y(new TestAngleFunc(0,b,0));
  VecTimeFunc*   rot_Z(new TestAngleFunc(0,0,c));
  VecTimeFunc*   rotAl(new TestAngleFunc(a,b,c));
  ForceDirField  f0(new SineFunc(force), norot, new QuadraticXFunc(shape));
  ForceDirField  fx(new SineFunc(force), rot_X, new QuadraticXFunc(shape));
  ForceDirField  fy(new SineFunc(force), rot_Y, new QuadraticXFunc(shape));
  ForceDirField  fz(new SineFunc(force), rot_Z, new QuadraticXFunc(shape));
  ForceDirField  fa(new SineFunc(force), rotAl, new QuadraticXFunc(shape));
  auto&& Shape = [&shape](Real x) { return shape(Vec3(x,0,0)); };

  Vec3 Zero, nVec(1,0,0);

  // Verify that all functions evaluate to zero at (X,t) = ({0,0,0},0)
  EXPECT_EQ(f0(Zero,nVec), Zero);
  EXPECT_EQ(fx(Zero,nVec), Zero);
  EXPECT_EQ(fy(Zero,nVec), Zero);
  EXPECT_EQ(fz(Zero,nVec), Zero);
  EXPECT_EQ(fa(Zero,nVec), Zero);

  // Check function evaluation at origin but with nonzero time (t=0.8)
  Vec4 X(nullptr,0.8);
  Real Fmax  = force(X.t);
  Vec3 Fvec(0,Fmax,0);
  Real alpha = a*X.t;
  Real beta  = b*X.t;
  Real gamma = c*X.t;
  Vec3 Frot  = Tensor(alpha,beta,gamma)*Fvec;
  EXPECT_EQ(f0(X,nVec), Fvec);
  EXPECT_EQ(fx(X,nVec), Vec3( 0, Fmax*cos(alpha), Fmax*sin(alpha) ));
  EXPECT_EQ(fy(X,nVec), Fvec);
  EXPECT_EQ(fz(X,nVec), Vec3(-Fmax*sin(gamma), Fmax*cos(gamma), 0 ));
  EXPECT_EQ(fa(X,nVec), Frot);

  // Check function evaluation at x=0.5*L ==> a scaling by 0.75
  Real x = 0.5*L;
  X.assign(Vec3(x,0,0));
  Fmax *= 0.75;
  EXPECT_EQ(shape(X),0.75);
  EXPECT_EQ(f0(X,nVec), Fvec * 0.75 );
  EXPECT_EQ(fx(X,nVec), Vec3( 0, Fmax*cos(alpha), Fmax*sin(alpha) ));
  EXPECT_EQ(fy(X,nVec), Fvec * Shape(x*cos(beta)) );
  Fmax *= Shape(x*cos(gamma))/0.75;
  EXPECT_EQ(fz(X,nVec), Vec3(-Fmax*sin(gamma), Fmax*cos(gamma), 0 ));
  Frot *= Shape(x*cos(beta)*cos(gamma));
  EXPECT_EQ(fa(X,nVec), Frot);
}
