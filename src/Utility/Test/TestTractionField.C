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

#include "Catch2Support.h"


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


TEST_CASE("TestForceDirField.Evaluate")
{
  const Real A = 10, freq = 0.5;        // load magnitude and frequency
  const Real a = 3.2, b = 7.4, c = 2.1; // load orientation angles (in degrees)
  const Real L = 30;                    // half of load domain length
  const Real Rad = M_PI/Real(180);      // Deg-to-Rad convertion

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
  ForceDirField  fY(new SineFunc(force), new LinearFunc(b), 'Y',
                    new QuadraticXFunc(shape));
  auto&& Shape = [&shape](Real x) { return shape(Vec3(x,0,0)); };

  Vec3 Zero, nVec(1,0,0);

  // Verify that all functions evaluate to zero at (X,t) = ({0,0,0},0)
  REQUIRE(f0(Zero,nVec) == Zero);
  REQUIRE(fx(Zero,nVec) == Zero);
  REQUIRE(fy(Zero,nVec) == Zero);
  REQUIRE(fz(Zero,nVec) == Zero);
  REQUIRE(fa(Zero,nVec) == Zero);
  REQUIRE(fY(Zero,nVec) == Zero);

  // Check function evaluation at origin but with nonzero time (t=0.8)
  Vec4 X(nullptr,0.8);
  Real Fmax  = force(X.t);
  Vec3 Fvec(0,Fmax,0);
  Real alpha = a*X.t*Rad;
  Real beta  = b*X.t*Rad;
  Real gamma = c*X.t*Rad;
  Vec3 Frot  = Tensor(alpha,beta,gamma)*Fvec;
  REQUIRE(A*sin(freq*X.t) ==  Fmax);
  REQUIRE(f0(X,nVec) ==  Fvec);
  REQUIRE(fx(X,nVec) == Vec3( 0, Fmax*cos(alpha), Fmax*sin(alpha) ));
  REQUIRE(fy(X,nVec) == Fvec);
  REQUIRE(fz(X,nVec) == Vec3(-Fmax*sin(gamma), Fmax*cos(gamma), 0 ));
  REQUIRE(fa(X,nVec) == Frot);
  REQUIRE(fY(X,nVec) == Fvec);

  // Check function evaluation at x=0.5*L ==> a scaling by 0.75
  Real x = 0.5*L;
  X.assign(Vec3(x,0,0));
  Fmax *= 0.75;
  REQUIRE(shape(X) == 0.75);
  REQUIRE(f0(X,nVec) == Fvec * 0.75 );
  REQUIRE(fx(X,nVec) == Vec3( 0, Fmax*cos(alpha), Fmax*sin(alpha) ));
  REQUIRE(fy(X,nVec) == Fvec * Shape(x*cos(beta)) );
  Fmax *= Shape(x*cos(gamma))/0.75;
  REQUIRE(fz(X,nVec) == Vec3(-Fmax*sin(gamma), Fmax*cos(gamma), 0 ));
  Frot *= Shape(x*cos(beta)*cos(gamma));
  REQUIRE(fa(X,nVec) == Frot);
  REQUIRE(fY(X,nVec) == Fvec * Shape(x*cos(beta)) );
}
