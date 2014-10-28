// $Id$
//==============================================================================
//!
//! \file NewmarkNLSIM.h
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#ifndef _NEWMARK_NL_SIM_H
#define _NEWMARK_NL_SIM_H

#include "NewmarkSIM.h"

class SystemVector;


/*!
  \brief Newmark-based solution driver for dynamic isogeometric FEM simulators.

  \details This class implements the Hilber-Hughes-Taylor (HHT) time integration
  algorithm, with a constant displacement predictor step. The algorithm relies 
  on the following three algorithm parameters:

  \f$\alpha_H\f$ = parameter controlling the amount of numerical damping<BR>
  \f$\alpha_1\f$ = mass-proportional damping factor<BR>
  \f$\alpha_2\f$ = stiffness-proportional damping factor<BR>

  For a given time discretization,
  <I>t<SUB>0</SUB></I> = 0 and <I>t<SUB>n</SUB></I> =
  <I>t<SUB>n-1</SUB></I> + \f$\Delta\f$<I>t<SUB>n</SUB></I>
  for <I>n</I>=1...<I>n<SUB>step</SUB></I>,
  the HHT time integration algorithm goes like this:

  <H3>Initialisation of time step loop:</H3>
  \f[ \beta=\frac{1}{4}(1.0-\alpha_H)^2,\; \gamma=\frac{1}{2}-\alpha_H \f]
  \f[ \begin{array}{lcll}
    t &=& 0 & \mbox{(initial time)} \\
    {\bf u}_0 &=& {\bf 0} & \mbox{(initial displacements)} \\
    \dot{\bf u}_0 &=& {\bf 0} & \mbox{(initial velocity)} \\
    \ddot{\bf u}_0 &=& {\bf 0} & \mbox{(initial acceleration)} \\
    \tilde{\bf R}_0 &=& {\bf 0} & \mbox{(initial actual inertia force)} \\
  \end{array} \f]

  <H3>FOR <I>n</I> = 1 TO <I>n<SUB>step</SUB></I> DO</H3>
  <OL>
    <LI>Initialisation of iteration loop:
      \f{eqnarray*}{
        i &=& 0 \quad\mbox{(iteration counter)}\\
        t &=& t + \Delta t_n \\
        {\bf u}_n^0 &=& {\bf u}_{n-1} \\
        \dot{\bf u}_n^0 &=& \dot{\bf u}_{n-1} \\
        \ddot{\bf u}_n^0 &=& \ddot{\bf u}_{n-1}
      \f}
    </LI>

    <LI>Predict new velocity and acceleration:
      \f{eqnarray*}{
        {\bf v}_n &=& (\frac{\gamma}{\beta}-1)\dot{\bf u}_{n-1} +
                      \Delta t_n(\frac{\gamma}{2\beta}-1)\ddot{\bf u}_{n-1} \\
        {\bf a}_n &=& (\frac{1}{2\beta}-1)\ddot{\bf u}_{n-1} +
                      \frac{1}{\Delta t_n\beta}\dot{\bf u}_{n-1}
      \f}\f[ \begin{array}{lcll}
        \dot{\bf u}_n^0 &=& {\bf v}_n & \mbox{(predicted velocity)} \\
        \ddot{\bf u}_n^0 &=& {\bf a}_n & \mbox{(predicted acceleration)} \\
        \Delta{\bf u}_n &=& {\bf 0} & \mbox{(displacement increment)}
        \end{array}
      \f]
    </LI>

    <LI>Assemble FE matrices and right-hand-side force vectors:
      \f[ \begin{array}{lcll}
        {\bf K}_n^0 &=& {\bf K}({\bf u}_n^0,t) &
        \mbox{(tangential stiffness matrix)} \\
        {\bf C}_n^0 &=& {\bf C}({\bf u}_n^0,t) &
        \mbox{(damping matrix)} \\
        {\bf M}_n^0 &=& {\bf M}({\bf u}_n^0,t) &
        \mbox{(mass matrix)} \\
        {\bf F}_n^{S,0} &=& {\bf F}^S({\bf u}_n^0,t) &
        \mbox{(internal stiffness forces)} \\
        {\bf F}_n^{I,0} &=& {\bf F}^I({\bf u}_n^0,\ddot{\bf u}_n^0,t) &
        \mbox{(internal inertia forces)} \\
        {\bf F}_n^{E,0} &=& {\bf F}^E({\bf u}_n^0,t) &
        \mbox{(external forces)}
        \end{array}
      \f]
    </LI>

    <LI>Compute Newton matrix and associated incremental load vector:
      \f{eqnarray*}{
        {\bf N}_n^0 &=& a_n{\bf M}_n^0 + b_n{\bf C}_n^0 + c_n{\bf K}_n^0 \\
        {\bf R}_n^0 &=& (1+\alpha_H)\left[{\bf F}_n^{E,0} - {\bf F}_n^{S,0} +
                        (\alpha_1{\bf M}_n^0+\alpha_2{\bf K}_n^0){\bf v}_n
                        \right] + {\bf F}_n^{I,0} - \alpha_H\tilde{\bf R}_{n-1}
      \f} where \f{eqnarray*}{
      a_n &=& \frac{1}{\Delta t_n^2\beta} +
              (1+\alpha_H)\frac{\alpha_1\gamma}{\Delta t_n\beta} \\
      b_n &=& (1+\alpha_H)\frac{\gamma}{\Delta t_n\beta} \\
      c_n &=& (1+\alpha_H)\left(1+\frac{\alpha_2\gamma}{\Delta t_n\beta}\right)
      \f} and \f[
        \tilde{\bf R}_n \;=\; {\bf F}_n^E - {\bf F}_n^S -
        (\alpha_1{\bf M}_n+\alpha_2{\bf K}_n)\dot{\bf u}_n
        \quad\mbox{(actual inertia force in time step $n > 0$)}
      \f]
    </LI>

    <LI>Solve for the incremental displacement:
      \f[ {\bf N}_n^0 \Delta{\bf u}_n^0 \;=\; {\bf R}_n^0
        \quad\Longrightarrow\; \Delta{\bf u}_n^0
      \f]
    </LI>

    <LI><H4>WHILE <I>i</I> < <I>MaxIt</I> AND
            \f$\|\Delta{\bf u}_n^i\|>\epsilon_{\rm tol}\f$ DO</H4>
    <OL>
      <LI>Update configuration:
        \f{eqnarray*}{
          i &=& i + 1 \\
          \Delta{\bf u}_n &=& \Delta{\bf u}_n + \Delta{\bf u}_n^{i-1} \\
          {\bf u}_n^i &=& {\bf u}_{n-1} \oplus \Delta{\bf u}_n \\
          \dot{\bf u}_n^i &=& \frac{\gamma}{\Delta t_n\beta}\Delta{\bf u}_n
                              - {\bf v}_n \\
          \ddot{\bf u}_n^i &=& \frac{1}{\Delta t_n^2\beta}\Delta{\bf u}_n
                              - {\bf a}_n
        \f}
      </LI>

      <LI>Assemble FE matrices and right-hand-side force vectors:
        \f[ \begin{array}{lll}
          {\bf K}_n^i \;=\; {\bf K}({\bf u}_n^i,t)\;, &
          {\bf C}_n^i \;=\; {\bf C}({\bf u}_n^i,t)\;, &
          {\bf M}_n^i \;=\; {\bf M}({\bf u}_n^i,t) \\[1mm]
          {\bf F}_n^{S,i} =\; {\bf F}^S({\bf u}_n^i,t)\;, &
          {\bf F}_n^{I,i} \;=\; {\bf F}^I({\bf u}_n^i,\ddot{\bf u}_n^i,t)\;, &
          {\bf F}_n^{E,i} \;=\; {\bf F}^E({\bf u}_n^i,t)
        \end{array} \f]
      </LI>

      <LI>Compute Newton matrix and associated incremental load vector:
        \f{eqnarray*}{
          {\bf N}_n^i &=& a_n{\bf M}_n^i + b_n{\bf C}_n^i + c_n{\bf K}_n^i \\
          {\bf R}_n^i &=& (1+\alpha_H)\left[{\bf F}_n^{E,i} - {\bf F}_n^{S,i} -
                          (\alpha_1{\bf M}_n^i+\alpha_2{\bf K}_n^i){\bf v}_n
                          \right] + {\bf F}_n^{I,i} -\alpha_H\tilde{\bf R}_{n-1}
        \f}
      </LI>

      <LI>Solve for the iterative displacement:
        \f[ {\bf N}_n^i \Delta{\bf u}_n^i \;=\; {\bf R}_n^i
          \quad\Longrightarrow\; \Delta{\bf u}_n^i
        \f]
      </LI>
    </OL>
    <H4>END DO</H4></LI>

    <LI>Update converged configuration:
      \f{eqnarray*}{
      \Delta{\bf u}_n &=& \Delta{\bf u}_n + \Delta{\bf u}_n^i \\
      {\bf u}_n &=& {\bf u}_{n-1} \oplus \Delta{\bf u}_n \\
      \dot{\bf u}_n&=&\frac{\gamma}{\Delta t_n\beta}\Delta{\bf u}_n-{\bf v}_n \\
      \ddot{\bf u}_n&=&\frac{1}{\Delta t_n^2\beta}\Delta{\bf u}_n - {\bf a}_n
      \f}
    </LI>
  </OL>
  <H3>END DO</H3>
*/

class NewmarkNLSIM : public NewmarkSIM
{
public:
  //! \brief The constructor initializes default solution parameters.
  NewmarkNLSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~NewmarkNLSIM() {}

  //! \brief Parses a data section from an XML document.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out problem-specific data to the given stream.
  virtual void printProblem(std::ostream& os) const;

  //! \brief Initializes primary solution vectors and integration parameters.
  virtual void init(size_t nSol = 3);

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Modifies the current solution vector (used by sub-iterations only).
  virtual void setSolution(const Vector& newSol, int idx);

protected:
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param, bool converged);
  //! \brief Finalizes the right-hand-side vector on the system level.
  virtual void finalizeRHSvector();

private:
  Vector incDis;  //!< Incremental displacements
  Vector predVel; //!< Predicted velocity vector
  Vector predAcc; //!< Predicted acceleration vector

  SystemVector* Finert; //!< Actual inertia forces in last converged time step
};

#endif
