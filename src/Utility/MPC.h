// $Id: MPC.h,v 1.4 2010-10-14 19:10:56 kmo Exp $
//==============================================================================
//!
//! \file MPC.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of multi-point constraint (MPC) equations.
//!
//==============================================================================

#ifndef _MPC_H
#define _MPC_H

#include <map>
#include <vector>
#include <iostream>


/*!
  \brief A class for representing a general multi-point constraint equation.

  \details A multi-point constraint (MPC) equation is used to introduce
  additional coupling between the degrees of freedom (DOF) in a FE grid.

  An MPC equation is in general a linear coupling of one (slave) dof to a set
  of master dofs, and may be expressed as

  \f[ s = c_0 + c_1 m_1 + c_2 m_2 + \ldots + c_n m_n \f]

  where
  <ul>
  <li> \f$s\f$ denotes the slave dof
  <li> \f$c_0\f$ is the value of the slave dof when all master dofs are zero,
       or when there are no masters in the equation (prescribed dof).
  <li> \f$m_i\f$ is a master dof
  <li> \f$c_i, i > 0\f$ is a constant scaling coefficient
       associated with the \a i'th master
  <li> \a n denotes the total number of master dofs in the equation
  </ul>
  When \a n = 0, the above equation represents a fixed (\f$c_0=0\f$) or
  prescribed (\f$c_0\ne0\f$) dof.

  One or more of the master dofs may also be a slave in another constraint
  equation (chained constrains). In this case the two equations are combined
  to eliminate the master dof that is constrained. This is done while
  preprocessing the model by the Grid::resolveMPCchains function.
*/

class MPC
{
public:

  /*!
    \brief A struct for representing one term (DOF number and associated
    coefficient) in a MPC equation.
  */
  struct DOF
  {
    //! \brief Default constructor with no arguments.
    DOF()
    {
      node = dof = 0;
      coeff = real(0);
    };

    //! \brief Convenience constructor creating a valid DOF object.
    //! \param[in] n Node number (1...NNOD)
    //! \param[in] d The local DOF number (1...3)
    //! \param[in] c Associated coefficient or constrained value
    DOF(int n, int d, real c = real(0))
    {
      node  = n;
      dof   = d;
      coeff = c;
    };

    //! \brief Global stream operator printing a DOF instance.
    friend std::ostream& operator<<(std::ostream& s, const DOF& dof)
    {
      return s <<"u_"<< char('w'+dof.dof) <<" in node "<< dof.node;
    }

    int  node;  //!< Node number identifying this DOF
    int  dof;   //!< Local DOF number within \a node
    real coeff; //!< The constrained value, or master DOF scaling coefficient
  };

  //! \brief Constructor creating a constraint for a specified slave DOF
  //! with no master DOFs.
  //! \param[in] n The node number of the slave DOF (1...NNOD)
  //! \param[in] d The local DOF number of the slave DOF (1...3)
  //! \param[in] c The actual value that this slave DOF is constrained to
  //! when there are no master DOFs, or all master DOFs are zero
  MPC(int n, int d, real c = real(0)) { iceq = -1; slave = DOF(n,d,c); }

  //! \brief Adds a master DOF to the constraint equation.
  //! \param[in] n The node number of the master DOF (1...NNOD)
  //! \param[in] d The local DOF number of the master DOF (1...3)
  //! \param[in] c The coefficient that this master should be scaled with
  //! \param[in] tol Tolerance for comparison with zero,
  //! if the coefficient \a c is zero, the master DOF is not added
  void addMaster(int n, int d, real c = real(1), real tol = real(1.0e-8))
  {
    if (c < -tol || c > tol) master.push_back(DOF(n,d,c));
  }

  //! \brief Updates the coefficient of the \a pos'th master DOF.
  void updateMaster(size_t pos, real c)
  {
    if (pos < master.size())
      master[pos].coeff = c;
  }

  //! \brief Removes the \a pos'th master DOF from the constraint equation.
  void removeMaster(size_t pos)
  {
    if (pos < master.size())
      master.erase(master.begin()+pos);
  }

  //! \brief Renumbers all node numbers in the constraint equation.
  //! \param[in] old2new Old-to-new node number mapping
  //! \param[in] msg If \e true, give error message for each invalid node number
  //! \return Number of invalid node numbers detected
  int renumberNodes(const std::map<int,int>& old2new, bool msg = false);

  //! \brief Increments the \a c0 coefficient by a given \a offset.
  void addOffset(real offset) { slave.coeff += offset; }

  //! \brief Assigns a new \a c0 coefficient to the constraint equation.
  void setSlaveCoeff(real c0) { slave.coeff = c0; }

  //! \brief Returns a reference to the slave DOF.
  const DOF& getSlave() const          { return slave; }
  //! \brief Returns a reference to the \a i'th master DOF.
  const DOF& getMaster(size_t i) const { return master[i]; }

  //! \brief Returns the number of master DOFs.
  size_t getNoMaster() const { return master.size(); }

  //! \brief Global stream operator printing a constraint equation.
  friend std::ostream& operator<<(std::ostream& s, const MPC& mpc);

  int              iceq;   //!< Global constraint equation identifier
private:
  DOF              slave;  //!< The slave DOF of this constraint equation
  std::vector<DOF> master; //!< The master DOFs of this constraint equation
};

#endif
