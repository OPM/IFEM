// $Id$
//==============================================================================
//!
//! \file MPC.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of multi-point constraint (MPC) equations.
//!
//==============================================================================

#include "MPC.h"
#include "MPCLess.h"
#include "Utilities.h"


bool MPCLess::compareSlaveDofOnly = false;


//! \brief Global operator for comparing two MPC::DOF objects.

bool operator< (const MPC::DOF& lhs, const MPC::DOF& rhs)
{
  if (lhs.node < rhs.node) return true;
  if (lhs.node > rhs.node) return false;
  if (lhs.dof  < rhs.dof)  return true;
  if (lhs.dof  > rhs.dof)  return false;

  if (MPCLess::compareSlaveDofOnly)
    return false; // ignore coefficient differences, if any
  else
    return lhs.coeff < rhs.coeff ? true : false;
}


bool MPCLess::operator() (const MPC* lhs, const MPC* rhs) const
{
  if (!rhs) return false;
  if (!lhs) return true;

  if (lhs->getSlave() < rhs->getSlave()) return true;
  if (rhs->getSlave() < lhs->getSlave()) return false;

  if (compareSlaveDofOnly) return false; // ignore master differences, if any

  size_t lMaster = lhs->getNoMaster();
  size_t rMaster = rhs->getNoMaster();
  for (size_t i = 0; i < lMaster && i < rMaster; i++)
  {
    if (lhs->getMaster(i) < rhs->getMaster(i)) return true;
    if (rhs->getMaster(i) < lhs->getMaster(i)) return false;
  }

  return lMaster < rMaster ? true : false;
}


void MPC::addMaster (const DOF& dof, Real tol)
{
  if (dof.coeff >= -tol && dof.coeff <= tol) return; // ignore if zero coeff.

  std::vector<DOF>::iterator it = std::find(master.begin(),master.end(),dof);
  if (it == master.end())
    master.push_back(dof);
  else
    it->coeff += dof.coeff;
}


bool MPC::merge (const MPC* mpc)
{
  if (!(this->slave == mpc->slave && this->slave.coeff == mpc->slave.coeff))
    return false; // the slave definitions did not match

  for (size_t i = 0; i < mpc->master.size(); i++)
    this->addMaster(mpc->master[i]);

  return true;
}


int MPC::renumberNodes (const std::map<int,int>& old2new, bool msg)
{
  int invalid = utl::renumber(slave.node,old2new,msg) ? 0 : 1;
  for (size_t i = 0; i < master.size(); i++)
    if (!utl::renumber(master[i].node,old2new,msg))
      invalid++;

  return invalid;
}


//! \brief Global stream operator printing a constraint equation.

std::ostream& operator<< (std::ostream& s, const MPC& mpc)
{
  s <<"Slave "<< mpc.slave.node <<","<< mpc.slave.dof;
  if (mpc.slave.coeff != Real(0))
    s <<" = "<< mpc.slave.coeff;
  else if (mpc.master.empty())
    return s <<" = 0"<< std::endl;
  else
    s <<" = ";

  for (size_t i = 0; i < mpc.master.size(); i++)
  {
    if (i == 0 && mpc.slave.coeff == Real(0))
      s << mpc.master[i].coeff;
    else if (mpc.master[i].coeff >= Real(0))
      s <<" + "<< mpc.master[i].coeff;
    else
      s <<" - "<< -mpc.master[i].coeff;
    s <<"*("<< mpc.master[i].node <<","<< mpc.master[i].dof <<")";
  }
  return s << std::endl;
}
