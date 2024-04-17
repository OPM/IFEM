// $Id$
//==============================================================================
//!
//! \file MPC.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of multi-point constraint (%MPC) equations.
//!
//==============================================================================

#include "MPC.h"
#include "MPCLess.h"
#include "Utilities.h"
#include <cmath>


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


size_t MPC::getNoMaster (bool recursive) const
{
  size_t nMaster = master.size();

  if (recursive)
    for (const DOF& mdof : master)
      if (mdof.nextc)
        nMaster += mdof.nextc->getNoMaster(true) - 1;

  return nMaster;
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

  // Do not merge if all masters are identical, which may happen if the same
  // MPC has been added for an interface slave node from multiple patches
  bool isIdentic = master.size() == mpc->master.size();
  for (size_t i = 0; i < master.size() && isIdentic; i++)
    isIdentic = (master[i] == mpc->master[i] &&
                 fabs(master[i].coeff - mpc->master[i].coeff) < 1.0e-8);
  if (isIdentic) return false;

  for (const DOF& mdof : mpc->master)
    this->addMaster(mdof);

  return true;
}


int MPC::renumberNodes (const std::map<int,int>& old2new, bool msg)
{
  int invalid = utl::renumber(slave.node,old2new,msg) ? 0 : 1;
  for (DOF& mdof : master)
    if (!utl::renumber(mdof.node,old2new,msg))
      invalid++;

  return invalid;
}


void MPC::shiftNodes (int nshift)
{
  slave.node += nshift;
  for (DOF& mdof : master)
    mdof.node += nshift;
}


void MPC::printMaster (std::ostream& os) const
{
  if (slave.coeff != Real(0))
    os << slave.coeff;
  else if (master.empty())
    os << Real(0);

  size_t count = 0;
  for (const MPC::DOF& mdof : master)
  {
    if (++count == 1 && slave.coeff == Real(0))
      os << mdof.coeff;
    else if (mdof.coeff >= Real(0))
      os <<" + "<< mdof.coeff;
    else
      os <<" - "<< -mdof.coeff;
    os <<"*("<< mdof.node <<","<< mdof.dof;
#if SP_DEBUG > 2
    if (mdof.nextc)
    {
      // Recursively print out the chained constraint
      os <<" (";
      mdof.nextc->printMaster(os);
      os <<")";
    }
#endif
    os <<")";
  }
}


//! \brief Global stream operator printing a constraint equation.

std::ostream& operator<< (std::ostream& s, const MPC& mpc)
{
  s <<"Slave "<< mpc.slave.node <<","<< mpc.slave.dof <<" = ";
  mpc.printMaster(s);
  return s << std::endl;
}
