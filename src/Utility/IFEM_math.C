// $Id$
//==============================================================================
//!
//! \file IFEM_math.C
//!
//! \date Oct 20 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various math utility methods.
//!
//==============================================================================

#include "IFEM_math.h"


size_t utl::Pascal (int p, unsigned short int nsd)
{
  size_t nM = 1;
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
      if (nsd == 2)
        nM++;
      else for (int j = q; j >= 0; j--)
        if (i+j <= q) nM++;

  return nM;
}


void utl::Pascal (int p, Real x, Real y, std::vector<Real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,2));
  phi.push_back(Real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
    {
      int k, j = q-i;
      Real a = Real(1);
      for (k = 0; k < i; k++) a *= x;
      for (k = 0; k < j; k++) a *= y;
      phi.push_back(a);
    }
}


void utl::Pascal (int p, Real x, Real y, Real z, std::vector<Real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,3));
  phi.push_back(Real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
      for (int j = q; j >= 0; j--)
        if (i+j <= q)
        {
          int l, k = q-i-j;
          Real a = Real(1);
          for (l = 0; l < i; l++) a *= x;
          for (l = 0; l < j; l++) a *= y;
          for (l = 0; l < k; l++) a *= z;
          phi.push_back(a);
        }
}
