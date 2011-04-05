// $Id$
//==============================================================================
//!
//! \file Utilities.C
//!
//! \date May 29 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility methods.
//!
//==============================================================================

#include "Utilities.h"
#include <ctype.h>


char* utl::readLine (std::istream& is)
{
  static char buf[1024];
  for (size_t n = 0; is.good(); n = strlen(buf))
  {
    for (size_t i = 0; i < n; i++)
      if (!isspace(buf[i]))
	if (buf[i] == '#')
	  break; // comment line - skip
	else
	  return buf+i; // first non-blank character

    is.getline(buf,sizeof(buf));
  }

  return 0;
}


bool utl::ignoreComments (std::istream& is)
{
  int c = ' ';
  while (isspace(c))
    c = is.get();

  while (c == '#')
  {
    is.ignore(256,'\n');
    c = is.get();
  }
  is.putback(c);

  return is.good();
}


bool utl::renumber (int& num, int& runner, std::map<int,int>& old2new)
{
  std::map<int,int>::iterator it = old2new.find(num);
  if (it == old2new.end())
    it = old2new.insert(std::make_pair(num,++runner)).first;

  if (num == it->second) return false;

  num = it->second;
  return true;
}


bool utl::renumber (int& num, const std::map<int,int>& old2new, bool msg)
{
  std::map<int,int>::const_iterator it = old2new.find(num);
  if (it == old2new.end())
  {
    if (msg)
      std::cerr <<" *** utl::renumber: Old value "<< num
		<<" does not exist in old2new mapping"<< std::endl;
    return false;
  }

  num = it->second;
  return true;
}


int utl::gather (const std::vector<int>& index, size_t nr,
		 const std::vector<real>& in, std::vector<real>& out,
		 size_t offset_in)
{
  int outside = 0;
  out.resize(nr*index.size());
  const real* data = &in.front() + offset_in;
  real* outVec = &out.front();
  for (size_t i = 0; i < index.size(); i++, outVec += nr)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      memcpy(outVec,data+nr*index[i],nr*sizeof(real));
    else
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t nr,
		 const utl::vector<real>& in, utl::matrix<real>& out,
		 size_t offset_in)
{
  int outside = 0;
  out.resize(nr,index.size());
  const real* data = &in.front() + offset_in;
  for (size_t i = 0; i < index.size(); i++)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      out.fillColumn(1+i,data+nr*index[i]);
    else
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t ir, size_t nr,
                 const std::vector<real>& in, std::vector<real>& out,
                 size_t offset_in, int shift_idx)
{
  if (ir >= nr) return index.size();

  int outside = 0;
  out.resize(index.size());
  offset_in += ir;
  for (size_t i = 0; i < index.size(); i++)
    if (index[i] >= shift_idx)
    {
      size_t ip = offset_in + nr*(index[i]-shift_idx);
      if (ip < in.size())
	out[i] = in[ip];
      else
	outside++;
    }
    else
      outside++;

  return outside;
}


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


void utl::Pascal (int p, real x, real y, std::vector<real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,2));
  phi.push_back(real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
    {
      int k, j = q-i;
      real a = real(1);
      for (k = 0; k < i; k++) a *= x;
      for (k = 0; k < j; k++) a *= y;
      phi.push_back(a);
    }
}


void utl::Pascal (int p, real x, real y, real z, std::vector<real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,3));
  phi.push_back(real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
      for (int j = q; j >= 0; j--)
	if (i+j <= q)
	{
	  int l, k = q-i-j;
	  real a = real(1);
	  for (l = 0; l < i; l++) a *= x;
	  for (l = 0; l < j; l++) a *= y;
	  for (l = 0; l < k; l++) a *= z;
	  phi.push_back(a);
	}
}


size_t utl::find_closest (const std::vector<real>& a, real v)
{
  // The lower_bound function uses binary search to find the index of the value.
  // Thus, this works only when the vector a is sorted in increasing order.
  size_t i = std::lower_bound(a.begin(),a.end(),v) - a.begin();
  if (i > 0 && i < a.size())
    return a[i]-v < v-a[i-1] ? i : i-1;
  else if (i == 0)
    return 0;
  else
    return i-1;
}
