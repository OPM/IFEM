// $Id$
//==============================================================================
//!
//! \file ThreadGroups.C
//!
//! \date May 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Threading group partitioning.
//!
//==============================================================================

#include "ThreadGroups.h"
#include <algorithm>
#include <numeric>
#include <iostream>
#ifdef USE_OPENMP
#include <omp.h>
#endif


void ThreadGroups::oneGroup (size_t nel)
{
  tg[0].resize(1);
  tg[1].resize(0);
  tg[0][0].resize(nel);
  std::iota(tg[0][0].begin(),tg[0][0].end(),0);
}


void ThreadGroups::oneStripe (size_t nel)
{
  tg[0].resize(nel);
  tg[1].resize(0);
  for (size_t iel = 0; iel < nel; iel++)
    tg[0][iel].resize(1,iel);
}


void ThreadGroups::calcGroups (const BoolVec& el1, const BoolVec& el2,
                               int p1, int p2)
{
#ifndef USE_OPENMP
  this->oneGroup(el1.size()*el2.size());
#else
  // Count the non-zero element in each direction, the zero-span elements
  // should not affect the partitioning as they don't involve any work
  size_t nel1 = 0, nel2 = 0;
  for (bool e : el1) if (e) nel1++;
  for (bool e : el2) if (e) nel2++;

  int threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  if (stripDir == ANY)
    stripDir = getStripDirection(nel1,nel2,parts);
  int mul = stripDir == U ? 1 : el1.size();
  int els = stripDir == U ? nel1 : nel2;

  // The minimum stripe size (width) depends on the polynomial degree
  // due to the overlapping support of the splines basis functions
  int stripsize = 0;
  int minsize = stripDir == U ? p1 : p2;
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

  int remainder = 0;
  if (threads > 1)
    remainder = els - stripsize*parts;
  else
    stripsize = els;

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available"
            << "\nnel1 " << nel1 << "\nnel2 " << nel2
            << "\nstripsize " << stripsize
            << "\n# of strips " << els/stripsize
            << "\nremainder " << remainder << std::endl;
#endif

  nel1 = el1.size();
  nel2 = el2.size();
  if (threads == 1)
    this->oneGroup(nel1*nel2);
  else
  {
    int i, j, t, zspan, offs = 0;
    IntVec stripsizes[2], startelms[2];
    const BoolVec& elz = stripDir == U ? el1 : el2;

    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    startelms[0].reserve(threads);
    startelms[1].reserve(threads);
    for (t = j = 0; t < threads; ++t) {
      startelms[0].push_back(offs*mul);
      for (i = zspan = 0; i < stripsizes[0][t]; offs++)
        if (elz[j++])
          i++;
        else
          zspan++;
      stripsizes[0][t] += zspan; // add zero-span elements to this thread

      startelms[1].push_back(offs*mul);
      for (i = zspan = 0; i < stripsizes[1][t]; offs++)
        if (elz[j++])
          i++;
        else
          zspan++;
      stripsizes[1][t] += zspan; // add zero-span elements to this thread
    }

    for (i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = stripDir == U ? stripsizes[i][t] : nel1;
        int maxy = stripDir == V ? stripsizes[i][t] : nel2;
        tg[i][t].reserve(maxx*maxy);
        for (int i2 = 0; i2 < maxy; ++i2)
          for (int i1 = 0; i1 < maxx; ++i1)
            tg[i][t].push_back(startelms[i][t]+i1+i2*nel1);
      }
#if SP_DEBUG > 1
      printGroup(tg[i],i);
#endif
    }
  }
#endif
}


void ThreadGroups::calcGroups (int nel1, int nel2, int minsize)
{
#ifndef USE_OPENMP
  this->oneGroup(nel1*nel2);
#else
  int threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  if (stripDir == ANY)
    stripDir = getStripDirection(nel1,nel2,parts);
  int mul = stripDir == U ? 1 : nel1;
  int els = stripDir == U ? nel1 : nel2;

  int stripsize = 0;
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

  int remainder = 0;
  if (threads > 1)
    remainder = els - stripsize*parts;
  else
    stripsize = els;

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available"
            << "\nnel1 " << nel1 << "\nnel2 " << nel2
            << "\nstripsize " << stripsize
            << "\n# of strips " << els/stripsize
            << "\nremainder " << remainder << std::endl;
#endif

  if (threads == 1)
    this->oneGroup(nel1*nel2);
  else
  {
    int i, offs = 0;
    IntVec stripsizes[2], startelms[2];

    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    for (i = 0; i < threads; ++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }

    for (i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = stripDir == U ? stripsizes[i][t] : nel1;
        int maxy = stripDir == V ? stripsizes[i][t] : nel2;
        for (int i2 = 0; i2 < maxy; ++i2)
          for (int i1 = 0; i1 < maxx; ++i1)
            tg[i][t].push_back(startelms[i][t]+i1+i2*nel1);
      }
#if defined(USE_OPENMP) && SP_DEBUG > 1
      printGroup(tg[i],i);
#endif
    }
  }
#endif
}


void ThreadGroups::calcGroups (const BoolVec& el1, const BoolVec& el2,
                               const BoolVec& el3, int p1, int p2, int p3)
{
#ifndef USE_OPENMP
  this->oneGroup(el1.size()*el2.size()*el3.size());
#else
  // Count the non-zero element in each direction, the zero-span elements
  // should not affect the partitioning as they don't involve any work
  size_t nel1 = 0, nel2 = 0, nel3 = 0;
  for (bool e : el1) if (e) nel1++;
  for (bool e : el2) if (e) nel2++;
  for (bool e : el3) if (e) nel3++;

  int threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  if (stripDir == ANY)
    stripDir = getStripDirection(nel1,nel2,nel3,parts);
  int mul = stripDir == U ? 1 : el1.size()*(stripDir == V ? 1 : el2.size());
  int els = stripDir == U ? nel1 : (stripDir == V ? nel2 : nel3);

  // The minimum stripe size (width) depends on the polynomial degree
  // due to the overlapping support of the splines basis functions
  int stripsize = 0;
  int minsize = stripDir == U ? p1 : (stripDir == V ? p2 : p3);
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

  int remainder = 0;
  if (threads > 1)
    remainder = els - stripsize*parts;
  else
    stripsize = els;

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available"
            << "\nnel1 " << nel1 << "\nnel2 " << nel2 << "\nnel3 " << nel3
            << "\nstripsize " << stripsize
            << "\n# of strips " << els/stripsize
            << "\nremainder " << remainder << std::endl;
#endif

  nel1 = el1.size();
  nel2 = el2.size();
  nel3 = el3.size();
  if (threads == 1)
    this->oneGroup(nel1*nel2*nel3);
  else
  {
    int i, j, t, zspan, offs = 0;
    IntVec stripsizes[2], startelms[2];
    const BoolVec& elz = stripDir == U ? el1 : (stripDir == V ? el2 : el3);

    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    startelms[0].reserve(threads);
    startelms[1].reserve(threads);
    for (t = j = 0; t < threads; ++t) {
      startelms[0].push_back(offs*mul);
      for (i = zspan = 0; i < stripsizes[0][t]; offs++)
        if (elz[j++])
          i++;
        else
          zspan++;
      stripsizes[0][t] += zspan; // add zero-span elements to this thread

      startelms[1].push_back(offs*mul);
      for (i = zspan = 0; i < stripsizes[1][t]; offs++)
        if (elz[j++])
          i++;
        else
          zspan++;
      stripsizes[1][t] += zspan; // add zero-span elements to this thread
    }

    for (i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = stripDir == U ? stripsizes[i][t] : nel1;
        int maxy = stripDir == V ? stripsizes[i][t] : nel2;
        int maxz = stripDir == W ? stripsizes[i][t] : nel3;
        tg[i][t].reserve(maxx*maxy*maxz);
        for (int i3 = 0; i3 < maxz; ++i3)
          for (int i2 = 0; i2 < maxy; ++i2)
            for (int i1 = 0; i1 < maxx; ++i1)
              tg[i][t].push_back(startelms[i][t]+i1+nel1*(i2+nel2*i3));
      }

#if SP_DEBUG > 1
      printGroup(tg[i],i);
#endif
    }
  }
#endif
}


void ThreadGroups::calcGroups (int nel1, int nel2, int nel3, int minsize)
{
#ifndef USE_OPENMP
  this->oneGroup(nel1*nel2*nel3);
#else
  int threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  if (stripDir == ANY)
    stripDir = getStripDirection(nel1,nel2,nel3,parts);
  int mul = stripDir == U ? 1 : nel1*(stripDir == V ? 1 : nel2);
  int els = stripDir == U ? nel1 : (stripDir == V ? nel2 : nel3);

  int stripsize = 0;
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

  int remainder = 0;
  if (threads > 1)
    remainder = els - stripsize*parts;
  else
    stripsize = els;

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available"
            << "\nnel1 " << nel1 << "\nnel2 " << nel2 << "\nnel3 " << nel3
            << "\nstripsize " << stripsize
            << "\n# of strips " << els/stripsize
            << "\nremainder " << remainder << std::endl;
#endif

  if (threads == 1)
    this->oneGroup(nel1*nel2*nel3);
  else
  {
    int offs, i;
    IntVec stripsizes[2];
    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    IntVec startelms[2];
    for (offs = i = 0; i < threads; ++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }

    for (i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = stripDir == U ? stripsizes[i][t] : nel1;
        int maxy = stripDir == V ? stripsizes[i][t] : nel2;
        int maxz = stripDir == W ? stripsizes[i][t] : nel3;
        for (int i3 = 0; i3 < maxz; ++i3)
          for (int i2 = 0; i2 < maxy; ++i2)
            for (int i1 = 0; i1 < maxx; ++i1)
              tg[i][t].push_back(startelms[i][t]+i1+nel1*(i2+nel2*i3));
      }

#if defined(USE_OPENMP) && SP_DEBUG > 1
      printGroup(tg[i],i);
#endif
    }
  }
#endif
}


ThreadGroups::StripDirection ThreadGroups::getStripDirection (int nel1,
                                                              int nel2,
                                                              int parts)
{
  int s1 = nel1 / parts;
  int s2 = nel2 / parts;
  int r1 = nel1 - s1*parts;
  int r2 = nel2 - s2*parts;

  if (r1*nel2 < r2*nel1)
    return U;
  else if (r1*nel2 > r2*nel1)
    return V;

  return nel1 > nel2 ? U : V;
}


ThreadGroups::StripDirection ThreadGroups::getStripDirection (int nel1,
                                                              int nel2,
                                                              int nel3,
                                                              int parts)
{
  int s1 = nel1 / parts;
  int s2 = nel2 / parts;
  int s3 = nel3 / parts;
  int r1 = nel1 - s1*parts;
  int r2 = nel2 - s2*parts;
  int r3 = nel3 - s3*parts;

  if (r1*nel2*nel3 < nel1*r2*nel3 && r1*nel2*nel3 < nel1*nel2*r3)
    return U;
  else if (nel1*r2*nel3 < r1*nel2*nel3 && nel1*r2*nel3 < nel1*nel2*r3)
    return V;
  else if (nel1*nel2*r3 < r1*nel2*nel3 && nel1*nel2*r3 < nel1*r2*nel3)
    return W;

  // The number of left-over elements is not smallest in one direction only
  if (r1*nel2*nel3 > nel1*r2*nel3)
    return nel2 > nel3 ? V : W;
  else if (nel1*r2*nel3 > nel1*nel2*r3)
    return nel1 > nel3 ? U : W;
  else if (nel1*nel2*r3 > r1*nel2*nel3)
    return nel1 > nel2 ? U : V;

  // The number of left-over elements is the same in all three directions
  if (nel1 >= nel2 && nel1 >= nel3)
    return U;
  else if (nel2 >= nel1 && nel2 >= nel3)
    return V;
  else
    return W;
}


void ThreadGroups::applyMap (const IntVec& map)
{
  for (size_t l = 0; l < 2; ++l)
    for (size_t k = 0; k < tg[l].size(); ++k)
      for (size_t j = 0; j < tg[l][k].size(); ++j)
        tg[l][k][j] = map[tg[l][k][j]];
}


void ThreadGroups::printGroup (const IntMat& group, int g)
{
  std::cout <<"group "<< g;
  for (size_t t = 0; t < group.size(); t++)
  {
    std::cout <<"\n\t thread "<< t <<" ("<< group[t].size() <<"):";
    for (int e : group[t]) std::cout <<" "<< e;
  }
  std::cout << std::endl;
}


ThreadGroups ThreadGroups::filter (const IntVec& elmList) const
{
  ThreadGroups filtered;
  const ThreadGroups& group = *this;
  for (size_t i = 0; i < 2; ++i) {
    filtered[i].resize(group[i].size());
    for (size_t j = 0; j < group[i].size(); ++j)
      for (size_t k = 0; k < group[i][j].size(); ++k)
        if (std::find(elmList.begin(),
                      elmList.end(), group[i][j][k]) != elmList.end())
          filtered[i][j].push_back(group[i][j][k]);
  }

  return filtered;
}
