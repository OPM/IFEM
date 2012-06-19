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
#if SP_DEBUG > 1
#include <iostream>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif


void ThreadGroups::calcGroups (const BoolVec& el1, const BoolVec& el2,
                               int p1, int p2)
{
  // Count the non-zero element in each direction, the zero-span elements
  // should not affect the partitioning as they don't involve any work
  size_t i, nel1 = 0, nel2 = 0;
  for (i = 0; i < el1.size(); i++)
    if (el1[i]) nel1++;
  for (i = 0; i < el2.size(); i++)
    if (el2[i]) nel2++;

  int threads=1;
  int stripsize=0;
  int remainder=0;
  int dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  dir = getStripDirection(nel1,nel2,parts);
  mul = dir == 0 ? 1 : el1.size();
  int els = dir == 0 ? nel1 : nel2;

  // The minimum strip size (with) depends on the polynomial degree
  // due to the overlapping support of the splines basis functions
  int minsize = dir == 0 ? p1 : p2;
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

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
#endif

  nel1 = el1.size();
  nel2 = el2.size();
  if (threads == 1)
  {
    tg[0].resize(1);
    tg[0][0].reserve(nel1*nel2);
    for (i = 0; i < nel1*nel2; ++i)
      tg[0][0].push_back(i);
    tg[1].clear();
  }
  else
  {
    int i, j, t, zspan, offs = 0;
    IntVec stripsizes[2], startelms[2];
    const BoolVec& elz = dir == 0 ? el1 : el2;

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
        int maxx = dir == 0 ? stripsizes[i][t] : nel1;
        int maxy = dir == 1 ? stripsizes[i][t] : nel2;
        tg[i][t].reserve(maxx*maxy);
        for (int i2 = 0; i2 < maxy; ++i2)
          for (int i1 = 0; i1 < maxx; ++i1)
            tg[i][t].push_back(startelms[i][t]+i1+i2*nel1);
      }
#if defined(USE_OPENMP) && SP_DEBUG > 1
      std::cout << "group " << i << std::endl;
      for (size_t j = 0; j < tg[i].size(); ++j) {
        std::cout << "\t thread " << j << ": ";
        for (size_t k = 0; k < tg[i][j].size(); ++k)
          std::cout << tg[i][j][k] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }
}


void ThreadGroups::calcGroups (int nel1, int nel2, int minsize)
{
  int threads=1;
  int stripsize=0;
  int remainder=0;
  int dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  dir = getStripDirection(nel1,nel2,parts);
  mul = dir == 0 ? 1 : nel1;
  int els = dir == 0 ? nel1 : nel2;

  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

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
#endif

  if (threads == 1)
  {
    tg[0].resize(1);
    tg[0][0].reserve(nel1*nel2);
    for (int i = 0; i < nel1*nel2; ++i)
      tg[0][0].push_back(i);
    tg[1].clear();
  }
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
        int maxx = dir == 0 ? stripsizes[i][t] : nel1;
        int maxy = dir == 1 ? stripsizes[i][t] : nel2;
        for (int i2 = 0; i2 < maxy; ++i2)
          for (int i1 = 0; i1 < maxx; ++i1)
            tg[i][t].push_back(startelms[i][t]+i1+i2*nel1);
      }
#if defined(USE_OPENMP) && SP_DEBUG > 1
      std::cout << "group " << i << std::endl;
      for (size_t j = 0; j < tg[i].size(); ++j) {
        std::cout << "\t thread " << j << ": ";
        for (size_t k = 0; k < tg[i][j].size(); ++k)
          std::cout << tg[i][j][k] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }
}


int ThreadGroups::getStripDirection (int nel1, int nel2, int parts)
{
  int s1 = nel1 / parts;
  int s2 = nel2 / parts;
  int r1 = nel1 - s1*parts;
  int r2 = nel2 - s2*parts;

  if (r1*nel2 < r2*nel1)
    return 0; // strips in u-direction
  else if (r1*nel2 > r2*nel1)
    return 1; // strips in v-direction

  if (nel1 > nel2)
    return 0; // strips in u-direction
  else
    return 1; // strips in v-direction
}


void ThreadGroups::calcGroups (const BoolVec& el1, const BoolVec& el2,
                               const BoolVec& el3, int p1, int p2, int p3)
{
  // Count the non-zero element in each direction, the zero-span elements
  // should not affect the partitioning as they don't involve any work
  size_t i, nel1 = 0, nel2 = 0, nel3 = 0;
  for (i = 0; i < el1.size(); i++)
    if (el1[i]) nel1++;
  for (i = 0; i < el2.size(); i++)
    if (el2[i]) nel2++;
  for (i = 0; i < el3.size(); i++)
    if (el3[i]) nel3++;

  int threads=1;
  int stripsize=0;
  int remainder=0;
  int dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  dir = getStripDirection(nel1,nel2,nel3,parts);
  mul = dir == 0 ? 1 : el1.size()*(dir == 1 ? 1 : el2.size());
  int els = dir == 0 ? nel1 : (dir == 1 ? nel2 : nel3);

  // The minimum strip size (with) depends on the polynomial degree
  // due to the overlapping support of the splines basis functions
  int minsize = dir == 0 ? p1 : (dir == 1 ? p2 : p3);
  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

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
#endif

  nel1 = el1.size();
  nel2 = el2.size();
  nel3 = el3.size();
  if (threads == 1)
  {
    tg[0].resize(1);
    tg[0].reserve(nel1*nel2*nel3);
    for (i = 0; i < nel1*nel2*nel3; ++i)
      tg[0][0].push_back(i);
    tg[1].clear();
  }
  else
  {
    int i, j, t, zspan, offs = 0;
    IntVec stripsizes[2], startelms[2];
    const BoolVec& elz = dir == 0 ? el1 : (dir == 1 ? el2 : el3);

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
        int maxx = dir == 0 ? stripsizes[i][t] : nel1;
        int maxy = dir == 1 ? stripsizes[i][t] : nel2;
        int maxz = dir == 2 ? stripsizes[i][t] : nel3;
        tg[i][t].reserve(maxx*maxy*maxz);
        for (int i3 = 0; i3 < maxz; ++i3)
          for (int i2 = 0; i2 < maxy; ++i2)
            for (int i1 = 0; i1 < maxx; ++i1)
              tg[i][t].push_back(startelms[i][t]+i1+nel1*(i2+nel2*i3));
      }

#if defined(USE_OPENMP) && SP_DEBUG > 1
      std::cout << "group " << i << std::endl;
      for (size_t j = 0; j < tg[i].size(); ++j) {
        std::cout << "\t thread " << j << " (" << tg[i][j].size() << "): ";
        for (size_t k = 0; k < tg[i][j].size(); ++k)
          std::cout << tg[i][j][k] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }
}


void ThreadGroups::calcGroups (int nel1, int nel2, int nel3, int minsize)
{
  int threads=1;
  int stripsize=0;
  int remainder=0;
  int i, dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  dir = getStripDirection(nel1,nel2,nel3,parts);
  mul = dir == 0 ? 1 : nel1*(dir == 1 ? 1 : nel2);
  int els = dir == 0 ? nel1 : (dir == 1 ? nel2 : nel3);

  while (threads > 1 && (stripsize = els/parts) < minsize) {
    threads --;
    parts -= 2;
  }

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
#endif

  if (threads == 1)
  {
    tg[0].resize(1);
    tg[0].reserve(nel1*nel2*nel3);
    for (i = 0; i < nel1*nel2*nel3; ++i)
      tg[0][0].push_back(i);
    tg[1].clear();
  }
  else
  {
    IntVec stripsizes[2];
    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    IntVec startelms[2];
    for (int offs = i = 0; i < threads; ++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }

    for (i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = dir == 0 ? stripsizes[i][t] : nel1;
        int maxy = dir == 1 ? stripsizes[i][t] : nel2;
        int maxz = dir == 2 ? stripsizes[i][t] : nel3;
        for (int i3 = 0; i3 < maxz; ++i3)
          for (int i2 = 0; i2 < maxy; ++i2)
            for (int i1 = 0; i1 < maxx; ++i1)
              tg[i][t].push_back(startelms[i][t]+i1+nel1*(i2+nel2*i3));
      }

#if defined(USE_OPENMP) && SP_DEBUG > 1
      std::cout << "group " << i << std::endl;
      for (size_t j = 0; j < tg[i].size(); ++j) {
        std::cout << "\t thread " << j << " (" << tg[i][j].size() << "): ";
        for (size_t k = 0; k < tg[i][j].size(); ++k)
          std::cout << tg[i][j][k] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }
}


int ThreadGroups::getStripDirection (int nel1, int nel2, int nel3, int parts)
{
  int s1 = nel1 / parts;
  int s2 = nel2 / parts;
  int s3 = nel3 / parts;
  int r1 = nel1 - s1*parts;
  int r2 = nel2 - s2*parts;
  int r3 = nel3 - s3*parts;

  if (r1*nel2*nel3 < nel1*r2*nel3 && r1*nel2*nel3 < nel1*nel2*r3)
    return 0; // strips along x axis
  else if (nel1*r2*nel3 < r1*nel2*nel3 && nel1*r2*nel3 < nel1*nel2*r3)
    return 1; // strips along y axis
  else if (nel1*nel2*r3 < r1*nel2*nel3 && nel1*nel2*r3 < nel1*r2*nel3)
    return 2; // strips along z axis

  // The number of left-over elements is not smallest in one direction only
  if (r1*nel2*nel3 > nel1*r2*nel3)
    return nel2 > nel3 ? 1 : 2;
  else if (nel1*r2*nel3 > nel1*nel2*r3)
    return nel1 > nel3 ? 0 : 2;
  else if (nel1*nel2*r3 > r1*nel2*nel3)
    return nel1 > nel2 ? 0 : 1;

  // The number of left-over elements is the same in all three directions
  if (nel1 >= nel2 && nel1 >= nel3)
    return 0;
  else if (nel2 >= nel1 && nel2 >= nel3)
    return 1;
  else
    return 2;
}


void ThreadGroups::applyMap (const IntVec& map)
{
  for (size_t l = 0; l < 2; ++l)
    for (size_t k = 0; k < tg[l].size(); ++k)
      for (size_t j = 0; j < tg[l][k].size(); ++j)
        tg[l][k][j] = map[tg[l][k][j]];
}
