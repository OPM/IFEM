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
#include "IFEM.h"

#include <algorithm>
#include <numeric>


namespace
{
  /*!
    \brief Assigns the knot spans in one parameter direction to tiles.
    \param[in] el Flags the non-zero knot spans
    \param[in] b Minimum number of non-zero knot spans in a tile
    \param[out] nTiles Number of tiles
    \return The tile index of each knot span

    \details The non-zero spans are spread evenly over as many tiles as have
    room for \a b of them, such that the tile widths differ by at most one.
    The zero spans stay with the tile of the non-zero span preceding them.
  */

  std::vector<int> tileIndex (const std::vector<bool>& el, int b, int& nTiles)
  {
    const int nzero = std::count(el.begin(),el.end(),true);
    nTiles = el.empty() ? 0 : std::max(1,nzero/std::max(b,1));

    std::vector<int> index(el.size(),0);
    for (size_t i = 0, k = 0; i < el.size(); i++)
      if (el[i])
        index[i] = (k++)*nTiles/nzero;
      else if (i > 0)
        index[i] = index[i-1];

    return index;
  }
}


void ThreadGroups::sequential (size_t nel, const IntVec& elms)
{
  tg.clear();
  if (elms.empty())
  {
    tg.resize(1,IntMat(1,IntVec(nel)));
    std::iota(tg[0][0].begin(),tg[0][0].end(),0);
  }
  else if (elms.front() >= 0)
    tg.resize(1,IntMat(1,elms));
}


void ThreadGroups::concurrent (size_t nel, const IntVec& elms)
{
  tg.clear();
  if (elms.empty())
  {
    tg.resize(1,IntMat(nel));
    for (size_t iel = 0; iel < nel; iel++)
      tg[0][iel].resize(1,iel);
  }
  else if (elms.front() >= 0)
  {
    tg.resize(1,IntMat(elms.size()));
    for (size_t i = 0; i < elms.size(); i++)
      tg[0][i].resize(1,elms[i]);
  }
}


void ThreadGroups::setColors (const IntMat& tasks, const IntVec& colors)
{
  tg.clear();
  for (size_t t = 0; t < tasks.size() && t < colors.size(); t++)
    if (!tasks[t].empty())
    {
      if (colors[t] >= static_cast<int>(tg.size()))
        tg.resize(colors[t]+1);
      tg[colors[t]].push_back(tasks[t]);
    }

  tg.erase(std::remove_if(tg.begin(),tg.end(),
                          [](const IntMat& color) { return color.empty(); }),
           tg.end());
}


std::vector<std::vector<int>> ThreadGroups::tiles (const BoolVec& el1,
                                                   const BoolVec& el2,
                                                   int b1, int b2,
                                                   IntVec* parity)
{
  return tiles(el1,el2,BoolVec(1,true),b1,b2,1,parity);
}


std::vector<std::vector<int>> ThreadGroups::tiles (const BoolVec& el1,
                                                   const BoolVec& el2,
                                                   const BoolVec& el3,
                                                   int b1, int b2, int b3,
                                                   IntVec* parity)
{
  int nt1, nt2, nt3;
  const IntVec t1 = tileIndex(el1,b1,nt1);
  const IntVec t2 = tileIndex(el2,b2,nt2);
  const IntVec t3 = tileIndex(el3,b3,nt3);

  IntMat result(nt1*nt2*nt3);
  for (size_t i3 = 0, iel = 0; i3 < el3.size(); i3++)
    for (size_t i2 = 0; i2 < el2.size(); i2++)
      for (size_t i1 = 0; i1 < el1.size(); i1++, iel++)
        result[t1[i1] + nt1*(t2[i2] + nt2*t3[i3])].push_back(iel);

  if (parity)
  {
    parity->resize(result.size());
    for (int k = 0, t = 0; k < nt3; k++)
      for (int j = 0; j < nt2; j++)
        for (int i = 0; i < nt1; i++, t++)
          (*parity)[t] = i%2 + 2*(j%2) + 4*(k%2);
  }

  return result;
}


void ThreadGroups::applyMap (const IntVec& map)
{
  for (IntMat& color : tg)
    for (IntVec& task : color)
      for (int& iel : task)
        iel = map[iel];
}


ThreadGroups ThreadGroups::filter (const IntVec& elmList) const
{
  IntVec sorted(elmList);
  std::sort(sorted.begin(),sorted.end());

  ThreadGroups filtered;
  for (const IntMat& color : tg)
  {
    IntMat tasks;
    for (const IntVec& task : color)
    {
      IntVec elms;
      for (int iel : task)
        if (std::binary_search(sorted.begin(),sorted.end(),iel))
          elms.push_back(iel);
      if (!elms.empty())
        tasks.push_back(elms);
    }
    if (!tasks.empty())
      filtered.tg.push_back(tasks);
  }

  return filtered;
}


void ThreadGroups::analyze (bool listAllSizes) const
{
  std::vector<size_t> sizes;
  sizes.reserve(tg.size());
  size_t nTasks = 0, nElms = 0;
  for (const IntMat& color : tg)
  {
    size_t n = 0;
    for (const IntVec& task : color)
      n += task.size();
    sizes.push_back(n);
    nTasks += color.size();
    nElms += n;
  }
  if (sizes.empty())
    return;

  const size_t num = sizes.size();
  const size_t min = *std::min_element(sizes.begin(),sizes.end());
  const size_t max = *std::max_element(sizes.begin(),sizes.end());
  std::vector<size_t> sorted(sizes);
  std::nth_element(sorted.begin(),sorted.begin()+num/2,sorted.end());

  IFEM::cout <<"\nElements are divided into "<< num <<" colors";
  if (nTasks < nElms)
    IFEM::cout <<" of "<< nTasks <<" tasks";
  IFEM::cout <<" (min = "<< min <<", max = "<< max
             <<", avg = "<< static_cast<double>(nElms) / static_cast<double>(num)
             <<", med = "<< sorted[num/2] <<").";
  if (listAllSizes)
    for (size_t i = 0; i < num; i++)
      IFEM::cout << (i%10 ? ' ' : '\n') << sizes[i];
  IFEM::cout << std::endl;
}
