#include "LagrangeInterpolator.h"
#include <cstddef>

static double Lagrange(double x, size_t nr,
                       const std::vector<double>& grid)
{
  double result = 1.0;
  for (size_t j = 0; j< grid.size(); ++j) {
    if (j != nr)
      result *= (x-grid[j])/(grid[nr]-grid[j]);
  }

  return result;
}


double LagrangeInterpolator::evaluate(double x,
                                      const std::vector<double>& data)
{
  double result = 0.0;
  for (size_t i = 0; i < data.size(); ++i)
    result += data[i]*Lagrange(x, i, grid);

  return result;
}

Matrix LagrangeInterpolator::get(const std::vector<double>& new_grid)
{
  Matrix result(grid.size(), new_grid.size());
  for (size_t i = 0; i < grid.size(); ++i)
    for (size_t j = 0; j < new_grid.size(); ++j)
      result(i+1, j+1) =  Lagrange(new_grid[j], i, grid);

  return result;
}
