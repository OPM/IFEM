#ifndef LAGRANGE_INTERPOLATOR_H_
#define LAGRANGE_INTERPOLATOR_H_

#include "MatVec.h"
#include <vector>

class LagrangeInterpolator {
  public:
    explicit LagrangeInterpolator(const std::vector<double>& grid_) :
      grid(grid_)
    {
    }

    double evaluate(double x, const std::vector<double>& data);

    Matrix get(const std::vector<double>& new_grid);

  protected:
    std::vector<double> grid; //!< Grid points.
};

#endif
