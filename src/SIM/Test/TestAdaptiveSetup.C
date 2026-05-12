//==============================================================================
//!
//! \file TestAdaptiveSetup.C
//!
//! \date May 12 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for AdaptiveSetup.
//!
//==============================================================================

#include "AdaptiveSetup.h"
#include "SIMdummy.h"
#include "SIMgeneric.h"

#include "Catch2Support.h"

#include <algorithm>
#include <numeric>
#include <string>

namespace {

class SIMA : public SIMdummy<SIMgeneric>
{
public:
  SIMA() {}
};

class TestAdaptiveSetupSimHolder
{
protected:
  SIMA simDummy;
};

class TestAdaptiveSetup : private TestAdaptiveSetupSimHolder,
                          public AdaptiveSetup
{
public:
  explicit TestAdaptiveSetup(size_t N)
    : AdaptiveSetup(simDummy)
  {
    for (size_t i = 0; i < N; ++i)
      errors.emplace_back(i+1.0, i);
    std::reverse(errors.begin(), errors.end());
  }

  double sum() const
  {
    return std::accumulate(errors.begin(), errors.end(), 0.0,
                           [](const double acc, const DblIdx& a)
                           { return acc + a.first; });
  }

  size_t size() const { return errors.size(); }

  using AdaptiveSetup::ErrorInfo;
  ErrorInfo errorInfo(const double beta,
                      const Threshold threshold)
  {
    return this->errorLimits(errors, beta, threshold);
  }

  std::pair<double,double>
  dorfelLimits(const double limit)
  {
    double currErr = 0.0;
    size_t i = 0;
    while (currErr < limit && i < errors.size()) {
      currErr += errors[i].first;
      ++i;
    }

    return i == errors.size() ? std::make_pair(errors.back().first, currErr)
                              : std::make_pair(errors[i].first, currErr);
  }

  using Threshold = AdaptiveSetup::Threshold;

private:
  std::vector<DblIdx> errors;
};

}


TEST_CASE("TestAdaptiveSetup.ErrorLimits")
{
  const size_t N = GENERATE(10, 13, 25);

  SECTION("N = " + std::to_string(N))
  {
    TestAdaptiveSetup setup(N);
    const double sum = setup.sum();
    const size_t size = setup.size();
    const double minErr = 1.0;

    const double beta = GENERATE(1.0, 10.0, 50.0, 101.0);

    SECTION("MAXIMUM, beta = " + std::to_string(beta))
    {
      const TestAdaptiveSetup::ErrorInfo info =
        setup.errorInfo(beta, TestAdaptiveSetup::Threshold::MAXIMUM);
      REQUIRE_THAT(info.limit, WithinRel(size*beta*0.01));
    }

    SECTION("AVERAGE, beta = " + std::to_string(beta))
    {
      const TestAdaptiveSetup::ErrorInfo info =
        setup.errorInfo(beta, TestAdaptiveSetup::Threshold::AVERAGE);
      REQUIRE_THAT(info.sumErr, WithinRel(sum));
      REQUIRE_THAT(info.limit, WithinRel(sum / size  * beta * 0.01));
    }

    SECTION("MINIMUM, beta = " + std::to_string(beta))
    {
      const TestAdaptiveSetup::ErrorInfo info =
        setup.errorInfo(beta, TestAdaptiveSetup::Threshold::MINIMUM);
      REQUIRE_THAT(info.limit, WithinRel(minErr * beta * 0.01));
    }

    SECTION("DÖRFEL, beta = " + std::to_string(beta))
    {
      const TestAdaptiveSetup::ErrorInfo info =
        setup.errorInfo(beta, TestAdaptiveSetup::Threshold::DORFEL);
      const auto& [limit, curr] = setup.dorfelLimits(sum * beta * 0.01);
      REQUIRE_THAT(info.limit, WithinRel(limit));
      REQUIRE_THAT(info.sumErr, WithinRel(sum * beta * 0.01));
      REQUIRE_THAT(info.curErr, WithinRel(curr));
    }
  }
}
