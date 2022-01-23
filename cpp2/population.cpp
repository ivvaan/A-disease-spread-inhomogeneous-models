#include "population.h"
#include "random.h"
#include "my_macros.h"
#include "utils.h"
#include <math.h>
#include <numeric> 
#include <cassert>

double* get_persons_distr(unsigned ps, double v)
{
    auto res = new double[ps];
    auto last = res + ps;
    auto m = -0.5 * v;
    double sv = sqrt(v);
    for (auto p = res; p < last; ++p)
        *p = exp(stdnormal() * sv + m);
    return res;
};

Population::Population(const Params &params)
    : pop_size(params.population_size)
    , vs(params.variance_of_social_activity)
{

    catchability = get_persons_distr(params.population_size, params.variance_of_catchability);
    spreadability = get_persons_distr(params.population_size, params.variance_of_spreadability);
    auto corrected_r0_log = log((params.initial_R0 + params.inverse_tau - 1.)/get_as0(params.variance_of_social_activity));

    auto get_sa = [m = corrected_r0_log - 0.5 * params.variance_of_social_activity, sv = sqrt(params.variance_of_social_activity)]() { return exp(stdnormal() * sv + m); };
    auto psa = get_sa();
    spreadability[0] *= psa;
    psa = (catchability[0] *= psa);
    for (unsigned int i = 1; i < params.population_size; ++i) {
        auto csa = get_sa();
        spreadability[i] *= csa;
        catchability[i] *= csa;
        psa = (catchability[i] += psa);
    }
    ca_sum = psa;
}

Population::~Population()
{
    MY_DEL_ARR(catchability);
    MY_DEL_ARR(spreadability);
};
unsigned int Population::candidate_to_infect() const
{
  int m, l = 0, r = pop_size - 1;
  double s = ca_sum * stduniform();
  while (r > l) {
    m = (l + r) / 2;
    if (s < catchability[m])
      r = m;
    else
      l = m + 1;
  }
  return l;
};


unsigned int Population::candidate_to_infect(unsigned f, unsigned t) const
{
  assert(t != 0);
  double cf = f ? catchability[f - 1] : 0;
  double s = cf+ stduniform()*(catchability[t-1]- cf);
  int m, l = f, r = t - 1;
  
  while (r > l) {
    m = (l + r) / 2;
    if (s < catchability[m])
      r = m;
    else
      l = m + 1;
  }
  return l;
};


CountryPopulation::CountryPopulation(const Params& params)
  :Population(),cities(params)
{
  pop_size = params.population_size;
  vs = params.variance_of_social_activity;
  catchability = get_persons_distr(params.population_size, params.variance_of_catchability);
  spreadability = get_persons_distr(params.population_size, params.variance_of_spreadability);
  cities.multiply_by_social_activity(params,catchability, spreadability);
  std::partial_sum(catchability, catchability + pop_size, catchability);
  ca_sum = catchability[params.population_size-1];
  cities.init_intercity_prob(params.intercity_spread_probability,catchability);
}

unsigned int CountryPopulation::candidate_to_infect(unsigned spreader) const
{
  auto spreader_city = cities.personc(spreader);
  if (stduniform() < cities.p_intercity_spread(spreader_city))
    return Population::candidate_to_infect();
  return Population::candidate_to_infect(cities.begin(spreader_city),cities.end(spreader_city));
};
