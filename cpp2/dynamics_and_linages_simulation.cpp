#include "simulation_funcs.h"

#include "compartment.h"
#include "lineages.h"
#include "population.h"
#include "person_status.h"
#include "random.h"
#include "trace.h"
#include "params.h"


#include <algorithm>
#include <cassert>
#include <vector>
#include <functional>


template <class TPopulation, class TLinages>
void SIR_dynamics_and_lineages_stat(const Params& params, TPopulation& population, TLinages& lineages,
  std::function<void(double, unsigned, unsigned)> emplace_back, std::vector<TL>& ttrace)
{
  auto r0 = params.initial_R0, gamma = params.inverse_tau,
    T_stop = params.final_time, external_rate = params.external_lineages_rate,
    fork_prop = params.new_lineage_prob, test_prop = params.test_prob;
  double T{ 0 };

  int N = population.get_size();
  ImmuneStatus immune_status(N);
  auto c_infected = population.candidate_to_infect(0);
  immune_status.set_immune(c_infected);
  InfectedCompartment active_infected(N);
  auto dA = population.get_spreadability(c_infected);
  active_infected.add(c_infected, lineages.add_new(dA, 0, c_infected, Lineages::external_lineage));

  unsigned I_active = 1, I_passive = 0;
  unsigned S = N - I_active;
  double gamma_passive = gamma / (gamma - 1.);

  while (T < T_stop) {
    /*dtrace.*/emplace_back(T, S, I_active + I_passive);
    double spread = external_rate + (I_active ? active_infected.A() : 0);
    double deactivate = gamma * I_active;
    double propensity = spread + deactivate + gamma_passive * I_passive;
    if (propensity == 0)
      break;
    T -= log(stduniform()) / propensity;
    double rnd = stduniform() * propensity;
    if (rnd < spread) {
      if (rnd < external_rate) {
        c_infected = population.candidate_to_infect();
        if (immune_status.is_immune(c_infected)) continue;
        dA = lineages.add_new(population.get_spreadability(c_infected), T, c_infected, Lineages::external_lineage);
      }
      else {
        assert(I_active > 0);
        auto cur_spreader = active_infected.select_spreader();
        c_infected = population.candidate_to_infect(cur_spreader);
        if (immune_status.is_immune(c_infected)) continue;
        dA = lineages.add(population.get_spreadability(c_infected), fork_prop * I_active, T, c_infected, lineages.personl(cur_spreader));
      }
      immune_status.set_immune(c_infected);
      active_infected.add(c_infected, dA);
      ++I_active;
      --S;
      //testing
      if (stduniform() < test_prop)
        ttrace.emplace_back(T, lineages.personl(c_infected));
      continue;
    }
    if (rnd < spread + deactivate) {
      assert(I_active);
      c_infected = active_infected.remove();
      --I_active;
      ++I_passive;
      lineages.deactivate(population.get_spreadability(c_infected), c_infected);
      continue;
    }
    assert(I_passive);
    --I_passive;
  }
  std::cout << "lineages " << lineages.lineage_numb();
  std::cout << "; ever infected " << population.get_size() - S << "; ";
};

/*template <class TPopulation, class TLinages>
void SIRS_dynamics_and_lineages_stat(const Params& params, TPopulation& population, TLinages& lineages,
*/  
template <class TPopulation, class TLinages>
void SIRS_dynamics_and_lineages_stat(const Params& params, TPopulation& population, TLinages& lineages,
  std::function<void(double, unsigned, unsigned)> emplace_back, std::vector<TL>& ttrace)
{
  auto r0 = params.initial_R0, gamma = params.inverse_tau,
    T_stop = params.final_time, external_rate = params.external_lineages_rate,
    fork_prop = params.new_lineage_prob, test_prop = params.test_prob;
  double T{ 0 };

  auto N = population.get_size();
  auto c_infected = population.candidate_to_infect(0);

  SIRStatus SIR_status(params);
  SIR_status.set_infected(c_infected);
  InfectedCompartment active_infected(N);
  Compartment passive_infected(N);
  auto dA = population.get_spreadability(c_infected);
  active_infected.add(c_infected, lineages.add_new(dA, 0, c_infected, Lineages::external_lineage));
  emplace_back(T, population.person_city(c_infected), SIRStatus::decSincI);

  unsigned I_active = 1, I_passive = 0;
  unsigned S = N - I_active;
  double gamma_passive = gamma / (gamma - 1.);

  while (T < T_stop) {
    /*dtrace.*/ //emplace_back(T, S, I_active + I_passive);
    double spread = external_rate + (I_active ? active_infected.A() : 0);
    double deactivate = gamma * I_active;
    double propensity = spread + deactivate + gamma_passive * I_passive;
    if (propensity == 0)
      break;
    T -= log(stduniform()) / propensity;
    double rnd = stduniform() * propensity;
    if (rnd < spread) {
      unsigned status;
      if (rnd < external_rate) {
        c_infected = population.candidate_to_infect();
        if ((status = lineages.try_to_infect(c_infected, Lineages::external_lineage,
          SIR_status.status_and_recovery_date(c_infected),T)) == SIRStatus::cant_be_infected) 
          continue;
        dA = lineages.add_new(population.get_spreadability(c_infected), T, c_infected, 0);
      }
      else {
        assert(I_active > 0);
        auto cur_spreader = active_infected.select_spreader();
        c_infected = population.candidate_to_infect(cur_spreader);
        auto spreader_lineage = lineages.personl(cur_spreader);
        if ((status = lineages.try_to_infect(c_infected, spreader_lineage,
          SIR_status.status_and_recovery_date(c_infected), T)) == SIRStatus::cant_be_infected)
          continue;
        dA = lineages.add(population.get_spreadability(c_infected), fork_prop * I_active, T, c_infected, spreader_lineage);
      }
      SIR_status.set_infected(c_infected);
      active_infected.add(c_infected, dA);
      ++I_active;
      S-= GET_UNSIGNED_DECREMENT(status);
      emplace_back(T, population.person_city(c_infected), status);

      //testing
      if (stduniform() < test_prop)
        ttrace.emplace_back(T, lineages.personl(c_infected));
      continue;
    }
    if (rnd < spread + deactivate) {
      assert(I_active);
      c_infected = active_infected.remove();
      --I_active; ++I_passive;
      passive_infected.add(c_infected);
      lineages.deactivate(population.get_spreadability(c_infected), c_infected);
      continue;
    }
    assert(I_passive);  
    c_infected = passive_infected.remove();
    emplace_back(T, population.person_city(c_infected), SIRStatus::decI);
    SIR_status.set_recovered(c_infected, T);
    --I_passive;
  }
  std::cout << "lineages " << lineages.lineage_numb();
  std::cout << "; ever infected " << population.get_size() - S << "; ";
};

template void SIR_dynamics_and_lineages_stat(const Params&, CountryPopulation&, LineageSet<TrivGenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIR_dynamics_and_lineages_stat(const Params&, Population&, LineageSet<TrivGenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIR_dynamics_and_lineages_stat(const Params&, CountryPopulation&, LineageSet<GenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIR_dynamics_and_lineages_stat(const Params&, Population&, LineageSet<GenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);

template void SIRS_dynamics_and_lineages_stat(const Params&, CountryPopulation&, LineageSet<TrivGenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIRS_dynamics_and_lineages_stat(const Params&, Population&, LineageSet<TrivGenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIRS_dynamics_and_lineages_stat(const Params&, CountryPopulation&, LineageSet<GenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
template void SIRS_dynamics_and_lineages_stat(const Params&, Population&, LineageSet<GenomeSet>&, std::function<void(double, unsigned, unsigned)>, std::vector<TL>&);
