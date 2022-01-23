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

void SIR_dynamics(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back)
{
  Population population{ params };
  auto gamma = params.inverse_tau, T_stop = params.final_time;


  double T{ 0 }, T_recover{ 0 };
  int I_active = 1, I_total = 1;
  int S = population.get_size() - I_active;

  auto c_infected = population.candidate_to_infect();
  ImmuneStatus immune_status(population.get_size());
  immune_status.set_immune(c_infected);
  Compartment infected(population.get_size());
  infected.add(c_infected);

  auto spread = population.get_spreadability(c_infected);
  while (T < T_stop) {
    /*dtrace.*/emplace_back(T, S, I_total);
    //double spread = A;
    double propensity = spread + gamma * I_active;
    T -= log(stduniform()) / propensity;
    if (stduniform() * propensity < spread) {
      c_infected = population.candidate_to_infect();
      if (immune_status.is_immune(c_infected))
        continue;
      immune_status.set_immune(c_infected);
      infected.add(c_infected);
      spread += population.get_spreadability(c_infected);
      ++I_active;
      ++I_total;
      --S;
      continue;
    }
    if (0 == --I_active)return; 
    spread -= population.get_spreadability(infected.remove());
    if (I_total > I_active)
      do
        T_recover -= log(stduniform()) / I_total;
      while ((T_recover < T) && (I_active < --I_total));
    T_recover = T;
  }
};

/*
Less stochastic but more simple and fast implementation.
In the implementation below a person's spreadability exponetially declines in time.
While in the implementation above, THE EXPECTATION of
a person's spreadability exponetially declines in time.
*/
void SIR_dynamics2(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back)
{
  Population population{ params };
  auto gamma = params.inverse_tau, T_stop = params.final_time;

  double T = 0;
  int I = 1;
  int S = population.get_size() - 1;

  auto c_infected = population.candidate_to_infect();
  ImmuneStatus immune_status(population.get_size());
  immune_status.set_immune(c_infected);
  auto A = population.get_spreadability(c_infected);
  double tau = 1.0 / gamma;
  while (T < T_stop) {
    /*dtrace.*/emplace_back(T, S, I);
    double delta = - log(stduniform()) * gamma;
    if (delta >= A) break;
    double dt = - tau * log(1.0 - delta / A);
    A -= delta;  //speadability exponential decline: A-delta=A*exp(-dt/tau)
    T += dt;
    while (I) {
      //recovery events
      dt += log(stduniform()) / I;
      if (dt < 0)  break;
      --I;
    }
    c_infected = population.candidate_to_infect();
    if (!immune_status.is_immune(c_infected)) {
      immune_status.set_immune(c_infected);
      A += population.get_spreadability(c_infected);
      ++I; --S;
    }
  }
  while (I && (T < T_stop)) {
    T -= log(stduniform()) / I--;
    /*dtrace.*/emplace_back(T, S, I);
  }
};

void SIRS_dynamics(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back)
{
  Population population{ params };
  auto  gamma = params.inverse_tau, T_stop = params.final_time,
    external_rate = params.external_lineages_rate, immunity_loss_rate = params.immunity_loss_rate;

  double T{ 0 };
  int I_active = 1, I_passive = 0, R = 0;
  int S = population.get_size() - I_active;

  Compartment active_infected(population.get_size());
  Compartment passive_infected(population.get_size());
  Compartment recovered(population.get_size());

  auto c_infected = population.candidate_to_infect();
  ImmuneStatus immune_status(population.get_size());
  immune_status.set_immune(c_infected);
  active_infected.add(c_infected);
  auto A = population.get_spreadability(c_infected);
  double gamma_passive = gamma / (gamma - 1.);
  while (T < T_stop) {
    /*dtrace.*/emplace_back(T, S, I_active + I_passive);
    double spread = A + external_rate;
    double deactivate = gamma * I_active;
    double recover = gamma_passive * I_passive;
    double become_s = immunity_loss_rate * R;
    double propensity = spread + deactivate + recover + become_s;
    if (propensity == 0)
      break;
    T -= log(stduniform()) / propensity;
    double rnd = stduniform() * propensity;
    if (rnd < spread) {
      c_infected = population.candidate_to_infect();
      if (immune_status.is_immune(c_infected))
        continue;
      immune_status.set_immune(c_infected);
      active_infected.add(c_infected);
      A += population.get_spreadability(c_infected);
      ++I_active;
      --S;
      continue;
    }
    if ((rnd -= spread) < deactivate) {
      assert(I_active);
      c_infected = active_infected.remove();
      --I_active;
      ++I_passive;
      passive_infected.add(c_infected);
      A -= population.get_spreadability(c_infected);
      if (I_active == 0)
        A = 0;
      continue;
    }
    if ((rnd -= deactivate) < recover) {
      assert(I_passive);
      recovered.add(passive_infected.remove());
      --I_passive;
      ++R;
      continue;
    }
    immune_status.set_susceptible(recovered.remove());
    --R;
    ++S;
  }
};

void SIRS_dynamics2(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back)
{
  Population population{ params };
  auto gamma = params.inverse_tau, T_stop = params.final_time,
    external_rate = params.external_lineages_rate, immunity_loss_rate = params.immunity_loss_rate;

  double T{ 0 };
  int I_active = 1, I_passive = 0;
  int S = population.get_size() - I_active;

  Compartment active_infected(population.get_size());
  Compartment passive_infected(population.get_size());

  auto c_infected = population.candidate_to_infect();
  SIRStatus SIR_status(params);
  SIR_status.set_infected(c_infected);
  active_infected.add(c_infected);
  auto A = population.get_spreadability(c_infected);
  double gamma_passive = gamma / (gamma - 1.);
  while (T < T_stop) {
    /*dtrace.*/emplace_back(T, S, I_active + I_passive);
    double spread = A + external_rate;
    double deactivate = gamma * I_active;
    double recover = gamma_passive * I_passive;
    //double become_s = immunity_loss_rate * R;
    double propensity = spread + deactivate + recover;// +become_s;
    if (propensity == 0)
      break;
    T -= log(stduniform()) / propensity;
    double rnd = stduniform() * propensity;
    if (rnd < spread) {
      c_infected = population.candidate_to_infect();
      unsigned status;
      if ((status = SIR_status.try_to_infect(c_infected, T)) == SIRStatus::cant_be_infected)
        continue;
        
      SIR_status.set_infected(c_infected);
      active_infected.add(c_infected);
      A += population.get_spreadability(c_infected);
      ++I_active;
      S-= GET_UNSIGNED_DECREMENT(status);
      continue;
    }
    if (rnd< deactivate+spread) {
      assert(I_active);
      c_infected = active_infected.remove();
      --I_active;
      ++I_passive;
      passive_infected.add(c_infected);
      A -= population.get_spreadability(c_infected);
      if (I_active == 0)
        A = 0;
      continue;
    }
    assert(I_passive);
    SIR_status.set_recovered(passive_infected.remove(),T);
    --I_passive;
  }
};
