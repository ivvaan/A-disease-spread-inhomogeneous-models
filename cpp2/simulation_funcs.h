#pragma once
#include "lineages.h"
#include "trace.h"
#include "params.h"

void SIR_dynamics(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back);

/*
Less stochastic but more simple and fast implementation.
In the implementation below a person's spreadability exponetially declines in time.
While in the implementation above, THE EXPECTATION of
a person's spreadability exponetially declines in time.
*/
void SIR_dynamics2(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back);

void SIRS_dynamics(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back);

void SIRS_dynamics2(const Params& params, std::function<void(double, unsigned, unsigned)> emplace_back);

template <class TPopulation, class TLinages>
void SIR_dynamics_and_lineages_stat(const Params& params, TPopulation& population, TLinages& lineages,
  std::function<void(double, unsigned, unsigned)> emplace_back, std::vector<TL>& ttrace);

template <class TPopulation, class TLinages>
void SIRS_dynamics_and_lineages_stat(const Params& params, TPopulation& population, TLinages& lineages,
  std::function<void(double, unsigned, unsigned)> emplace_back, std::vector<TL>& ttrace);
