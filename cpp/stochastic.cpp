// SIR.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <iostream>

#include "compartment.h"
#include "lineages.h"
#include "population.h"
#include "immune.h"
#include "random.h"
#include "trace.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <functional>

namespace fs = std::filesystem;

void get_lineages_stat(ClineageSet &lineages,
    double r0, double gamma, double T_stop,
    double external_rate, double fork_prop, double test_prop,
    CPopulation& population,
    std::vector<TSI>& dtrace, std::vector<TL>& ttrace)
{
    double T{ 0 };

    auto cs = population.candidate_to_infect();
    CImmuneStatus immune_status(population.get_size());
    immune_status.set_immune(cs);
    CCompartment active_infected(population.get_size());
    active_infected.add(cs);
    //CCompartment passive_infected(population.get_size());

    auto A = lineages.add_new(population.get_spreadability(cs), 0, cs, 0);
    int I_active = 1, I_passive = 0;
    int N = population.get_size();
    int S = N - I_active;

    dtrace.emplace_back(T, S, I_active + I_passive);
    r0 = (r0 + gamma - 1.0) / population.get_as0();
    double gamma_passive = gamma / (gamma - 1.);
    while (T < T_stop) {
        double spread = r0 * A + external_rate;
        double deactivate = gamma * I_active;
        double propensity = spread + deactivate + gamma_passive * I_passive;
        if (propensity == 0)
            break;
        T -= log(stduniform()) / propensity;
        double rnd = stduniform() * propensity;
        if (rnd < spread) {
            cs = population.candidate_to_infect();
            if (immune_status.is_immune(cs))
                continue;
            immune_status.set_immune(cs);
            active_infected.add(cs);
            if (rnd < external_rate)
                A += lineages.add_new(population.get_spreadability(cs), T, cs, 0);
            else 
                A += lineages.add(population.get_spreadability(cs), fork_prop * I_active, A, T, cs);
                //A += lineages.add(population.get_spreadability(cs), fork_prop, A, T, cs);

            ++I_active;
            --S;
            //testing
            if (stduniform() < test_prop)
                ttrace.emplace_back(T, lineages.personl(cs));

        } else {
            if ( rnd < spread+deactivate) {
                assert(I_active);
                cs = active_infected.remove();
                --I_active;
                ++I_passive;
                //passive_infected.add(cs);
                A -= lineages.deactivate(population.get_spreadability(cs), cs);
                if (I_active == 0)
                    A = 0;

            } else  {
                assert(I_passive);
                //cs = passive_infected.remove();
                //lineages.recover(cs);
                --I_passive;
            } 
        }
        dtrace.emplace_back(T, S, I_active + I_passive);
    }
    std::cout << "lineages " << lineages.lineage_numb();
    std::cout << "; ever infected " << population.get_size() - S << "; ";
};

void get_SIR_dynamics(double r0, double gamma, double T_stop, CPopulation& population, std::vector<TSI>& trace)
{
    double T{ 0 }, T_recover{ 0 };
    int I_active = 1, I_total = 1;
    int S = population.get_size() - I_active;

    auto cs = population.candidate_to_infect();
    CImmuneStatus immune_status(population.get_size());
    immune_status.set_immune(cs);
    CCompartment infected(population.get_size());
    infected.add(cs);

    auto A = population.get_spreadability(cs);
    trace.emplace_back(T, S, I_total);
    r0 = (r0 + gamma - 1.0) / population.get_as0();
    while (T < T_stop) {
        double spread = r0 * A;
        double propensity = spread + gamma * I_active;
        T -= log(stduniform()) / propensity;
        if (stduniform() * propensity < spread) {
            cs = population.candidate_to_infect();
            if (immune_status.is_immune(cs))
                continue;
            immune_status.set_immune(cs);
            infected.add(cs);
            A += population.get_spreadability(cs);
            ++I_active;
            ++I_total;
            --S;
        } else {
            if (0 == --I_active) {
                trace.emplace_back(T_stop, S, 0);
                return;
            }
            cs = infected.remove();
            A -= population.get_spreadability(cs);
            if (I_total > I_active)
                do
                    T_recover -= log(stduniform()) / I_total;
                while ((T_recover < T) && (I_active < --I_total));
            T_recover = T;
        }
        trace.emplace_back(T, S, I_total);
    }
};

/*
Less stochastic but more simple and fast implementation.
In the implementation below a person's spreadability exponetially declines in time.
While in the implementation above, THE EXPECTATION of 
a person's spreadability exponetially declines in time.
*/
void get_SIR_dynamics2(double r0, double gamma, double T_stop, CPopulation& population, std::vector<TSI>& trace)
{
    double T=0;
    int I = 1;
    int S = population.get_size() - 1;

    auto cs = population.candidate_to_infect();
    CImmuneStatus immune_status(population.get_size());
    immune_status.set_immune(cs);
    auto A = population.get_spreadability(cs);
    trace.emplace_back(T, S, I);
    double tau = 1.0 / gamma;
    double minus_r0tau = - tau * (r0 + gamma - 1.0) / population.get_as0();
    while (T < T_stop) {
        double delta = log(stduniform())/minus_r0tau;
        if (delta >= A) break;
        double dt=-tau*log(1.0-delta/A);
        A -= delta;  //speadability exponential decline: A-delta=A*exp(-dt/tau)
        T += dt;
        while (I) {
            //recovery events
            dt += log(stduniform()) / I;
            if (dt < 0)  break;
            --I;
        }
        cs = population.candidate_to_infect();
        if (!immune_status.is_immune(cs)){
            immune_status.set_immune(cs);
            A += population.get_spreadability(cs);
            ++I; --S;
        }
        trace.emplace_back(T, S, I);
    }
    while (I && (T < T_stop)){
        T -= log(stduniform()) / I--;
        trace.emplace_back(T, S, I);
    }
};

/*void get_SIR_dynamics3(double r0, double gamma, double T_stop, CPopulation& population, std::vector<TSI>& trace)
{
    double T{ 0 };
    int I_active = 1, I_passive = 0;
    int S = population.get_size() - I_active;

    auto cs = population.candidate_to_infect();
    CImmuneStatus immune_status(population.get_size());
    immune_status.set_immune(cs);
    CCompartment active_infected(population.get_size());
    active_infected.add(cs);

    auto A = population.get_spreadability(cs);
    trace.emplace_back(T, S, I_active + I_passive);
    r0 = (r0 + gamma - 1.0) / population.get_as0();
    double gamma_passive = gamma / (gamma - 1.);
    while (T < T_stop) {
        double spread = r0 * A;
        double deactivate = gamma * I_active;
        double propensity = spread + deactivate + gamma_passive * I_passive;
        if (propensity == 0)
            break;
        T -= log(stduniform()) / propensity;
        double rnd = stduniform() * propensity;
        if (rnd < spread) {
            cs = population.candidate_to_infect();
            if (immune_status.is_immune(cs))
                continue;
            immune_status.set_immune(cs);
            active_infected.add(cs);
            A += population.get_spreadability(cs);
            ++I_active;
            --S;
        } else {
            if (rnd < deactivate + spread) {
                assert(I_active);
                cs = active_infected.remove();
                --I_active;
                ++I_passive;
                A -= population.get_spreadability(cs);
                if (I_active == 0)
                    A = 0;

            } else  {
                assert(I_passive);
                --I_passive;
            }
        }
        trace.emplace_back(T, S, I_active + I_passive);
    }
}; */



void get_SIRS_dynamics(double r0, double gamma, double external_rate, double s_rate, double T_stop, CPopulation& population, std::vector<TSI>& trace)
{
    double T{ 0 };
    int I_active = 1, I_passive = 0, R=0;
    int S = population.get_size() - I_active;

    CCompartment active_infected(population.get_size());
    CCompartment passive_infected(population.get_size());
    CCompartment recovered(population.get_size());

    auto cs = population.candidate_to_infect();
    CImmuneStatus immune_status(population.get_size());
    immune_status.set_immune(cs);
    active_infected.add(cs);
    auto A = population.get_spreadability(cs);
    trace.emplace_back(T, S, I_active+I_passive);
    r0 = (r0 + gamma - 1.0) / population.get_as0();
    double gamma_passive = gamma / (gamma - 1.);
    while (T < T_stop) {
        double spread = r0 * A + external_rate;
        double deactivate = gamma * I_active;
        double recover = gamma_passive * I_passive;
        double become_s = s_rate * R;
        double propensity = spread + deactivate + recover + become_s;
        if (propensity == 0)
            break;
        T -= log(stduniform()) / propensity;
        double rnd = stduniform() * propensity;
        if (rnd < spread) {
            cs = population.candidate_to_infect();
            if (immune_status.is_immune(cs))
                continue;
            immune_status.set_immune(cs);
            active_infected.add(cs);
            A += population.get_spreadability(cs);
            ++I_active;
            --S;
        } else {
            if ((rnd -= spread) < deactivate) {
                assert(I_active);
                cs = active_infected.remove();
                --I_active;
                ++I_passive;
                passive_infected.add(cs);
                A -= population.get_spreadability(cs);
                if (I_active == 0)
                    A = 0;

            } else if ((rnd -= deactivate) < recover) {
                assert(I_passive);
                recovered.add(passive_infected.remove());
                --I_passive;
                ++R;
            } else  {
                immune_status.set_susceptible(recovered.remove());
                --R;
                ++S;
            } 
        }
        trace.emplace_back(T, S, I_active+I_passive);
    }
};


#define ADD_INT_PARAM(par_name, def_val) \
    unsigned int par_name = def_val;\
    readers.push_back([&](const char* s) { sscanf_s(s, #par_name "%d", &par_name); });

#define ADD_DOUBLE_PARAM(par_name, def_val) \
    double par_name = def_val;\
    readers.push_back([&](const char* s) { sscanf_s(s, #par_name "%lf", &par_name); });


const int max_exclude = 1000;

int main(int argc, char* argv[])
{
    std::vector<std::function<void(const char*)>> readers;
    ADD_DOUBLE_PARAM(final_time, 12);
    ADD_INT_PARAM(simulations_number, 50);
    ADD_INT_PARAM(exec_type, 0);
    //ADD_INT_PARAM(test_type, 0);
    ADD_INT_PARAM(initial_seed, 317);
    ADD_INT_PARAM(min_trace_len, 10);
    ADD_DOUBLE_PARAM(resample_step, 0.003);

    ADD_INT_PARAM(population_size, 100000);
    ADD_DOUBLE_PARAM(initial_R0, 5.2);
    ADD_DOUBLE_PARAM(variance_of_social_activity, 0.5);
    ADD_DOUBLE_PARAM(variance_of_catchability, 0.4);
    ADD_DOUBLE_PARAM(variance_of_spreadability, 0.4);
    ADD_DOUBLE_PARAM(inverse_tau, 2.5); // 1/tau
    ADD_DOUBLE_PARAM(new_lineage_prob, 0.25); // 1/tau
    ADD_DOUBLE_PARAM(test_prob, 0.25);
    ADD_DOUBLE_PARAM(external_lineages_rate, 100);
    ADD_DOUBLE_PARAM(s_rate, 0.06);

    const char* fname = "C:\\tmp";

    if (argc == 1) {
        printf("usage: SIR [-d[directory]]\n");
        //return 0;
    } else 
        for (int an = 1; an < argc; an++)
            if ((argv[an][0] == '-') && (argv[an][1] == 'd'))
                fname = argv[an] + 2;
   

      fs::path path(fname);

      std::ifstream conff(path / "config.txt");
      char s[1024];
      while (!conff.eof()) {
          conff.getline(s, 1024);
          for (auto r : readers)
              r(s);
      };

      set_seed(initial_seed);
      std::ofstream dynamicsf(path / "dynamics.txt");
      std::ofstream testf(path / "tests.txt");

      std::vector<TSI> dynamicst;
      dynamicst.reserve(70000);
      std::vector<TL> testt;
      testt.reserve(50000);
      for (unsigned int i = 0, N = 0; N < simulations_number; i++)
      {
          dynamicst.clear();
          testt.clear();
          CPopulation population{ population_size, variance_of_social_activity, 
            variance_of_catchability, variance_of_spreadability};
          ClineageSet lineages(population.get_size());
          switch (exec_type) {
          case 1:
              get_SIR_dynamics(initial_R0, inverse_tau, final_time, population, dynamicst);
              break;
          case 2:
              get_SIR_dynamics2(initial_R0, inverse_tau, final_time, population, dynamicst);
              break;
          case 3:
              get_SIRS_dynamics(initial_R0, inverse_tau, external_lineages_rate, s_rate, final_time, population, dynamicst);
              break;
          default:
              get_lineages_stat(lineages, initial_R0, inverse_tau, final_time,
                  external_lineages_rate, new_lineage_prob, test_prob, population,
                  dynamicst, testt);
              break;
          }
          std::cout << "length of #" << i + 1 << " trace=" << dynamicst.size();
          if (dynamicst.size() < min_trace_len) {
              std::cout << " excluded \n"; 
              if (i > N + max_exclude)
                  break;
              else
                  continue;
          }
          std::cout << "\n";

          print_resampled(dynamicst, resample_step, N, dynamicsf);
          if (exec_type==0)
              for (auto cur : testt)
                  cur.print_test(N,lineages, testf);
          N++;
      }
      return 0;
 }
 
