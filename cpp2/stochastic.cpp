// SIR.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include "compartment.h"
#include "lineages.h"
#include "population.h"
#include "person_status.h"
#include "random.h"
#include "trace.h"
#include "params.h"
#include "simulation_funcs.h"

#include <algorithm>
#include <iostream>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <functional>

namespace fs = std::filesystem;

const int max_exclude = 1000;

int main(int argc, char* argv[])
{
    const char* fname = "C:\\tmp";

    if (argc == 1) {
        printf("usage: stochastic [-d[directory]]\n");
        //return 0;
    } else 
        for (int an = 1; an < argc; an++)
            if ((argv[an][0] == '-') && (argv[an][1] == 'd'))
                fname = argv[an] + 2;
    fs::path path(fname);
    std::ifstream params_file(path / "config.txt");
    Params params(params_file);
    set_seed(params.initial_seed);
    std::ofstream dynamicsf(path / "dynamics.txt");
    std::ofstream testf(path / "tests.txt");
    std::ofstream qdynamicsf;
    if (params.additional_output == 1)
      qdynamicsf.open(path / "qdynamics.txt");
    std::ofstream nri_dynamicsf;
    if (params.additional_output == 2)
      nri_dynamicsf.open(path / "nri_dynamics.txt");

    std::vector<TSI> dynamicst;
    dynamicst.reserve(40000);
    std::vector<TSIQ> qdynamicst;
    qdynamicst.reserve(40000);    std::vector<TNRI> nri_dynamicst;
    nri_dynamicst.reserve(40000);

    std::vector<TL> testt;
    testt.reserve(50000);

    for (unsigned int i = 0, N = 0; N < params.simulations_number; i++)
    {
      dynamicst.clear();
      qdynamicst.clear();
      testt.clear();
      nri_dynamicst.clear();
          
      auto exec_type = params.exec_type;
      if (exec_type == 0)
        exec_type = 3 + params.is_country_level;
        switch (exec_type) {
          case 1: {
            SIRS_dynamics2(params, get_resampler(dynamicst, params.resample_step));
          }
              break;
          case 2: {
            CountryPopulation population(params);
            LineageSet<TrivGenomeSet> lineages(params);
            SIR_dynamics_and_lineages_stat(params, 
              population, lineages, get_resampler(dynamicst, params.resample_step), testt);
              for (auto cur : testt)
                cur.print_test(N, lineages, testf);
          }
              break;
          case 3: {
            Population population(params);
            LineageSet<GenomeSet> lineages(params);
            auto resampler = get_resampler_ex(params,population,
              dynamicst, qdynamicst,nri_dynamicst);
            SIRS_dynamics_and_lineages_stat(params,
              population, lineages,
              resampler, testt);
            for (auto cur : testt)
              cur.print_test(N, lineages, testf);
          }
                break;
          case 4: {
            CountryPopulation population(params);
            LineageSet<GenomeSet> lineages(params);
            auto resampler = get_resampler_ex(params, population,
              dynamicst, qdynamicst, nri_dynamicst);
            SIRS_dynamics_and_lineages_stat(params,
              population, lineages,
              resampler, testt);
            for (auto cur : testt)
              cur.print_test(N, lineages, testf);
          }
                break;
          default: {};
            break;
        }

        std::cout << "length of #" << i + 1 << " trace=" << dynamicst.size();
        if (dynamicst.size() < params.min_trace_len) {
            std::cout << " excluded \n"; 
            if (i > N + max_exclude)
                break;
            else
                continue;
        }
        std::cout << "\n";
        for(auto &cur: dynamicst)
          dynamicsf << N << " " << cur << "\n";
        if (params.additional_output == 1)
          for (auto& cur : qdynamicst)
            qdynamicsf << N << " " << cur << "\n";
        if (params.additional_output == 2)
          for (auto& cur : nri_dynamicst)
            nri_dynamicsf << N << " " << cur << "\n";
        //print_resampled(dynamicst, params.resample_step, N, dynamicsf);
        N++;
    }
    return 0;
 }
 
