#pragma once
#include <iostream>
#include <functional>

#define FLD_INT(act,var_name,def_val) FLD_INT_##act(var_name,def_val)
#define FLD_DOUBLE(act,var_name,def_val) FLD_DOUBLE_##act(var_name,def_val)

#define PARAMS_LIST(act)\
FLD_INT(act,number_of_cities,10000) \
FLD_INT(act,population_size, 500000) \
FLD_INT(act,simulations_number, 50) \
FLD_INT(act,exec_type, 0) \
FLD_INT(act,is_country_level, 0) \
FLD_INT(act,additional_output, 0) \
FLD_INT(act,initial_seed, 317) \
FLD_INT(act,min_trace_len, 10) \
FLD_DOUBLE(act,final_time, 12) \
FLD_DOUBLE(act,resample_step, 0.003) \
FLD_DOUBLE(act,initial_R0, 5.2) \
FLD_DOUBLE(act,variance_of_social_activity, 0.5) \
FLD_DOUBLE(act,variance_of_catchability, 0.4) \
FLD_DOUBLE(act,variance_of_spreadability, 0.4) \
FLD_DOUBLE(act,inverse_tau, 2.5) \
FLD_DOUBLE(act, new_lineage_prob, 0.5) \
FLD_DOUBLE(act, test_prob, 0.05) \
FLD_DOUBLE(act, external_lineages_rate, 100) \
FLD_DOUBLE(act, immunity_loss_rate, 0.06) \
FLD_DOUBLE(act, immunity_loss_power, 2.) \
FLD_DOUBLE(act, protective_immunity_interval, 4.) \
FLD_DOUBLE(act,theta_city_attraction,0.1) \
FLD_DOUBLE(act,theta_city_R0,0.1) \
FLD_DOUBLE(act,theta_city_sa_variance,0.05) \
FLD_DOUBLE(act,zipf_exponent,1.07) \
FLD_DOUBLE(act, intercity_spread_probability, 1./128.) \
FLD_DOUBLE(act, nonneutral_mutation_prob, 1. / 8.)




#define FLD_INT_DECL(var_name,def_val) unsigned var_name = def_val;
#define FLD_DOUBLE_DECL(var_name,def_val) double var_name = def_val;
#define FLD_INT_READ(var_name,def_val) readers.push_back([&](const char* s) { sscanf_s(s, #var_name "%d", &var_name); });
#define FLD_DOUBLE_READ(var_name,def_val) readers.push_back([&](const char* s) { sscanf_s(s, #var_name "%lf", &var_name); });

struct Params {
  PARAMS_LIST(DECL)
  Params(std::istream& params_file)
  {
    std::vector<std::function<void(const char*)>> readers;
    PARAMS_LIST(READ);
    char s[1024];
    if (params_file)
      while (!params_file.eof()) {
        params_file.getline(s, 1024);
        for (auto& r : readers)
          r(s);
      };
  };

};

