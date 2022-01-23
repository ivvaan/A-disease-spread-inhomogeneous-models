#pragma once

#include "my_macros.h"
#include "params.h"
#include <algorithm>
#include <limits>

class ImmuneStatus {
  char* set=nullptr;
public:
  
  ImmuneStatus(unsigned int ps) {
    set = new char[ps];
    std::fill_n(set, ps, 0);
  };

  ~ImmuneStatus()
  {
      MY_DEL_ARR(set);
  };
  bool is_immune(unsigned int person) const { return set[person]; };
  void set_immune(unsigned int person) { set[person] = 1; };
  void set_susceptible(unsigned int person) { set[person] = 0; };
};


#define GET_UNSIGNED_DECREMENT(s) (unsigned)s
#define GET_SIGNED_DECREMENT(s) -(int)s
class SIRStatus {
public:
  static constexpr unsigned recovered = 0;
  static constexpr unsigned naive = 1;
  static constexpr unsigned cant_be_infected = 2;
  static constexpr unsigned incI = recovered;
  static constexpr unsigned decSincI = naive;
  static constexpr unsigned decI = 2;

  //static constexpr double end_of_time = std::numeric_limits<double>::max();
  struct ST {
    unsigned status = naive;
    double recovery_date;
  } *info=nullptr;

private:
  double immunity_loss_rate = 0.04;
  double immunity_loss_power = 2.0;
  double protective_immunity_interval = 4.;
public:
  SIRStatus(const Params& params):
  immunity_loss_rate(params.immunity_loss_rate),
  immunity_loss_power(params.immunity_loss_power),
  protective_immunity_interval(params.protective_immunity_interval)
  {
    info = new ST[params.population_size];
  };

  ~SIRStatus()
  {
    MY_DEL_ARR(info);
  };
  int try_to_infect(unsigned person, double T) const {
    auto& p_info = info[person];
    if (p_info.status == recovered) {
      auto dt = T - p_info.recovery_date - protective_immunity_interval;
      if (dt < 0)return cant_be_infected;
      if(
        //stduniform() < exp(-pow(dt * immunity_loss_rate, immunity_loss_power))
        stduniform() * (1.+pow(dt * immunity_loss_rate, immunity_loss_power)) < 1.0
        ) return cant_be_infected;
    }; 
    return p_info.status;
  };
  void set_infected(unsigned int person) {
    info[person].status = cant_be_infected; 
  };
  void set_recovered(unsigned int person,double T) { 
    info[person].status = recovered;
    info[person].recovery_date = T;
  };

  const ST& status_and_recovery_date(unsigned  person) { return info[person]; };
  //void set_susceptible(unsigned int person) { set[person] = 0; };
};

