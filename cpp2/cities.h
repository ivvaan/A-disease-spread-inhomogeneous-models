#pragma once
#include "random.h"
#include "my_macros.h"
#include "params.h"
#include "utils.h"


#include <math.h>
#include <algorithm>

/*
class GroupStatBase {

public:
  enum class Events { infect, deactivate, cure, reinfect };
  void count(double T, unsigned city, Events action) {};
};

template<class GroupStatIpml>
class GroupStat:public GroupStatBase {
protected:
  GroupStatIpml* parent = nullptr;
public:
  void count(double T,unsigned city, Events action) { 
    if (parent)parent->count(T, city, action); 
  };
  void set_parent(GroupStatIpml* p) { parent = p; };
};

class GroupStatIpml :public GroupStat<GroupStatIpml>{
  unsigned N = 0, S = 0, I = 0;
public:
  void count(double T, unsigned city, Events action) {
    if (parent)parent->count(T,city, action);
    if (action == Events::infect) { --S; ++I; return; }
    if (action == Events::cure) { --I; return; }
    if (action == Events::reinfect) { ++I; return; }
  };

};
*/

struct City {
public:
  unsigned int ID = 0;
  unsigned int begin_from = 0;
  unsigned int size = 0;
  unsigned int I_active = 0;
  double expected_size = 1;
  double A = 0;
  double intercity_prob = 1.;
};



class CitySet {

public:
  static constexpr unsigned big_cities = 0;
  static constexpr unsigned mid_cities = 1;
  static constexpr unsigned small_cities = 2;
  static constexpr unsigned tiny_cities = 3;

  struct QuartileBounds {
    unsigned big, mid, small;
    QuartileBounds():
      big(0), mid(0), small(0) {};
    QuartileBounds(unsigned b, unsigned m, unsigned s):
    big(b),mid(m),small(s){};
    unsigned classify(unsigned city)const {
      return (city < mid) ? (city < big ? big_cities : mid_cities) :
        (city < small ? small_cities : tiny_cities);
    }
  };



  unsigned N=0;
  unsigned total_population=0;
  City* cities=nullptr;
  unsigned int* person_city=nullptr;
  double  theta_i=0.1, theta_r=0.1, theta_v=0.05;
  double c_v=1.;
  CitySet(const Params& params)
    : N{ params.number_of_cities },
    total_population(params.population_size),
    theta_i{ params.theta_city_attraction }, 
    theta_r{ params.theta_city_R0 }, 
    theta_v{ params.theta_city_sa_variance} 
  {
    cities = new City[params.number_of_cities];

    double zipf_coef = 1;
    for (unsigned i = 1; i < params.number_of_cities; i++)zipf_coef += pow(1.0 + i, -params.zipf_exponent);
    zipf_coef = (double)params.population_size / zipf_coef;
    {
      DynamicCategiricalDistr<2> distr(params.number_of_cities);
      for (unsigned i = 0; i < params.number_of_cities; i++) {
        cities[i].ID = i;
        auto es = pow(1.0 + i, -params.zipf_exponent) * zipf_coef;
        distr.inc_weight(i, es*exp(stdnormal()/8.0));
        cities[i].expected_size = es;
      };
      c_v = (exp(params.variance_of_social_activity) - 1.) * pow(zipf_coef, -theta_v);
      for (unsigned i = 0; i < params.population_size; i++) {
        ++(cities[distr.rvs()].size);
      };
    }

    person_city = new unsigned int[params.population_size];
    for (unsigned i = 0, from = 0; i < params.number_of_cities; ++i) {
      cities[i].begin_from = from;
      std::fill_n(person_city + from, cities[i].size, i);
      from += cities[i].size;
    }

 }

  ~CitySet()
  {
    N = 0;
    MY_DEL_ARR(cities);
    MY_DEL_ARR(person_city);
  }

  City& at(unsigned int c) const
  {
    return cities[c];
  };

  //void inc_active(unsigned p) { ++(cities[personc(p)].I_active); };
  //void dec_active(unsigned p) { --(cities[personc(p)].I_active); };

  void init_intercity_prob(double p_intercity_spr,double ca_cum[]) {
    auto coef = ca_cum[total_population - 1];
    double cur = 0;
    for (unsigned c = 0; c < N; ++c) {
      auto e = ca_cum[end(c) - 1];
      cities[c].intercity_prob = p_intercity_spr*coef/(coef+cur-e);
      cur=e;
    }

  };

  double p_intercity_spread(unsigned c) const {
    return   cities[c].intercity_prob;
  }
  unsigned personc(unsigned int p) const { return person_city[p]; };

  /*double get_R0_coef_log(unsigned c) {
    return  theta_r * log(cities[c].expected_size / cities[0].expected_size);
  }; */

  double get_vs_exp(unsigned c) {
    return (1.+c_v* pow(cities[c].size, theta_v));
  };

  void multiply_by_social_activity(const Params &params,double catchability[],double spreadability[]) {
    double main_city_size = cities[0].size;
    auto main_city_tau_corrected_r0 = params.initial_R0 + params.inverse_tau - 1.;
    for (unsigned c = 0; c < N; ++c) {
      auto exp_v = get_vs_exp(c);
      auto v = log(exp_v);
      auto m = -0.5 * v, sv = sqrt(v);
      auto city_tau_corrected_r0 = main_city_tau_corrected_r0 * pow(cities[c].size / main_city_size, theta_r);
      auto spr_coef = city_tau_corrected_r0/get_as0(v);
      auto cat_coef = pow(cities[c].size / main_city_size, theta_i);
      for (unsigned person = begin(c); person < end(c); ++person) {
        auto sa = exp(stdnormal() * sv + m);
        spreadability[person] *= sa * spr_coef;
        catchability[person] *= sa * cat_coef;
      }
    }

  };

  QuartileBounds calc_quartiles_bounds(unsigned n[]) {
    unsigned b = 0;
    unsigned sum = cities[b].size;
    while (sum < total_population / 4)
      sum += cities[++b].size;
    n[0] = sum - cities[b].size;
    unsigned m = b;
    while (sum < total_population / 2)
      sum += cities[++m].size;
    n[1] = sum - n[0] - cities[m].size;
    unsigned s = m;
    while (sum < 3*total_population / 4)
      sum += cities[++s].size;
    n[2] = sum - n[0] - n[1] - cities[s].size;
    n[3] = total_population - sum + cities[s].size;
    return QuartileBounds(b, m, s);
  };

  unsigned begin(unsigned c) const { return cities[c].begin_from; }
  unsigned end(unsigned c) const { return cities[c].begin_from + cities[c].size; }

 
};


