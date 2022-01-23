#pragma once
#include "cities.h"
#include "params.h"
#include <memory>
#include <climits>

class Population {
public:
    unsigned int pop_size=0;
    double ca_sum=1, vs=1;
    double *catchability=nullptr, *spreadability=nullptr;
    Population(const Params &params);
    Population() {};
    ~Population();
    auto get_size() const { return pop_size; };
    unsigned person_city(unsigned person)const { return 0; };

    unsigned int candidate_to_infect() const;

    unsigned int candidate_to_infect(unsigned spreader) const {
      return candidate_to_infect();
    };

    unsigned int candidate_to_infect(unsigned f, unsigned t) const;

    double get_spreadability(int person) const { return spreadability[person]; };
    
    CitySet::QuartileBounds calc_quartiles_bounds(unsigned n[]) {
      n[0] = pop_size;
      return { pop_size ,pop_size + 1,pop_size + 2 };
    };

};

class CountryPopulation:public Population {
  CitySet cities;
public:

  using Population::candidate_to_infect;
  using Population::get_size;

  CountryPopulation(const Params& params);
  ~CountryPopulation() {};
 /* auto get_size() const { return pop_size; };
  */
  unsigned int candidate_to_infect(unsigned spreader) const;
  unsigned person_city(unsigned person)const { return cities.personc(person); };
  CitySet::QuartileBounds calc_quartiles_bounds(unsigned n[]) {
    return cities.calc_quartiles_bounds(n); 
  }

 
};
