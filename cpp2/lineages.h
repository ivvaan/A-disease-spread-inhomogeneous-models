#pragma once
#include "random.h"
#include "my_macros.h"
#include "params.h"
#include "mutations.h"
#include "person_status.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>




struct Lineage {
public:
    unsigned int ID;
    unsigned int I_active;
    unsigned int ancestor, generation;
    double T_birth;
    double contagiosity;

    Lineage(){};
 
    Lineage(const Lineage& l)
        : ID(l.ID)
        , I_active(l.I_active)
        , ancestor(l.ancestor)
        , generation(l.generation)
        , T_birth(l.T_birth){};

    ~Lineage(){};

    void print_details(std::ostream& os) {
      os << " " << ancestor << " " << generation << " " << T_birth << " " << contagiosity;
    };

    double add_first(double cont, double dA, double T, unsigned int id, unsigned int a, unsigned int g)
    {
      ID = id;
      //I_total = 1;
      I_active = 1;
      T_birth = T;
      ancestor = a;
      generation = g;
      contagiosity = cont;
      return dA*cont;
    };

    unsigned int get_generation() const { return generation; };

    double add(double dA)
    {
        ++I_active;
        //++I_total;
        return contagiosity*dA;
    };

    void recover()
    {
        //--I_total;
        //assert((I_total >= I_active));
    }

    double deactivate(double dA)
    {
        --I_active;
        return dA;
    };
   /* std::ostream& print_n(int n, std::ostream& os) const
    {
        return os << n << " " << ID << " " << I_total << " " << T_birth << " " << ancestor << " " << generation << "\n";
    }; */
};


class Lineages {
protected:
  
  int N;
  double T_prev;


  Lineage* lineages;
  unsigned int* person_lineage;

  double immunity_loss_rate = 0.04;
  double immunity_loss_power = 2.0;
  double protective_immunity_interval = 4.;
public: 
  static constexpr unsigned external_lineage = 0;

  Lineages(const Params& params)
    : N(1),
    T_prev(0),
    immunity_loss_rate(params.immunity_loss_rate),
    immunity_loss_power(params.immunity_loss_power),
    protective_immunity_interval(params.protective_immunity_interval)

  {
    auto ps = params.population_size;
    lineages = new Lineage[ps];

    lineages->ID = 0;
    lineages->I_active = 0;
    lineages->ancestor = 0;
    lineages->generation = 0;
    lineages->T_birth = 0;
    lineages->contagiosity = 1.;

    person_lineage = new unsigned int[ps];
    std::fill_n(person_lineage, ps, 0);
  }

  ~Lineages()
  {
    N = 0;
    MY_DEL_ARR(lineages);
    MY_DEL_ARR(person_lineage);
  }

  Lineage& at(unsigned int l) const
  {
    return lineages[l];
  };

  int personl(unsigned int p) const { return person_lineage[p]; };

  unsigned int lineage_numb() const { return N; }

  int try_to_infect(unsigned candidate_to_infect,unsigned spreader_lineage, const SIRStatus::ST& p_info, double T) const {
    if (p_info.status == SIRStatus::recovered) {
      auto dt = T - p_info.recovery_date - protective_immunity_interval;
      if (dt < 0)return SIRStatus::cant_be_infected;
      
      //it could be something dependant on spreder_lineage and recovered_person_lineage
      //auto recovered_person_lineage = personl(candidate_to_infect);
      
      if (
        stduniform() < exp(-pow(dt * immunity_loss_rate, immunity_loss_power))
        //stduniform() * (1.+pow(dt * immunity_loss_rate, immunity_loss_power)) < 1.0
        ) return SIRStatus::cant_be_infected;
    };
    return p_info.status;
  };


};

template<class GS>
class LineageSet:public Lineages {
  GS genomes;
public:

  LineageSet(const Params &params)
    :Lineages(params),
    genomes(params)
  {  }

  ~LineageSet()
  {
  }
 double add_new(double dA, double T, unsigned int cs, unsigned int a) {
    person_lineage[cs] = N;
    dA=lineages[N].add_first(genomes.originate_from(N,a),dA, T, N, a, lineages[a].get_generation() + 1);
    ++N;
    return dA;
  };

  double add(double dA, double new_prop, double T, unsigned int c_infected, unsigned int lineage)
  {
    new_prop *= T_prev - T; //negative!!!!!
    T_prev = T;

    if (exp(new_prop) + stduniform() < 1.)
      return add_new(dA, T, c_infected, lineage);

    person_lineage[c_infected] = lineage;
    return lineages[lineage].add(dA);
  };

  void recover(unsigned int cs)
  {
    lineages[personl(cs)].recover();
  };

  double deactivate(double dA, unsigned int c_infected)
  {
    auto lineage = personl(c_infected);
    return lineages[lineage].deactivate(dA);
  };
};

