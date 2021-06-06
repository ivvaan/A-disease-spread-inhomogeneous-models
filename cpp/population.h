#pragma once
#include <memory>
#include <climits>
class CPopulation {
public:
    unsigned int pop_size;
    double ca_sum, vs;
    double *catchability, *spreadability;
    CPopulation(unsigned int ps, double v_sa, double v_bc, double v_bs);
    ~CPopulation();
    auto get_size() const { return pop_size; };
    auto get_as0() const { return 1.0 + vs * (1.0 + 0.5 * vs); };
    //auto get_as0() const { return exp(vs); };

    unsigned int candidate_to_infect() const;

    double get_spreadability(int person) const { return spreadability[person]; };


};
