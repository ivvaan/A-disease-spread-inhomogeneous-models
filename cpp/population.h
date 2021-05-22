#pragma once
#include <memory>
class CPopulation {
public:
    unsigned int pop_size;
    char* is_infected;
    double ca_sum, vs;
    double *catchability, *spreadability;
    void refresh() { memset(is_infected, 0, pop_size * sizeof(is_infected[0])); };
    double count_infected() const
    {
        double res = 0;
        for (char* cur = is_infected; cur < is_infected + pop_size; ++cur)
            res += *cur;
        return res;
    }
    CPopulation(unsigned int ps, double v_sa, double v_bc, double v_bs);
    ~CPopulation();
    auto get_size() const { return pop_size; };
    auto get_as0() const { return 1.0 + vs * (1.0 + 0.5 * vs); };
    //auto get_as0() const { return exp(vs); };

    double get_spreadability(int person) const { return spreadability[person]; };

    int spreader_candidate() const;
 
    int first_speader();

    int next_speader();

    void make_susceptible(unsigned int person) { is_infected[person] = 0; };
};
