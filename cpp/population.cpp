#include "population.h"
#include "random.h"
#include <math.h>


#define MY_DEL_ARR(a) \
    if (a)            \
        delete[] a;   \
    a = nullptr;

double* get_persons_distr(int ps, double v)
{
    auto res = new double[ps];
    auto last = res + ps;
    auto m = -0.5 * v;
    double sv = sqrt(v);
    for (auto p = res; p < last; ++p)
        *p = exp(stdnormal() * sv + m);
    return res;
};

    CPopulation::CPopulation(unsigned int ps, double v_sa, double v_bc, double v_bs)
        : pop_size(ps)
        , vs(v_sa)
    {
        is_infected = new char[ps];
        //memset(is_infected, 0, ps * sizeof(is_infected[0]));

        catchability = get_persons_distr(ps, v_bc);
        spreadability = get_persons_distr(ps, v_bs);
        auto get_sa = [m = -0.5 * v_sa, sv = sqrt(v_sa)]() { return exp(stdnormal() * sv + m); };
        auto psa = get_sa();
        spreadability[0] *= psa;
        psa = (catchability[0] *= psa);
        for (unsigned int i = 1; i < ps; ++i) {
            auto csa = get_sa();
            spreadability[i] *= csa;
            catchability[i] *= csa;
            psa = (catchability[i] += psa);
        }
        ca_sum = psa;
    }

    CPopulation::~CPopulation()
    {
        MY_DEL_ARR(is_infected);
        MY_DEL_ARR(catchability);
        MY_DEL_ARR(spreadability);
    };


    int CPopulation::spreader_candidate() const
    {
        int m, l = 0, r = pop_size - 1;
        double s = ca_sum * stduniform();
        while (r > l) {
            m = (l + r) / 2;
            if (s < catchability[m])
                r = m;
            else
                l = m + 1;
        }
        return l;
    };

    int CPopulation::first_speader()
    {
        auto sc = spreader_candidate();
        is_infected[sc] = 1;
        return sc;
    };

    int CPopulation::next_speader()
    {
        auto sc = spreader_candidate();
        if (is_infected[sc])
            return -1;
        is_infected[sc] = 1;
        return sc;
    };
