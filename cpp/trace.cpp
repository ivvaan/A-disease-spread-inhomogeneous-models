#include "trace.h"
#include "lineages.h"

std::ostream& print_resampled(std::vector<TSI>& trace, double step, int n, std::ostream& os)
{
    auto cur = trace.begin();
    cur->print_nt(n, 0, os);
    auto T = step;
    auto prev = cur;
    for (++cur; cur != trace.end(); prev = cur++)
        for (int counter = 0; cur->T > T; T += step)
            //if ((counter++) % 10 == 0)
                prev->print_nt(n, T, os);
    return os;
};

//#define PRINT_lineage_INFO

std::ostream& TL::print_test(int n, ClineageSet& lineages, std::ostream& os) const
{
    os << n << " " << T << " " << L;
#ifdef PRINT_lineage_INFO
    Clineage& l = lineages.at(L);
    os << " " << l.ancestor << " " << l.generation << " " << l.T_birth;
#endif
    os << "\n";
    return os;
};



