#include "trace.h"
#include "lineages.h"


#define PRINT_lineage_INFO

std::ostream& TL::print_test(int n, Lineages& lineages, std::ostream& os) const
{
  //if (lineages.at(L).contagiosity <= 1.) return os;
  //if (L == 1)return os;
    os << n << " " << T << " " << L;
#ifdef PRINT_lineage_INFO
    lineages.at(L).print_details(os);
#endif
    os << "\n";
    return os;
};

std::ostream& operator<<(std::ostream& os, const TSI& tsi) {
  os << tsi.T << " " << tsi.S << " " << tsi.I;
  return os;
};

std::ostream& operator<<(std::ostream& os, const TSIQ& tsiq) {
  os << tsiq.T << " " <<
    tsiq.S[0] << " " << tsiq.S[1] << " " << tsiq.S[2] << " " << tsiq.S[3] <<
    " "<< tsiq.I[0] << " " << tsiq.I[1] << " " << tsiq.I[2] << " " << tsiq.I[3];
    return os;
};

std::ostream& operator<<(std::ostream& os, const TNRI& tsi) {
  os << tsi.T << " " <<tsi.NI<<" "<<tsi.RI<<" "<<tsi.I;
  return os;
};


/*

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


*/