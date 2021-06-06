#include "lineages.h"

unsigned int ClineageSet::CReorder::select(double A, Clineage* lineages)
{
    A *= stduniform();
    auto o = order;
    double lA = lineages[o[0]].A;
    if (A < lA)
        return o[0];
    double pA = lA;
    A -= lA;
    for (unsigned int i = 1; i < nzN; ++i) {
        unsigned int lineage = o[i];
        if (lineages[lineage].I_active == 0) {
            do {
                lineage = o[--nzN];
                assert((i < nzN));
            } while (lineages[lineage].I_active == 0);
            o[i] = lineage;
        }
        lA = lineages[lineage].A;
        if (A < lA) {
            if (pA < lA) // lazy stochastic sorting
                std::swap(o[i - 1], o[i]);
            return lineage;
        }
        A -= pA = lA;
    }
    assert(("lineage to select is not found!", false));
    return 0;
};


