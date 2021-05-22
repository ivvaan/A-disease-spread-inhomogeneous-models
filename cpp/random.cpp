#include <random>
#include "random.h"

#ifdef NDEBUG
std::random_device rd_gen;
std::default_random_engine drd_gen;
#else
std::default_random_engine rd_gen;
std::default_random_engine drd_gen;
#endif // NDEBUG




double stduniform1()
{
    static std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
    return dist(rd_gen);
}
double stdnormal1()
{
    static std::normal_distribution<double> dist{ 0.0, 1.0 };
    return dist(rd_gen);
}

double stduniform2()
{
    static std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
    return dist(drd_gen);
}
double stdnormal2()
{
    static std::normal_distribution<double> dist{ 0.0, 1.0 };
    return dist(drd_gen);
}

double (*stduniform)() = stduniform1;

double (*stdnormal)() = stdnormal1;

//zero seed means use of random gen
//nonzero - use of pseudo random setting seed
void set_seed(unsigned seed)
{
    if (seed) {
        drd_gen.seed(seed);
        stduniform = stduniform2;
        stdnormal = stdnormal2;
        return;
    }
    stduniform = stduniform1;
    stdnormal = stdnormal1;
};