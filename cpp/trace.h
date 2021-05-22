#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>


struct TSI {
    double T;
    int S, I;
    TSI()
        : T(0)
        , S(0)
        , I(0){};
    TSI(double t, int s, int i)
        : T(t)
        , S(s)
        , I(i){};
    std::ostream& print_nt(int n, double t, std::ostream& os) const { 
    return os << n << " " << t << " " << S << " " << I<<"\n"; 
    };
};


/*struct NL {
    int N;
    int L;
    NL()
        : N(0)
        , L(0){};
    NL(int n, int l)
        : N(n)
        , L(l){};
    std::ostream& print_n(int n, std::ostream& os) const
    {
        return os << n << " " << N << " " << L << "\n";
    };
};*/

class ClineageSet;

struct TL {
    double T;
    int L;
    TL()
        : T(0)
        , L(0){};
    TL(double t, int l)
        : T(t)
        , L(l){};
    std::ostream& print_test(int n, ClineageSet& lineages, std::ostream& os) const;
};

std::ostream& print_resampled(std::vector<TSI>& trace,double step, int n, std::ostream& os);

