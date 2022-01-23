#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <type_traits>
#include <utility>
#include "cities.h"


struct TSI {
    double T;
    unsigned S, I;
    TSI(): T(0), S(0), I(0){};
    TSI(double t, unsigned s, unsigned i)
        : T(t), S(s), I(i){};
    std::ostream& print_nt(int n, double t, std::ostream& os) const { 
    return os << n << " " << t << " " << S << " " << I<<"\n"; 
    };
    friend std::ostream& operator<<(std::ostream& os, const TSI& tsi);
};

struct TSIQ {
  double T;
  std::array<unsigned, 4> S;
  std::array<unsigned, 4> I;
  TSIQ()
    : T(0)
    , S{ 0, 0, 0, 0 }
  , I{ 0, 0, 0, 0 } {};
  TSIQ(double t,const std::array<unsigned, 4>&  s, const std::array<unsigned, 4>& i)
    : T(t)
    , S(s)
  , I{ i } {};
  friend std::ostream& operator<<(std::ostream& os, const TSIQ& tsiq);

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

class Lineages;

struct TL {
    double T;
    int L;
    TL() : T(0) , L(0){};
    TL(double t, int l)
        : T(t) , L(l){};
    std::ostream& print_test(int n, Lineages& lineages, std::ostream& os) const;
};

//std::ostream& print_resampled(std::vector<TSI>& trace, double step, int n, std::ostream& os);

struct TNRI {
  double T;
  float NI, RI, I; //NI-naive infected, RI-recovered infected
  TNRI(): T(0), NI(0), RI(0), I(0) {};
  TNRI(double t, float n, float r, float i)
    : T(t), NI(n), RI(r),I(i) {};

  friend std::ostream& operator<<(std::ostream&, const TNRI&);
};

template <class Trace>
auto get_resampler(Trace& trace,double step_)  {
  return [t = (double)0.0,step=step_, &trace](double T, unsigned S, unsigned I) mutable {
    while(T>t){
      trace.emplace_back(t, S, I);
      t += step;
    }
  };
};


template <class Population>
std::function<void(double, unsigned, unsigned)> get_resampler_ex(const Params &params,
  Population& population,
  std::vector<TSI>& main_trace, 
  std::vector<TSIQ>& qtrace, 
  std::vector<TNRI>& trace_nri) {
  switch (params.additional_output)
  {
  case 1:
  {
    unsigned N[4];
    CitySet::QuartileBounds qb = population.calc_quartiles_bounds(N);
    return[s = N[0] + N[1] + N[2] + N[3], i = 0, S = std::array<unsigned, 4>{N[0], N[1], N[2], N[3]}, I = std::array<unsigned, 4>{}, t = 0.0, step = params.resample_step, & main_trace, & qtrace, qb](double T, unsigned c, unsigned op) mutable {
      auto quartile = qb.classify(c);
      unsigned dS = op & 1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;;
      int dI = (int)((op + 2) & 2) - 1; //op==incI->(op+2)&2-1==1;op==decSincI->(op+2)&2-1==1;op==decI->(op+2)&2-1==-1;
      S[quartile] -= dS;
      I[quartile] += dI;
      s -= dS;
      i += dI;
      while (T > t) {
        main_trace.emplace_back(t, s, (unsigned)i);
        qtrace.emplace_back(t, S, I);
        t += step;
      }
    };
  };
  case 2:
    return[N = params.population_size, R = 0, NI = 0u, RI = 0u, 
      t = 0.0, step = params.resample_step, &main_trace, &trace_nri](double T, unsigned c, unsigned op) mutable {
      unsigned dS = op & 1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;
      //S -= dS;
      R += (int)op - 1; //op==incI->op - 1==-1;op==decSincI->op - 1==0;op==decI->op - 1==1;
      NI += dS;
      RI += (dS | (op >> 1)) ^ 1;
      while (T > t) {
        unsigned S = N - NI;
        unsigned I = NI - R;
        main_trace.emplace_back(t, S, I);
        
        trace_nri.emplace_back(t, NI/ (double)N, RI/ (double)N,I/(double)S);
        t += step;
      }
    };
  default:
    return[S = params.population_size, I = 0, t = 0.0, step = params.resample_step, &main_trace](double T, unsigned c, unsigned op) mutable {
      S -= op & 1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;
      I += (int)((op + 2) & 2) - 1; //op==incI->(op+2)&2-1==1;op==decSincI->(op+2)&2-1==1;op==decI->(op+2)&2-1==-1;
      while (T > t) {
        main_trace.emplace_back(t, S, (unsigned)I);
        t += step;
      }
    };    break;
  }

};




/*
auto get_resampler2(std::vector<TSI>& trace, double step_,unsigned N) {
  return[S=N,I=0,t = 0.0, step = step_, &trace](double T,unsigned c, unsigned op) mutable {
    S -= op&1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;
    I += (int)((op + 2) & 2) - 1; //op==incI->(op+2)&2-1==1;op==decSincI->(op+2)&2-1==1;op==decI->(op+2)&2-1==-1;
    while (T > t) {
      trace.emplace_back(t, S, (unsigned)I);
      t += step;
    }
  };
};


template <class Population>
auto get_city_quartiles_resampler(std::vector<TSI>& main_trace, std::vector<TSIQ>& qtrace, double step_, Population &population) 
{
  unsigned N[4];
  CitySet::QuartileBounds qb= population.calc_quartiles_bounds(N);
  return[s= N[0]+N[1]+N[2]+N[3],i=0,S = std::array<unsigned, 4>{N[0], N[1], N[2], N[3]}, I = std::array<unsigned, 4>{}, t = 0.0, step = step_,&main_trace, & qtrace, qb](double T, unsigned c, unsigned op) mutable {
    auto quartile = qb.classify(c);
    unsigned dS = op & 1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;;
    int dI = (int)((op + 2) & 2) - 1; //op==incI->(op+2)&2-1==1;op==decSincI->(op+2)&2-1==1;op==decI->(op+2)&2-1==-1;
    S[quartile] -= dS;
    I[quartile] += dI;
    s -= dS;
    i += dI;
    while (T > t) {
      main_trace.emplace_back(t, s, (unsigned)i);
      qtrace.emplace_back(t, S, I);
      t += step;
    }
  };
};




auto get_ie_resampler(std::vector<TNRI>& trace, double step_, unsigned N) {
  return[S = N, I = 0,NI=0u,RI=0u, t = 0.0, step = step_, &trace](double T, unsigned c, unsigned op) mutable {
    unsigned dS= op & 1; //op==incI->op&1==0;op==decSincI->op&1==1;op==decI->op&1==0;
    S -= dS;
    I += (int)((op + 2) & 2) - 1; //op==incI->(op+2)&2-1==1;op==decSincI->(op+2)&2-1==1;op==decI->(op+2)&2-1==-1;
    NI += dS;
    RI += (dS|(op>>1))^1;
    while (T > t) {
      trace.emplace_back(t, S, (unsigned)I, NI, RI);
      t += step;
    }
  };
};



template <typename... Funcs>
auto combine_resamplers(Funcs &... funcs) {
  return [&](double T, unsigned c, int op) {
    ((void)funcs(T, c, op), ...);
  };
} */
