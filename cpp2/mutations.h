#pragma once
#include "random.h"
#include "params.h"

template <unsigned dim>
double rosenbrock(double x[]) {
  double (*sq)(double) = [](double x) { return x * x; };
  auto res = sq(1. - x[0]) + 100. * sq(x[1] - x[0] * x[0]);
  for(unsigned i=1;i<dim-1;++i)
    res+= sq(1. - x[i]) + 100. * sq(x[i+1] - x[i] * x[i]);
  return res;
}

struct Genome {
  static constexpr unsigned landscape_dim = 3;
  static double get_unnorm_R0(double DNA[]) {
    return landscape_dim / (rosenbrock<landscape_dim>(DNA) + 1.0);
  };
  double DNA[landscape_dim];
  double originate_from(double nonneutral_prob, Genome& g) {
    memcpy(DNA, g.DNA, sizeof(DNA));
    auto rnd = stduniform();
    auto step = rnd < nonneutral_prob ? 0.25 : 1. / 256.;
    for (unsigned i = 0; i < landscape_dim; ++i) {
      DNA[i] += step * stdnormal();
/*    if (rnd < nonneutral_prob) {
      //if (rnd < 0.00003)return 5.;
      for (unsigned i = 0; i < landscape_dim; ++i) {
        DNA[i] += 0.25*stdnormal();
      } */
    }
    return get_unnorm_R0(DNA);
  }

};

class GenomeSet {
  double norm_mult = 1., nonneutral_prob = 0.01;
  Genome* genomes = nullptr;
public:
  GenomeSet(const Params& params):
    nonneutral_prob(params.nonneutral_mutation_prob)
  {
    genomes = new Genome[params.population_size];
    std::fill_n(genomes->DNA, Genome::landscape_dim,0);
    norm_mult = 1. / Genome::get_unnorm_R0(genomes->DNA);
  };
  ~GenomeSet() { MY_DEL_ARR(genomes); };
  double originate_from(unsigned g, unsigned ancestor) {
    return norm_mult*genomes[g].originate_from( nonneutral_prob, genomes[ancestor]);
  }
};

struct TrivGenomeSet {
  TrivGenomeSet(const Params& params) {};
  double originate_from(unsigned g, unsigned ancestor) { return 1; }

};
