#pragma once
#include <algorithm>
extern double (*stduniform)();
extern double (*stdnormal)();
void set_seed(unsigned seed);

//Heap like 2^P-ary tree (binary for P=1, quaternary for P=2 ets.).
//For millions of categories the structure does not fit to the processor cache,
//so each the tree level needs expensive reading from main memory.
//To reduce the tree height in P times we use 2^P-ary trees.
template <unsigned P> 
class DynamicCategiricalDistr {

  constexpr static unsigned S = (1 << P) - 1;//number of values to compare in a node
  int Nt = 0, Nv = 0, size = 0;
  double(*tree)[S] = nullptr;
  double* values = nullptr;

  static void add2tree(unsigned p, const double v, double(*tree)[S]) {
    while (p) {
      auto q = p & S;
      p >>= P;
      if (q != S)tree[p][q] += v;
    }
  };

  void Init(unsigned s) {
    Nv = s;
    for (Nt = 1; s; s >>= P) Nt <<= P;
    if ((Nv << P) == Nt)
      Nt = Nv;
    size = ((Nt + Nv - 1) >> P) + 1;
    values = new double[Nv];
    std::fill_n(values, Nv, 0);
    tree = new double[size][S];
    std::fill_n(reinterpret_cast<double*>(tree), size * S, 0);
  }

public:

  DynamicCategiricalDistr(int s) {
    if (s <= 0)return;
    Init(s);
  };

  DynamicCategiricalDistr(int s, double p[]) {
    if (s <= 0)return;
    Init(s);
    for (int pos = 0; pos < s; ++pos)
      inc_weight(pos, p[pos]);
  };

  ~DynamicCategiricalDistr() {
    if (tree) { delete[] tree; tree = nullptr; }
    if (values) { delete[] values; values = nullptr; }
    Nv = Nt = size = 0;
  };

  double inc_weight(int pos, double v) { //adding v>0
    add2tree(pos + Nt, v, tree);
    values[pos] += v;
    return v;
  };

  double dec_weight(int pos, double v) { //subtracting v>0
    if (v > values[pos])
      v = values[pos];
    add2tree(pos + Nt, -v, tree);
    values[pos] -= v;
  };

  void set_weight(int pos, double v) {
    if (v < 0)
      v = 0;
    add2tree(pos + Nt, v - values[pos], tree);
    values[pos] = v;
  };

 unsigned int get_cat(double v) {
    auto T = tree;
    unsigned pos = 1;
    do {
      auto& t = T[pos];
      pos <<= P;
      for (auto c : t) {
        if (c > v)break;
        v -= c; ++pos;
      }
    } while (pos < Nt);
    return pos - Nt;
  };

  double weights_sum() {
    return tree[0][1];
  };

  int rvs() {
    return get_cat(stduniform()* weights_sum());
  };

};

template<>
class DynamicCategiricalDistr<1> {
  //public:
  unsigned int Nt = 0, Nv = 0, size = 0;
  double* tree = nullptr;
  double* values = nullptr;

  static void add2tree(unsigned int p, const double v, double* tree) {
/*    for (; p; tree[p >>= 1] += v)
      while (p & 1)
        p >>= 1; */
    do {
      auto w = ((p & 1) ^ 1) * v;
      tree[p >>= 1] += w;
    } while (p);
    tree[0] += v;
  };

  void Init(unsigned int s) {
    Nv = s;
    for (Nt = 1; s; s >>= 1) Nt <<= 1;
    if (2 * Nv == Nt)
      Nt = Nv;
    size = (Nt + Nv - 1) / 2 + 1;
    values = new double[Nv];
    std::fill_n(values, Nv, 0);
    tree = new double[size];
    std::fill_n(tree, size, 0);
  };

public:

  DynamicCategiricalDistr(unsigned int s) {
    if (s <= 0)return;
    Init(s);
  };

  DynamicCategiricalDistr(unsigned int s, double p[]) {
    if (s <= 0)return;
    Init(s);
    for (unsigned int pos = 0; pos < s; ++pos)
      inc_weight(pos, p[pos]);
  };

  ~DynamicCategiricalDistr() {
    if (tree) { delete[] tree; tree = nullptr; }
    if (values) { delete[] values; values = nullptr; }
    Nv = Nt = size = 0;
  };

  double inc_weight(unsigned int pos, double v) { //adding v>0
    add2tree(pos + Nt, v, tree);
    values[pos] += v;
    return v;
  };

  double dec_weight(unsigned int pos, double v) { //subtracting v>0
    if (v > values[pos])
      v = values[pos];
    add2tree(pos + Nt, -v, tree);
    values[pos] += v;
    return v;
  };

  void set_weight(unsigned int pos, double v) {
    if (v < 0)
      v = 0;
    add2tree(pos + Nt, v - values[pos], tree);
    values[pos] = v;
  };

  unsigned int rvs() {
    return get_cat(stduniform() * weights_sum());
  };

  double weights_sum() {
    return tree[0];
  };

  unsigned int get_cat(double v) {
    auto T = tree, _T = tree - 1, C = tree + 1, E = tree + Nt;
    //E - end position aka &tree[Nt]
    //C - current position aka &tree[pos]
    do {
      if (*C < v) {
        v -= *C;
        C += C - _T; //pos=pos*2+1;    right branch
      }
      else
        C += C - T;  //pos=pos*2;      left branch
    } while (C < E);
    return C - E;     //pos - Nt
  };
};

template<>
unsigned int DynamicCategiricalDistr<2>::get_cat(double v);

template<>
unsigned int DynamicCategiricalDistr<3>::get_cat(double v);


