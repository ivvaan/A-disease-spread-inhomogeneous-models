#pragma once

#include "my_macros.h"
#include <memory>

class CImmuneStatus {
public:
  char* set;
  CImmuneStatus(unsigned int ps) {
    set = new char[ps];
    memset(set, 0, ps * sizeof(set[0]));
  };
  ~CImmuneStatus()
  {
      MY_DEL_ARR(set);
  };
  bool is_immune(unsigned int person) const { return set[person]; };
  void set_immune(unsigned int person) { set[person] = 1; };
  void set_susceptible(unsigned int person) { set[person] = 0; };
};
