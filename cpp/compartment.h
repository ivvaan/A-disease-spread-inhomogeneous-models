#pragma once
#include <climits>
class CCompartment {
    static const unsigned int removed = UINT_MAX;

public:
    unsigned int size, n_removed;
    unsigned int* members;
    auto get_numb() { return size - n_removed; }
    CCompartment(unsigned int ps);
    ~CCompartment();
    void add(unsigned int person);
    void clear_removed();
    unsigned int remove();
};
