#include <cassert>
#include "compartment.h"
#include "random.h"
#include "my_macros.h"

Compartment::Compartment(unsigned int ps)
    : size{ 0 }
    , n_removed{ 0 }
 {
    members = new unsigned int[ps];
 };


Compartment::~Compartment()
{
    MY_DEL_ARR(members);
};

void Compartment::add(unsigned int person)
{
    members[size++] = person;
};

void Compartment::clear_removed()
{
    auto isn = members;
    auto isl = isn + size;
    for (auto iso = isn; iso < isl; ++iso)
        if (auto el = *iso; el != removed)
            *(isn++) = el;
    size = isn - members;
    n_removed = 0;
}

unsigned Compartment::remove()
{
    assert(n_removed < size);
    double dsize = size;
    unsigned int* removing;
    do 
        removing = &members[unsigned int(stduniform() * dsize)];
    while (*removing == removed);
    auto res = *removing;
    *removing = removed;
    if (2 * ++n_removed > size)  clear_removed();
    return res;
}


