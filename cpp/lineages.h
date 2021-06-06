#pragma once
#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>

#include "random.h"
#include "my_macros.h"

struct Clineage {
public:
    unsigned int ID;
    //unsigned int I_total,
    unsigned int I_active;
    double A;
    unsigned int ancestor, generation;
    double T_birth;

    Clineage(){};
 
    Clineage(const Clineage& l)
        : ID(l.ID)
        //, I_total(l.I_total)
        , I_active(l.I_active)
        , A(l.A)
        , ancestor(l.ancestor)
        , generation(l.generation)
        , T_birth(l.T_birth){};

    ~Clineage(){};

    double add_first(double dA, double T, unsigned int id, unsigned int a, unsigned int g)
    {
        A = dA;
        ID = id;
        //I_total = 1;
        I_active = 1;
        T_birth = T;
        ancestor = a;
        generation = g;
        return dA;
    };

    unsigned int get_generation() const { return generation; };

    double add(double dA)
    {
        A += dA;
        ++I_active;
        //++I_total;
        return dA;
    };

    void recover()
    {
        //--I_total;
        //assert((I_total >= I_active));
    }

    double deactivate(double dA)
    {
        A -= dA;
        --I_active;
        return dA;
    };
   /* std::ostream& print_n(int n, std::ostream& os) const
    {
        return os << n << " " << ID << " " << I_total << " " << T_birth << " " << ancestor << " " << generation << "\n";
    }; */
};

class ClineageSet {

    class CReorder {
    public:

        unsigned int nzN;
        unsigned int* order;
        CReorder(int ps)
            : nzN{ 0 }
        {
            order = new unsigned int[ps];
        }
        ~CReorder()
        {
            nzN = 0;
            MY_DEL_ARR(order);
        }
        void add(int N) { order[nzN++] = N; }

        unsigned int select(double A, Clineage* lineages);
    };

    CReorder spawn_order;

public:
    int N;
    double T_prev;
    Clineage* lineages;
    unsigned int* person_lineage;
    ClineageSet(unsigned int ps)
        : N{ 1 }
        , T_prev{ 0 }
        , spawn_order(ps)
    //, recover_order(ps)
    //, test_order(ps)

    {
        lineages = new Clineage[ps];
        
        lineages->ID=0;
        //lineages->I_total=0;
        lineages->I_active=0;
        lineages->A=0;
        lineages->ancestor=0;
        lineages->generation=0;
        lineages->T_birth=0;

        person_lineage = new unsigned int[ps];
    }

    ~ClineageSet()
    {
        N = 0;
        MY_DEL_ARR(lineages);
        MY_DEL_ARR(person_lineage);
    }

    Clineage& at(unsigned int l) const
    { 
      return lineages[l]; 
    };

    int personl(unsigned int p) const { return person_lineage[p]; };

    unsigned int lineage_numb() const { return N; }

    //int get_active(int l) const { return lineages[l].I_active; }
    //int get_total(int l) const { return lineages[l].I_total; }

    double add_new(double dA, double T, unsigned int cs, unsigned int a)
    {
        spawn_order.add(N);
        person_lineage[cs] = N;
        lineages[N].add_first(dA, T, N, a, lineages[a].get_generation() + 1);
        ++N;
        return dA;
    };

    double add(double dA, double new_prop, double A, double T, unsigned int cs)
    {
        new_prop *= T_prev - T; //negative!!!!!
        T_prev = T;
        auto lineage = spawn_order.select(A, lineages);
        if (exp(new_prop) + stduniform() < 1.)
            return add_new(dA, T, cs, lineage);
        person_lineage[cs] = lineage;
        return lineages[lineage].add(dA);
    };

    void recover(unsigned int cs)
    {
        lineages[person_lineage[cs]].recover();
    };

    double deactivate(double dA, unsigned int cs)
    {
        return lineages[person_lineage[cs]].deactivate(dA);
    };
};
