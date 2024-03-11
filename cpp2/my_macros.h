#pragma once
#define SQURE(x) ((x)*(x))
#define MY_DEL_ARR(a) \
    if (a)            \
        delete[] a;   \
    a = nullptr;
