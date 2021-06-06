#pragma once
#define MY_DEL_ARR(a) \
    if (a)            \
        delete[] a;   \
    a = nullptr;
