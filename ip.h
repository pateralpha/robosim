#pragma once

#include <stdio.h>
#include <vector>
#include <array>
#include <cmath>
#include "la.h"


struct xd_t{
    fvector<3> xd;
    enum type_t{
        Line,
        Arc,
    } type;

    fvector<5> arc;
    float v_m, v_e;

    operator fvector<3> &() { return xd; }
    void print(void) const{
        xd.print();
        if(this->type == Arc) arc.print();
    }
};

int make_trj(std::vector<fvector<4>> _xd, std::vector<xd_t> &_v, float _vm, float _ve = 100);
int make_trj(std::vector<fvector<3>> _xd, std::vector<xd_t> &_v, float _vm, float _ve = 0);

fvector<3> ragrange(std::vector<xd_t> _xd_vec, float _u);
float drdu(std::vector<xd_t> _xd_vec, float _u);
float arg_rag(std::vector<xd_t> _xd_vec, float _u);
float v_max(float vm, float ddx, std::vector<xd_t> _xd_vec, float _u);