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

int make_trj(std::vector<fvector<4>> _xd, std::vector<xd_t> &_v);
int make_trj(std::vector<fvector<3>> _xd, std::vector<xd_t> &_v);