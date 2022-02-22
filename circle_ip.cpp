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

    operator fvector<3> &() { return xd; }
    void print(void) const{
        xd.print();
        if(this->type == Arc) arc.print();
    }
};

fmatrix<2, 2> Rl90(0, -1, 1, 0);

int make_trj(std::vector<fvector<4>> _xd, std::vector<xd_t> &_v){
    _v.clear();
    int num = _xd.size();
    for(int k = 1;k < num;k++){
        fvector<2> s(_xd[k - 1][0], _xd[k - 1][1]), e(_xd[k][0], _xd[k][1]);
        float s_ang = _xd[k - 1][3], e_ang = _xd[k][3];

        xd_t val;

        if(s_ang == e_ang){
            val.xd = {e[0], e[1], _xd[k][2]};
            val.type = xd_t::Line;
            _v.push_back(val);
        }
        else {

            fmatrix<2, 2> A(cosf(s_ang), cosf(e_ang), sinf(s_ang), sinf(e_ang));
            fvector<2> b = A.solve(e - s);

            fvector<2> x_cross = s + b[0] * fvector<2>(cos(s_ang), sin(s_ang));

            if(fabsf((x_cross - s).norm()-(x_cross - e).norm())< 0.001){
                A = fmatrix<2, 2>(
                    - sinf(e_ang) + sinf(s_ang), cosf(e_ang),
                    cosf(e_ang) - cosf(s_ang), sinf(e_ang)
                );
                b = A.solve(s - x_cross);

                fvector<2> x_o = s + b[0] * fvector<2>(- sinf(s_ang), cosf(s_ang));
                x_cross += b[1] * fvector<2>{cosf(e_ang), sinf(e_ang)};
                
                val.xd = {x_cross[0], x_cross[1], _xd[k][2]};
                val.type = xd_t::Arc;
                val.arc = {x_o[0], x_o[1], fabsf(b[0]), s_ang - pi/2, e_ang - pi/2};
                _v.push_back(val);
            }
            else if((x_cross - s).norm()<(x_cross - e).norm()){
                A = fmatrix<2, 2>(
                    - sinf(e_ang) + sinf(s_ang), cosf(e_ang),
                    cosf(e_ang) - cosf(s_ang), sinf(e_ang)
                );
                b = A.solve(s - x_cross);

                fvector<2> x_o = s + b[0] * fvector<2>(- sinf(s_ang), cosf(s_ang));
                x_cross += b[1] * fvector<2>{cosf(e_ang), sinf(e_ang)};
                
                val.xd = {x_cross[0], x_cross[1], _xd[k][2]};
                val.type = xd_t::Arc;
                val.arc = {x_o[0], x_o[1], fabsf(b[0]), s_ang - pi/2, e_ang - pi/2};
                _v.push_back(val);

                val.xd = {e[0], e[1], _xd[k][2]};
                val.type = xd_t::Line;
                _v.push_back(val);
            }
            else {
                A = fmatrix<2, 2>(
                    sinf(e_ang) - sinf(s_ang), cosf(s_ang),
                    - cosf(e_ang) + cosf(s_ang), sinf(s_ang)
                );
                b = A.solve(e - s);

                fvector<2> x_o = e + b[0] * fvector<2>(- sinf(e_ang), cosf(e_ang));
                x_cross = s + b[1] * fvector<2>{cosf(s_ang), sinf(s_ang)};
                
                val.xd = {x_cross[0], x_cross[1], _xd[k][2]};
                val.type = xd_t::Line;
                _v.push_back(val);

                val.xd = {e[0], e[1], _xd[k][2]};
                val.type = xd_t::Arc;
                val.arc = {x_o[0], x_o[1], fabsf(b[0]), s_ang - pi/2, e_ang - pi/2};
                _v.push_back(val);
            }
        }
    }
}

int main(void){

    std::vector<fvector<4>> xd = {
        {0, 0, 0, 0},
        {2, 1.5, pi/4, pi/2},
        {2, 2.5, pi/4, pi/2},
        {1, 4, 0, 3*pi/4}
    };

    std::vector<xd_t> ret;
    printf("%d\r\n", make_trj(xd, ret));

    for(auto v : ret)v.print();

    FILE *fp = fopen("c.dat", "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", "c.dat");
    fprintf(fp, "# control points\r\n");

    for(auto v : ret){
        if(v.type == xd_t::Line){
            fprintf(fp, "0 %6.3f %6.3f 0 0 0 ", v.xd[0], v.xd[1]);
        }
        else {
            fprintf(fp, "1 %6.3f %6.3f %6.3f %6.3f %6.3f ", v.arc[0], v.arc[1], v.arc[2], v.arc[3], v.arc[4]);
        }
        fprintf(fp, "\r\n");
    }

    fclose(fp);
    printf("file closed.\r\n");

    return 0;
}