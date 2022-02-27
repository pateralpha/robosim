#include "circle.h"



int make_trj(std::vector<fvector<4>> _xd, std::vector<xd_t> &_v, float _vm, float _ve){
    _v.clear();
    int num = _xd.size();
    for(int k = 1;k < num;k++){
        fvector<2> s(_xd[k - 1][0], _xd[k - 1][1]), e(_xd[k][0], _xd[k][1]);
        float s_ang = _xd[k - 1][3], e_ang = _xd[k][3];

        float sign = s_ang > e_ang ? 1 : - 1;

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
                val.arc = {x_o[0], x_o[1], fabsf(b[0]), s_ang + sign * pi/2, e_ang + sign * pi/2};
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
                val.arc = {x_o[0], x_o[1], fabsf(b[0]), s_ang + sign * pi/2, e_ang + sign * pi/2};
                _v.push_back(val);
            }
        }
    }

    for(auto &x : _v){
        x.v_m = _vm;
        x.v_e = std::min(_ve, _vm);
    }
    (_v.end() - 1)->v_e = 0;

    return _v.size();
}


int make_trj(std::vector<fvector<3>> _xd, std::vector<xd_t> &_v, float _vm, float _ve){
    _v.clear();
    int num = _xd.size();
    for(int k = 1;k < num;k++){
        fvector<2> e(_xd[k][0], _xd[k][1]);

        xd_t val;

        val.xd = {e[0], e[1], _xd[k][2]};
        val.type = xd_t::Line;
        val.v_m = _vm;
        val.v_e = std::min(_ve, _vm);
        _v.push_back(val);
    }
    (_v.end() - 1)->v_e = 0;
    return _v.size();
}

fvector<3> ragrange(std::vector<xd_t> _xd_vec, float _u){

    int num = _xd_vec.size();
    float l[num];
    fvector<3> ret;

    _xd_vec.insert(_xd_vec.begin(), xd_t{});

    for(int i = 0;i < num + 1;i++){
        l[i] = 1;
        for(int j = 0;j < num + 1;j++){
            if(i == j) continue;
            l[i] *= (_u * num - j)/(i - j);
        }

        ret += l[i] * _xd_vec[i].xd;
    }

    return ret;
}

float drdu(std::vector<xd_t> _xd_vec, float _u){
    float du = 0.001;

    fvector<3> r = ragrange(_xd_vec, _u);
    fvector<3> dr = ragrange(_xd_vec, _u + du);

    return (r - dr).norm() / du;
}

float arg_rag(std::vector<xd_t> _xd_vec, float _u){
    float du = 0.001;

    fvector<3> r = ragrange(_xd_vec, _u);
    fvector<3> dr = ragrange(_xd_vec, _u + du);

    return atan2((dr - r)[1], (dr - r)[0]);
}

float v_max(float vm, float ddx, std::vector<xd_t> _xd_vec, float _u){
    float du = 0.001;

    float t = arg_rag(_xd_vec, _u);
    float dt = arg_rag(_xd_vec, _u + du);

    float r = 1000;
    if(fabsf(dt - t)/ du > 0.001) r = drdu(_xd_vec, _u) / fabsf(dt - t) * du;

    return std::min(vm, sqrtf(ddx * r));
}