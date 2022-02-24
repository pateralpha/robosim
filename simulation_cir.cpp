#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>

#include "la.h"
#include "circle.h"
#include "wheel.h"

using std::array;


inline const float adj_pi(float _t){
    while(_t > pi) _t -= 2 * pi;
    while(_t < - pi) _t += 2 * pi;
    return _t;
}
constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }

constexpr float r = 0.05; // radius of the wheel
constexpr float R = 0.5; // radius of the machine / distance between CoG and wheels

wheel w1(0.0083 * 14, 0.22, 0.33, 0.8);
wheel w2(0.0090 * 14, 0.25, 0.25, 0.8);
wheel w3(0.0075 * 14, 0.22, 0.40, 0.7);

auto M = diag(fvector<3>(15, 15, inertia_of_cylinder(15, R)));
auto I = diag(fvector<3>(inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r)));

fmatrix<3, 3> Juw = fmatrix<3, 3>{
        cosf(2 * pi / 3), sinf(2 * pi / 3), R,
        cosf(4 * pi / 3), sinf(4 * pi / 3), R,
        cosf(0), sinf(0), R
} / r;

fmatrix<3, 3> Rux(float _theta){
    return fmatrix<3, 3>{
            cosf(_theta), - sinf(_theta), 0, 
            sinf(_theta), cosf(_theta), 0,
            0, 0, 1
    };
}

fmatrix<3, 3> dJux(float _theta, float _dtheta){
    return fmatrix<3, 3>{
            - sinf(_theta) * _dtheta, - cosf(_theta) * _dtheta, 0, 
            cosf(_theta) * _dtheta, - sinf(_theta) * _dtheta, 0,
            0, 0, 0
    };
}

fmatrix<3, 3> J(float _theta){ return Rux(_theta) * Juw.inv(); }
fmatrix<3, 3> dJ(float _theta, float _dtheta){ return dJux(_theta, _dtheta) * Juw.inv(); }

fvector<6> f(fvector<6> _x, fvector<3> _e){
    fvector<3> x = top_v<6, 3>(_x), dx = bottom_v<6, 3>(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> w = J(theta).inv() * dx;
    fvector<3> tau(w1.torque(w[0], _e[0]), w2.torque(w[1], _e[1]), w3.torque(w[2], _e[2]));

    fvector<3> ddx = (M + J(theta).trans().inv() * I * J(theta).inv()).inv() * J(theta).trans().inv() * (tau - I * dJ(theta, dtheta) * dx);

    return merge_v(dx, ddx);
}

fvector<3> du(fvector<3> _dx, fvector<3> _u){
    return Rux(_u[2]).inv() * _dx;
}

void state(fvector<6> _x, fvector<3> &_tau, fvector<3> &_i, fvector<3> _e){
    fvector<3> x = top_v<6, 3>(_x), dx = bottom_v<6, 3>(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> w = J(theta).inv() * dx;

    _tau = fvector<3>(w1.torque(w[0], _e[0]), w2.torque(w[1], _e[1]), w3.torque(w[2], _e[2]));
    _i = fvector<3>(w1.current(w[0], _e[0]), w2.current(w[1], _e[1]), w3.current(w[2], _e[2]));
}


constexpr int N = 700;
constexpr float dt = 0.01;
constexpr int skip = 5;

// #define POS_PD
#define VEL_PD
// #define PURE_PURSUIT

#ifdef POS_PD
    auto Kp = diag<3>(1, 5, 2);
    auto Kd = diag<3>(0.3, 1.5, 0.6);
#elif defined VEL_PD
    auto Kp = diag<3>(1, 1, 1);
    auto Kd = diag<3>(0.4, 0.4, 0.4);
#elif defined PURE_PURSUIT
    auto Kp = diag<3>(1, 5, 2);
    auto Kd = diag<3>(0.3, 1.5, 0.6);
#endif

struct output{

    float t;

    fvector<3> x, dx;
    fvector<3> xd, dxd;
    fvector<3> x_cal, dx_cal;

    fvector<3> e, tau, i;

} out_v[N];

output log(float _t, fvector<6> _x, fvector<3> _xd, fvector<3> _dxd, 
        fvector<3> _x_cal, fvector<3> _dx_cal, 
        fvector<3> _e, fvector<3> _tau, fvector<3> _i){
    output ret;

    state(_x, ret.tau, ret.i, _e);

    ret.t = _t;
    ret.x = top_v<6, 3>(_x);
    ret.dx = bottom_v<6, 3>(_x);
    ret.xd = _xd;
    ret.dxd = _dxd;
    ret.x_cal = _x_cal;
    ret.dx_cal = _dx_cal;
    ret.e = _e;

    return ret;
}

template <typename V>
V rk4(V _v, std::function<V(V)> _f, float _dt){
    V k1 = _f(_v);
    V k2 = _f(_v + k1 *(_dt / 2));
    V k3 = _f(_v + k2 *(_dt / 2));
    V k4 = _f(_v + k3 * _dt);
    return _v +(_dt / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
}

int main(void){

    char fname[] = "a.dat";
    char fname2[] = "b.dat";
    char fname3[] = "c.dat";

    fvector<3> machine_pos[6];
    fvector<3> tire_pos[3][4];
    fmatrix<3, 3> Rr = Rux(2*pi/3), Rl = Rux(-2*pi/3);

    machine_pos[4] = {-0.08, -0.45, 0};
    machine_pos[5] = {0.08, -0.45, 0};

    machine_pos[0] = Rr * machine_pos[4];
    machine_pos[1] = Rr * machine_pos[5];

    machine_pos[2] = Rl * machine_pos[4];
    machine_pos[3] = Rl * machine_pos[5];

    tire_pos[2][0] = {-0.05, -0.46, 0};
    tire_pos[2][1] = { 0.05, -0.46, 0};
    tire_pos[2][2] = { 0.05, -0.49, 0};
    tire_pos[2][3] = {-0.05, -0.49, 0};

    for(int k = 0;k < 4;k++){
        tire_pos[0][k] = Rr * tire_pos[2][k];
        tire_pos[1][k] = Rl * tire_pos[2][k];
    }


    constexpr float max_voltage = 18; // muximum voltage
    constexpr float sat_time = 0.2; // minimum time from 0 V to maximum
    constexpr float de = max_voltage / sat_time * dt;

    constexpr float m_per_pulse = 0.2355e-3;

    /* parameters of the controller */
    constexpr float ke = 0.009;
    constexpr float eps = 0;

    constexpr float ddx = 5;
    constexpr float ddx_slow = 3;
    constexpr float ddx_buffer = 0.3;
    constexpr float ddx_rel = ddx * 1.2;

    constexpr float dx_d_max = 5;
    /* */

    /* list of desired points */
    // std::vector<fvector<3>> xd_vec = {
    //     {0, 0, 0},
    //     {2, 1, pi/4},
    //     {2, 4, pi/4},
    //     {1, 4, 0},
    // };

    // std::vector<fvector<3>> xd_vec = {
    //     {0, 0, 0}.
    //     {0.5, 0.15, pi/4},
    //     {1, 0.4, pi/4},
    //     {1.5, 0.75, pi/4},
    //     {2, 1.2, pi/4},
    //     {2.2, 1.5, pi/4},
    //     {2.3, 2, pi/4},
    //     {2.3, 2.5, pi/4},
    //     {2.2, 3, pi/4},
    //     {2, 3.3, pi/4},
    //     {1.5, 3.7, pi/4},
    //     {1, 4, 0},
    // };

    // std::vector<fvector<4>> xd_vec = {
    //     {0, 0, 0, 0},
    //     {2, 1.5, 0, pi/2},
    //     {2, 2.5, 0, pi/2},
    //     {1, 4, 0, 3*pi/4}
    // };

    std::vector<fvector<4>> xd_vec = {
        {0, 0, 0, 0},
        {2, 1, pi/4, pi/2},
        {2, 2.5, pi/4, pi/2},
        {3, 4, 0, pi/4}
    };

    // std::vector<fvector<4>> xd_vec = {
    //     {0, 0, 0, 0},
    //     {1, 0, 0, 0},
    //     {3, 2, 0, pi/2},
    //     {1, 4, 0, pi},
    //     {-1, 2, 0, -pi/2},
    //     {1, 0, 0, 0}
    // };

    std::vector<xd_t> xdc_arr;

    int n_xd = make_trj(xd_vec, xdc_arr);
    printf("Control point maked:: %d\r\n", n_xd);

    for(auto &v : xdc_arr){
        v.v_m = dx_d_max;
        v.v_e = dx_d_max;
        v.xd.print();
    }
    (xdc_arr.end() - 1)->v_e = 0;

    fvector<3> xd;
    /* */


    int index_xd = 0;

    /* variables for circular following */
    xd_t::type_t type;
    fvector<5> xd_arc;
    float arc_dir;
    /* */

    /* velocity control */
    fvector<3> dx_d, p_dx_d;
    float vm, ve;
    /* */

    /* dead reckoning */
    int enc_u = 0, enc_v = 0;
    int p_enc_u = enc_u, p_enc_v = enc_v;
    fvector<3> delta_u;

    fvector<3> x_cal = xd, p_x_cal;
    fvector<2> x_cal_pol;

    fvector<3> dx_cal, p_dx_cal;
    /* */
    
    /* rotational coordinate */
    fvector<3> x_org;
    fvector<3> uv, duv;
    float arg_uv;
    /* */

    /* variables for pure-pursuit */
    fvector<2> xd_tmp, p_xd_tmp;
    float arg, vel_tmp = 0;
    /* */

    bool over = true;
    fvector<3> input;


    fvector<6> v;
    fvector<3> x, dx;
    fvector<3> e, pe;

    fvector<3> u;

    fvector<3> tau, i;




    float t = 0;

    /* simulation for loop starts here */
    for(int n = 0;n < N;n++, t += dt){

        x = top_v<6, 3>(v);
        dx = bottom_v<6, 3>(v);

#ifndef PURE_PURSUIT
        out_v[n] = log(t, v, xd, dx_d, x_cal, dx_cal, e, tau, i);
#else
        out_v[n] = log(t, v, fvector<3>{xd_tmp[0], xd_tmp[1], xd[2]}, dx_d, x_cal, dx_cal, e, tau, i);
#endif

        /**** controller ****/
        /* dead reckoning */
        enc_u = u[0] / m_per_pulse;
        enc_v = u[1] / m_per_pulse;

        delta_u = fvector<3>(
                (enc_u - p_enc_u)* m_per_pulse, 
                (enc_v - p_enc_v)* m_per_pulse, 
                0
        );

        p_x_cal = x_cal;

        x_cal[2] = x[2];
        x_cal = x_cal + Rux(x_cal[2]) * delta_u;

        p_enc_u = enc_u;
        p_enc_v = enc_v;

        
        p_dx_cal = dx_cal;
        dx_cal = 0.6 * dx_cal + 0.4 * (x_cal - p_x_cal)/ dt;

        /* control input calculation */
        // check convergence and update xd if x_cal converged to xd
        if((xd - x_cal).norm() <= eps || over){

            over = false;
            if(index_xd < n_xd){
                if(xd != xdc_arr[index_xd].xd) x_org = xd;
                
                type = xdc_arr[index_xd].type;
                xd_arc = xdc_arr[index_xd].arc;
                vm = xdc_arr[index_xd].v_m;
                ve = xdc_arr[index_xd].v_e;
                xd = xdc_arr[index_xd++].xd;
                if(vm > dx_d_max) vm = dx_d_max;
                if(ve > dx_d_max) ve = dx_d_max;

                arg = xd_arc[3];

                if(type == xd_t::Arc){
                    vm = sqrtf(xd_arc[2] * ddx);
                }

                if(index_xd < n_xd){
                    if(xdc_arr[index_xd].type == xd_t::Arc){
                        ve = std::min(sqrtf(xdc_arr[index_xd].arc[2] * ddx), ve);
                    }
                }

                if((x_org[0] != xd[0])||(x_org[1] != xd[1])){
                    arg_uv = atan2(xd[1] - x_org[1], xd[0] - x_org[0]);
                }
                else arg_uv = 0;
            }
        }

        float dist;

        if(type == xd_t::Line){
            uv = Rux(arg_uv).trans()*(x_cal - xd);
            duv = Rux(arg_uv).trans()* dx_cal;

            arc_dir = - 1;
            dist = (Rux(arg_uv).trans()*(merge_v(xd_tmp, 0) - xd)).abs()[0];
        }
        else {
            fvector<2> x_cal_tmp = top_v<3, 2>(x_cal) - top_v<5, 2>(xd_arc);
            x_cal_pol = {x_cal_tmp.norm(), atan2(x_cal_tmp[1], x_cal_tmp[0])};

            arc_dir = xd_arc[4] > xd_arc[3] ? -1 : 1;
            arg_uv = ((arc_dir*(x_cal_pol[1] - xd_arc[3]))> 0 ? xd_arc[3] : x_cal_pol[1]) - pi / 2 * arc_dir;
            arg_uv = adj_pi(arg_uv);

            float x_t_diff = adj_pi(- arc_dir *(x_cal_pol[1] - xd_arc[4]));
            
            uv = fvector<3>{xd_arc[2]*(x_t_diff), arc_dir*(x_cal_pol[0] - xd_arc[2]), x_cal[2] - xd[2]};
            duv = Rux(arg_uv).trans()* dx_cal;

            dist = fabsf(xd_arc[2]*x_t_diff);
        }

#ifdef POS_PD
        /* PD feedback in uv-coordinate */
        duv += Rux(arg_uv).trans() * delta_dx_d;

        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Rux(arg_uv) * input;

#elif defined VEL_PD

        fvector<3> duv_d;

        if(fabsf(uv[0])< ddx_buffer) duv_d[0] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[0])+ ve * ve), uv[0]);
        else duv_d[0] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[0]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer + ve * ve), uv[0]);
        if(fabsf(uv[1])< ddx_buffer) duv_d[1] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[1])), uv[1]);
        else duv_d[1] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[1]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer), uv[1]);
        duv_d[2] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[2])), uv[2]);

        duv_d = duv_d.each([=](float _f){ return fabsf(_f) > vm ? std::copysign(vm, _f) : _f; });
        if(duv_d.norm() > vm) duv_d *= vm / duv_d.norm();

        p_dx_d = dx_d;
        dx_d = Rux(arg_uv) * duv_d;

        fvector<3> delta_dx = dx_d - p_dx_d;
        float dx_max = delta_dx.abs().max();

        for(int i = 0;i < 3;i++){
            if(dx_max > ddx_rel * dt){
                delta_dx[i] = ddx_rel * dt * delta_dx[i] / dx_max;
                if(fabsf(p_dx_d[i] - dx_d[i])> fabsf(delta_dx[i])) dx_d[i] = p_dx_d[i] + delta_dx[i];
            }
        }

        input = ke * 14 * J(x_cal[2]).inv()* dx_d
                + Kp * J(x_cal[2]).inv()*(dx_d - dx_cal)
                + Kd * J(x_cal[2]).inv()*(((dx_d - p_dx_d) - (dx_cal - p_dx_cal))/ dt - dJ(x_cal[2], dx_cal[2])* J(x_cal[2])*(dx_d - dx_cal));
        e = input;

        if(uv[0] > 0) over = true;

#elif defined PURE_PURSUIT

        fvector<2> xd_2dim = top_v<3, 2>(xd);

        vel_tmp = std::min(sqrtf(2 * ddx * dist + ve * ve), vel_tmp + ddx * dt);
        vel_tmp = std::min(vel_tmp, vm);

        p_xd_tmp = xd_tmp;

        if(type == xd_t::Line){
            if((xd_2dim - xd_tmp).norm()){
                xd_tmp += (xd_2dim - xd_tmp)/(xd_2dim - xd_tmp).norm()* vel_tmp * dt;
            }
        }
        else {
            arg += std::copysign(vel_tmp * dt / xd_arc[2], (xd_arc[4] - xd_arc[3]));
            xd_tmp = fvector<2>{xd_arc[0], xd_arc[1]} + xd_arc[2] * fvector<2>{cosf(arg), sinf(arg)};
        }
        
        if((p_xd_tmp - xd_2dim)*(xd_tmp - xd_2dim)<= 0){
            xd_tmp = xd_2dim;
            over = true;
        }

        dx_d = merge_v((xd_tmp - p_xd_tmp)/ dt, 0);

        uv = Rux(arg_uv).trans()*(x_cal - merge_v(xd_tmp, xd[2]));
        duv = Rux(arg_uv).trans()*(dx_cal - dx_d);

        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Rux(arg_uv) * input;


#endif

        float e_max = e.abs().max();
        if(e_max > max_voltage){
            for(int i = 0;i < 3;i++) e[i] *= max_voltage / e_max;
        }

        fvector<3> diff_e = e - pe;
        float de_max = diff_e.abs().max();

        for(int i = 0;i < 3;i++){
            diff_e[i] = de * diff_e[i] / de_max;
            if(fabsf(pe[i] - e[i])>fabsf(diff_e[i]))e[i] = pe[i] + diff_e[i];
        }
        pe = e;

        /**** end of calc. for controller ****/



        v = rk4<fvector<6>>(v, [=](fvector<6> _v){ return f(_v, e); }, dt);
        u = rk4<fvector<3>>(u, [=](fvector<3> _u){ return du(dx, _u); }, dt);

    }

    /* simulation completed */
    printf("calculation finish.\r\n");
    



    FILE *fp = fopen(fname, "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", fname);
    fprintf(fp, "# time x y theta e1 e2 e3 mx1 my1...\r\n");

    fvector<3> machine;
    fvector<3> tire;

    int cnt = 0;
    fvector<3> i_ave, i_peak;

    for(int n = 0;n < N;n += skip){
        fprintf(fp, "%5.2f %6.3f %6.3f %6.3f ",
                out_v[n].t, out_v[n].x[0], out_v[n].x[1], out_v[n].x[2]);
        fprintf(fp, "%6.3f %6.3f %6.3f ",
                out_v[n].dx[0], out_v[n].dx[1], out_v[n].dx[2]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].xd[0], out_v[n].xd[1]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].dxd[0], out_v[n].dxd[1]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].x_cal[0], out_v[n].x_cal[1]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].dx_cal[0], out_v[n].dx_cal[1]);
        fprintf(fp, "%6.2f %6.2f %6.2f ", out_v[n].e[0], out_v[n].e[1], out_v[n].e[2]);
        fprintf(fp, "%6.2f %6.2f %6.2f ", out_v[n].i[0], out_v[n].i[1], out_v[n].i[2]);

        for(int k = 0;k < 6;k++){
            machine = out_v[n].x + Rux(out_v[n].x[2]) * machine_pos[k];
            fprintf(fp, "%6.3f %6.3f ", machine[0], machine[1]);
        }
        for(int l = 0;l < 4*3;l++){
            tire = out_v[n].x + Rux(out_v[n].x[2]) * tire_pos[l / 4][l % 4];
            fprintf(fp, "%6.3f %6.3f ", tire[0], tire[1]);
        }
        fprintf(fp, "\r\n");

        cnt++;
        i_ave += out_v[n].i.abs();
        i_peak += (out_v[n].i.abs() - i_peak).each([](float _f){ return _f > 0 ? _f : 0; });
    }

    fclose(fp);

    fp = fopen(fname3, "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", fname3);
    fprintf(fp, "# control points\r\n");

    for(auto v : xdc_arr){
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
