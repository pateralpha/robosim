#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>
#include <vector>

#include "la.h"
#include "ip.h"
#include "machine.h"

using std::array;


inline const float adj_pi(float _t){
    while(_t > pi) _t -= 2 * pi;
    while(_t < - pi) _t += 2 * pi;
    return _t;
}

const wheel w1(0.0083 * 14, 0.22, 0.33, 0.8);
const wheel w2(0.0090 * 14, 0.25, 0.25, 0.8);
const wheel w3(0.0075 * 14, 0.22, 0.40, 0.7);

const fmatrix<3, 3> M = diag(fvector<3>(15, 15, inertia_of_cylinder(15, R)));
const fmatrix<3, 3> I = diag(fvector<3>(inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r)));

const fmatrix<3, 3> Juw = fmatrix<3, 3>{
        cosf(2 * pi / 3), sinf(2 * pi / 3), R,
        cosf(4 * pi / 3), sinf(4 * pi / 3), R,
        cosf(0), sinf(0), R
} / r;


void state(fvector<6> _x, fvector<3> &_tau, fvector<3> &_i, fvector<3> _e){
    fvector<3> x = top_v<6, 3>(_x), dx = bottom_v<6, 3>(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> w = J(theta).inv() * dx;

    _tau = fvector<3>(w1.torque(w[0], _e[0]), w2.torque(w[1], _e[1]), w3.torque(w[2], _e[2]));
    _i = fvector<3>(w1.current(w[0], _e[0]), w2.current(w[1], _e[1]), w3.current(w[2], _e[2]));
}

struct output{

    float t;

    fvector<6> x, xd, x_cal;
    fvector<3> e, tau, i;
};

output log(float _t, fvector<6> _x, fvector<6> _xd, 
        fvector<6> _x_cal, fvector<3> _e){
    
    output ret;
    state(_x, ret.tau, ret.i, _e);

    ret.t = _t;
    ret.x = _x;
    ret.xd = _xd;
    ret.x_cal = _x_cal;
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




constexpr int N = 500;
constexpr float dt = 0.01;
constexpr int skip = 5;

#define POS_PD
// #define VEL_PD
// #define PURE_PURSUIT

// #define LAGRANGE


#ifdef POS_PD
    auto Kp = diag<3>(1, 5, 2);
    auto Kd = diag<3>(0.3, 1.5, 0.6);
    float eps = 1;
#elif defined VEL_PD
    auto Kp = diag<3>(1, 1, 1);
    auto Kd = diag<3>(0.4, 0.4, 0.4);
    float eps = 0;
#elif defined PURE_PURSUIT
    auto Kp = diag<3>(1, 5, 2);
    auto Kd = diag<3>(0.3, 1.5, 0.6);
    float eps = 0;
#endif


void write(output out_v[], std::vector<xd_t> xdc_arr);

int main(void){

    constexpr float max_voltage = 18; // muximum voltage
    constexpr float sat_time = 0.2; // minimum time from 0 V to maximum
    constexpr float de = max_voltage / sat_time * dt;

    /* parameters of the controller */
    constexpr float ke = 0.01;

    constexpr float ddx = 5;
    constexpr float ddx_slow = 2;
    constexpr float ddx_rel = ddx * 1.5;
    constexpr float ddx_buffer = 0.3;

    constexpr float dx_max = 5;
    /* */

    /* list of desired points */
    // std::vector<fvector<3>> xd_vec = {
    //     {0, 0, 0},
    //     {2, 1, pi/4},
    //     {2, 4, pi/4},
    //     {1, 4, 0},
    // };

    std::vector<fvector<3>> xd_vec = {
        {0, 0, 0},
        {2, 1.5, pi/4},
        {2, 3, pi/4},
        {1, 4, 0},
    };

    // std::vector<fvector<3>> xd_vec = {
    //     {0, 0, 0},
    //     {2, 1.5, pi/4},
    //     {2, 3, pi/4},
    //     {3, 4, 0},
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
    //     {2, 1, pi/4, pi/2},
    //     {2, 2.5, pi/4, pi/2},
    //     {3, 4, 0, pi/4}
    // };

    // std::vector<fvector<4>> xd_vec = {
    //     {0, 0, 0, 0},
    //     {1, 0, 0, 0},
    //     {3, 2, 0, pi/2},
    //     {1, 4, 0, pi},
    //     {-1, 2, 0, -pi/2},
    //     {1, 0, 0, 0}
    // };

#ifdef VEL_PD
    if(typeid(xd_vec[0]) == typeid(fvector<3>)) eps = 1;
#endif

    std::vector<xd_t> xdc_arr;

#ifdef PURE_PURSUIT
    int n_xd = make_trj(xd_vec, xdc_arr, dx_max, 2);
#else
    int n_xd = make_trj(xd_vec, xdc_arr, dx_max);
#endif
    printf("Control point created:: %d\r\n", n_xd);

    for(const auto &x : xdc_arr) x.xd.print();

    xd_t xd = {}, p_xd;

    bool over = true;
    int index_xd = 0;
    float arc_dir = - 1;
    /* */

    /* velocity control */
    fvector<3> dx_d, p_dx_d;
    float dist;
    /* */

    /* dead reckoning */
    fvector<2> enc, p_enc;

    fvector<3> x_cal = xd.xd, p_x_cal;
    fvector<3> dx_cal, p_dx_cal;
    /* */
    
    /* the machine coordinate */
    fvector<3> uv, duv;
    float arg_uv;
    auto Rux = [&](){ return Rotr(arg_uv); };
    auto Rxu = [&](){ return Rux().trans(); };
    /* */

    /* variables for pure-pursuit */
    fvector<2> xd_tmp, p_xd_tmp;
    float arg_tmp, vel_tmp = 0;
    /* */

    /* lagrange's interpolation */
    float p = 0;
    /* */

    fvector<6> v;
    fvector<3> e, pe;

    fvector<3> u;
    output  out_v[N];



    float t = 0;

    /* simulation for loop starts here */
    for(int n = 0;n < N;n++, t += dt){

        fvector<3> x = top_v<6, 3>(v);
        fvector<3> dx = bottom_v<6, 3>(v);

#ifndef PURE_PURSUIT
        out_v[n] = log(t, v, merge_v(xd.xd, dx_d), merge_v(x_cal, dx_cal), e);
#else
        out_v[n] = log(t, v, merge_v(merge_v(xd_tmp, xd.xd[2]), dx_d), merge_v(x_cal, dx_cal), e);
#endif /* PURE_PURSUIT */

        /**** controller ****/
        /* dead reckoning */
        p_enc = enc;
        enc = (top_v<3, 2>(u) / m_per_pulse).each([](float _f){ return (int) _f; });

        fvector<3> delta_u = merge_v((enc - p_enc)* m_per_pulse, 0);

        p_x_cal = x_cal;
        x_cal = merge_v(top_v<3, 2>(x_cal + Rotr(x_cal[2]) * delta_u), x[2]);
        
        p_dx_cal = dx_cal;
        dx_cal = 0.6 * dx_cal + 0.4 * (x_cal - p_x_cal)/ dt;

        /* control input calculation */
        
#ifndef LAGRANGE
        // check convergence and update xd if x_cal converged to xd
        if((xd.xd - x_cal).norm() <= eps || over){

            over = false;
            if(index_xd < n_xd){

                p_xd = xd;
                xd = xdc_arr[index_xd++];
                if(xd.v_m > dx_max) xd.v_m = dx_max;
                if(xd.v_e > dx_max) xd.v_e = dx_max;

                arg_tmp = xd.arc[3];

                if(xd.type == xd_t::Arc){
                    xd.v_m = sqrtf(xd.arc[2] * ddx);
                }

                if(index_xd < n_xd){
                    if(xdc_arr[index_xd].type == xd_t::Arc){
                        xd.v_e = std::min(sqrtf(xdc_arr[index_xd].arc[2] * ddx), xd.v_e);
                    }
                }

                if(top_v<3, 2>(p_xd.xd - xd.xd).norm()){
                    arg_uv = atan2(xd.xd[1] - p_xd.xd[1], xd.xd[0] - p_xd.xd[0]);
                }
            }
        }

        // calcurate position in the local coordinate
        if(xd.type == xd_t::Line){
            uv = Rxu()*(x_cal - xd.xd);
            duv = Rxu()* dx_cal;

            arc_dir = - 1;
            dist = (xd_tmp - top_v<3, 2>(xd.xd)).norm();
        }
        else {
            fvector<2> x_cal_tmp = top_v<3, 2>(x_cal) - top_v<5, 2>(xd.arc);
            fvector<3> x_cal_pol = {x_cal_tmp.norm(), atan2(x_cal_tmp[1], x_cal_tmp[0])};

            arc_dir = xd.arc[4] > xd.arc[3] ? -1 : 1;
            arg_uv = ((arc_dir*(x_cal_pol[1] - xd.arc[3]))> 0 ? xd.arc[3] : x_cal_pol[1]) - pi / 2 * arc_dir;
            arg_uv = adj_pi(arg_uv);

            float x_t_diff = adj_pi(- arc_dir *(x_cal_pol[1] - xd.arc[4]));
            
            uv = fvector<3>{xd.arc[2] * x_t_diff, arc_dir*(x_cal_pol[0] - xd.arc[2]), x_cal[2] - xd.xd[2]};
            duv = Rxu()* dx_cal;

            dist = fabsf(xd.arc[2] * adj_pi(arg_tmp - xd.arc[4]));
        }

#else

        arg_uv = arg_rag(xdc_arr, p);

        xd.xd = ragrange(xdc_arr, p);
        uv = Rxu()*(x_cal - xd.xd);

        index_xd = p * n_xd;
        if(index_xd >= n_xd) index_xd = n_xd - 1;

#ifdef VEL_PD
        dist = (xdc_arr[index_xd].xd - xd.xd).norm();
#else
        dist = (top_v<3, 2>(xdc_arr[index_xd].xd) - xd_tmp).norm();
#endif /* VEL_PD */

        for(int i = index_xd + 1;i < n_xd;i++){
            dist += (xdc_arr[i].xd - xdc_arr[i - 1].xd).norm();
        }

        if(dist > 1)uv[0] = - dist;
        else uv[0] = (Rxu()*(x_cal - (xdc_arr.end() - 1)->xd))[0];

        xd.v_e = 0;
        xd.v_m = v_max(dx_max, ddx, xdc_arr, p);
#endif /* LAGRANGE */



#ifdef POS_PD
        /* PD feedback in uv-coordinate */
        fvector<3> input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Rux() * input;

#elif defined VEL_PD

#ifdef LAGRANGE

        if(drdu(xdc_arr, p) > 0.0001){
            p += (dx_cal.norm() * dt + (Rxu()*(x_cal - xd.xd))[0]) / drdu(xdc_arr, p);
        }
        else p += 0.0001;

        if(p < 0) p = 0;
        if(p > 1) p = 1;
#endif

        fvector<3> duv_d;

        if(fabsf(uv[0])< ddx_buffer) duv_d[0] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[0])+ xd.v_e * xd.v_e), uv[0]);
        else duv_d[0] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[0]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer + xd.v_e * xd.v_e), uv[0]);
        if(fabsf(uv[1])< ddx_buffer) duv_d[1] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[1])), uv[1]);
        else duv_d[1] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[1]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer), uv[1]);
        duv_d[2] = - std::copysign(sqrtf(2 * ddx * fabsf(uv[2])), uv[2]);

        duv_d = duv_d.each([=](float _f){ return std::copysign(std::min(fabsf(_f), xd.v_m), _f); });
        if(duv_d.norm() > xd.v_m) duv_d *= xd.v_m / duv_d.norm();

        p_dx_d = dx_d;
        dx_d = Rux() * duv_d;

        fvector<3> delta_dx = dx_d - p_dx_d;
        float dx_max = delta_dx.abs().norm();

        if(dx_max > ddx * dt){
            for(int i = 0;i < 3;i++){
                if(delta_dx[i] * dx_d[i] < 0) delta_dx[i] = ddx_rel * dt * delta_dx[i] / dx_max;
                else delta_dx[i] = ddx * dt * delta_dx[i] / dx_max;
            }
        }
        dx_d = p_dx_d + delta_dx;

        fvector<3> input = ke * 14 * J(x_cal[2]).inv()* dx_d
                + Kp * J(x_cal[2]).inv()*(dx_d - dx_cal)
                + Kd * J(x_cal[2]).inv()*(((dx_d - p_dx_d) - (dx_cal - p_dx_cal))/ dt - dJ(x_cal[2], dx_cal[2])* J(x_cal[2])*(dx_d - dx_cal));
        e = input;

        if(uv[0] > 0) over = true;

#elif defined PURE_PURSUIT

        fvector<2> xd_2dim = top_v<3, 2>(xd.xd);

        vel_tmp = std::min(sqrtf(2 * ddx * dist + xd.v_e * xd.v_e), vel_tmp + ddx * dt);
        vel_tmp = std::min(vel_tmp, xd.v_m);

        p_xd_tmp = xd_tmp;

#ifdef LAGRANGE

        if(drdu(xdc_arr, p) > 0.0001) p += vel_tmp * dt / drdu(xdc_arr, p);
        else p += 0.0001;

        if(p > 1)p = 1;
        if(p < 0)p = 0;
        

        xd_tmp = top_v<3, 2>(ragrange(xdc_arr, p));
        uv = Rxu()*(x_cal - ragrange(xdc_arr, p));
#else

        if(xd.type == xd_t::Line){
            if((xd_2dim - xd_tmp).norm()){
                xd_tmp += (xd_2dim - xd_tmp)/(xd_2dim - xd_tmp).norm()* vel_tmp * dt;
            }
        }
        else {
            arg_tmp += std::copysign(vel_tmp * dt / xd.arc[2], (xd.arc[4] - xd.arc[3]));
            xd_tmp = fvector<2>{xd.arc[0], xd.arc[1]} + xd.arc[2] * fvector<2>{cosf(arg_tmp), sinf(arg_tmp)};
        }

        
        if((p_xd_tmp - xd_2dim)*(xd_tmp - xd_2dim)<= 0){
            xd_tmp = xd_2dim;
            over = true;
        }

        uv = Rxu()*(x_cal - merge_v(xd_tmp, xd.xd[2]));
#endif /* LAGRANGE */

        dx_d = merge_v((xd_tmp - p_xd_tmp)/ dt, 0);
        duv = Rxu()*(dx_cal - dx_d);

        fvector<3> input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Rux() * input;


#endif /* POS_PD */

        float e_max = e.abs().max();
        if(e_max > max_voltage) e *= max_voltage / e_max;

        fvector<3> diff_e = e - pe;
        float de_max = diff_e.abs().max();

        for(int i = 0;i < 3;i++) diff_e[i] = de * diff_e[i] / de_max;

        e = pe + diff_e;
        pe = e;

        /**** end of calc. for controller ****/



        v = rk4<fvector<6>>(v, [=](fvector<6> _v){ return f(_v, e); }, dt);
        u = rk4<fvector<3>>(u, [=](fvector<3> _u){ return du(dx, _u); }, dt);

    }/* for */

    /* simulation completed */
    printf("calculation finish.\r\n");
    
    write(out_v, xdc_arr);

    return 0;
}





void write(output out_v[], std::vector<xd_t> xdc_arr){

    char fname[] = "a.dat";
    char fname2[] = "c.dat";

    fvector<3> machine_pos[6];
    fvector<3> tire_pos[3][4];
    fmatrix<3, 3> Rr = Rotr(2*pi/3), Rl = Rotr(-2*pi/3);

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
                out_v[n].x[3], out_v[n].x[4], out_v[n].x[5]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].xd[0], out_v[n].xd[1]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].xd[3], out_v[n].xd[4]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].x_cal[0], out_v[n].x_cal[1]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].x_cal[3], out_v[n].x_cal[4]);
        fprintf(fp, "%6.2f %6.2f %6.2f ", out_v[n].e[0], out_v[n].e[1], out_v[n].e[2]);
        fprintf(fp, "%6.2f %6.2f %6.2f ", out_v[n].i[0], out_v[n].i[1], out_v[n].i[2]);

        for(int k = 0;k < 6;k++){
            machine = top_v<6, 3>(out_v[n].x) + Rotr(out_v[n].x[2]) * machine_pos[k];
            fprintf(fp, "%6.3f %6.3f ", machine[0], machine[1]);
        }
        for(int l = 0;l < 4*3;l++){
            tire = top_v<6, 3>(out_v[n].x) + Rotr(out_v[n].x[2]) * tire_pos[l / 4][l % 4];
            fprintf(fp, "%6.3f %6.3f ", tire[0], tire[1]);
        }
        fprintf(fp, "\r\n");

        cnt++;
        i_ave += out_v[n].i.abs();
        i_peak += (out_v[n].i.abs() - i_peak).each([](float _f){ return _f > 0 ? _f : 0; });
    }

    fclose(fp);

    fp = fopen(fname2, "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", fname2);
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

}