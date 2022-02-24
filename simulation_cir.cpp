#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>

#include "la.h"
#include "circle.h"

using std::array;

constexpr float th_omega = 1e-2; // threshold value for determining stop or not
constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }


template <std::size_t N>
inline fvector<N + 1> merge_v(float _f, fvector<N> _v){
    fvector<N + 1> ret;
    ret[0] = _f;
    for(int i = 1;i <= N;i++) ret[i] = _v[i];
    return ret;
}

template <std::size_t N>
inline fvector<N + 1> merge_v(fvector<N> _v, float _f){
    fvector<N + 1> ret;
    for(int i = 0;i < N;i++) ret[i] = _v[i];
    ret[N] = _f;
    return ret;
}

template <std::size_t N, std::size_t M>
inline fvector<N + M> merge_v(fvector<N> _v1, fvector<M> _v2){
    fvector<N + 1> ret;
    for(int i = 0;i < N;i++) ret[i] = _v1[i];
    for(int i = N;i < N + M;i++) ret[i] = _v2[i];
    return ret;
}

template <std::size_t N, std::size_t M>
inline fvector<M> top_v(fvector<N> _v){
    fvector<M> ret;
    for(int i = 0;i < N && i < M;i++) ret[i] = _v[i];
    return ret;
}

class wheel{
public:
    wheel(
        float _k,   // toruque/reverse voltage coefficient
        float _R,   // internal resistance
        float _f0,  // maximum static friction torque
        float _gamma// transfer efficiency
    ) : k(_k), R(_R), f0(_f0), gamma(_gamma){}

    float torque(float _omega, float _E) const{
        float tau = torque_ac(_omega, _E); // actuation torque
        float tauf; // friction torque or starting torque

        if(fabsf(_omega)> th_omega) tauf = std::copysign(f0, _omega);
        else tauf = std::copysign(std::min(fabsf(tau), f0), tau);

        return tau - tauf;
    }

    float torque_ac(float _omega, float _E) const{
        float i = current(_omega, _E);
        return i * k * gamma;
    }

    float current(float _omega, float _E) const{
        return (_E - _omega*k)/R;
    }

private:
    float k, R;
    float f0, gamma;
};

constexpr float r = 0.05; // radius of the wheel
constexpr float R = 0.5; // radius of the machine / distance between CoG and wheels

wheel w1(0.0083 * 14, 0.22, 0.33, 0.8);
wheel w2(0.0090 * 14, 0.25, 0.25, 0.8);
wheel w3(0.0075 * 14, 0.22, 0.40, 0.7);

auto M = diag(fvector<3>(15, 15, inertia_of_cylinder(15, R)));
auto I = diag(fvector<3>(inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r)));

fmatrix<3, 3> Juw = fmatrix<3, 3>{
        cosf(2 * pi / 3)/ r, cosf(4 * pi / 3)/ r, cosf(0)/ r,
        sinf(2 * pi / 3)/ r, sinf(4 * pi / 3)/ r, sinf(0)/ r,
        R / r, R / r, R / r
}.trans();

fmatrix<3, 3> Jux(float _theta){
    return fmatrix<3, 3>{
            cosf(_theta), sinf(_theta), 0, 
            - sinf(_theta), cosf(_theta), 0,
            0, 0, 1
    }.trans();
}

fmatrix<3, 3> dJux(float _theta, float _dtheta){
    return fmatrix<3, 3>{
            - sinf(_theta) * _dtheta, cosf(_theta) * _dtheta, 0, 
            - cosf(_theta) * _dtheta, - sinf(_theta) * _dtheta, 0,
            0, 0, 0
    }.trans();
}

fmatrix<3, 3> J(float _theta){
    return Jux(_theta) * Juw.inv();
}

fmatrix<3, 3> dJ(float _theta, float _dtheta){
    return dJux(_theta, _dtheta) * Juw.inv();
}

const fvector<3> upper(const fvector<6> &_v){ return fvector<3>(_v[0], _v[1], _v[2]); }
const fvector<3> lower(const fvector<6> &_v){ return fvector<3>(_v[3], _v[4], _v[5]); }

fvector<6> f(fvector<6> _x, fvector<3> _e){
    fvector<3> x = upper(_x), dx = lower(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> omega = J(theta).inv() * dx;

    fvector<3> tau(
            w1.torque(omega[0], _e[0]), 
            w2.torque(omega[1], _e[1]), 
            w3.torque(omega[2], _e[2])
    );

    fvector<3> ddx = (M + J(theta).trans().inv() * I * J(theta).inv()).inv() * J(theta).trans().inv() * (tau - I * dJ(theta, dtheta) * dx);

    return fvector<6>(dx[0], dx[1], dx[2], ddx[0], ddx[1], ddx[2]);
}

fvector<3> du(fvector<3> _dx, fvector<3> _u){
    return Jux(_u[2]).inv() * _dx;
}

void state(fvector<6> &_x, fvector<3> &_tau, fvector<3> &_i, fvector<3> _e){
    fvector<3> x = upper(_x), dx = lower(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> omega = J(theta).inv() * dx;

    _tau = fvector<3>(
            w1.torque(omega[0], _e[0]), 
            w2.torque(omega[1], _e[1]), 
            w3.torque(omega[2], _e[2])
    );
    _i = fvector<3>(
            w1.current(omega[0], _e[0]), 
            w2.current(omega[1], _e[1]),
            w3.current(omega[2], _e[2])
    );
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

    fvector<3> x;
    fvector<3> dx;

    fvector<3> xd;
    fvector<3> dxd;

    fvector<3> x_cal;
    fvector<3> dx_cal;

    fvector<3> e;
    fvector<3> tau;
    fvector<3> i;
} out_v[N];

output log(float _t, fvector<6> _x, fvector<3> _xd, fvector<3> _dxd, 
        fvector<3> _x_cal, fvector<3> _dx_cal, 
        fvector<3> _e, fvector<3> _tau, fvector<3> _i){
    output ret;
    ret.t = _t;
    ret.x = upper(_x);
    ret.dx = lower(_x);
    ret.xd = _xd;
    ret.dxd = _dxd;
    ret.x_cal = _x_cal;
    ret.dx_cal = _dx_cal;
    ret.e = _e;
    ret.tau = _tau;
    ret.i = _i;
    return ret;
}

template <typename v>
v rk4(v _v, std::function<v(v)> _f, float _dt){
    v k1 = _f(_v);
    v k2 = _f(_v + k1 *(_dt / 2));
    v k3 = _f(_v + k2 *(_dt / 2));
    v k4 = _f(_v + k3 * _dt);

    return _v +(_dt / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
}

int main(void){

    char fname[] = "a.dat";
    char fname2[] = "b.dat";
    char fname3[] = "c.dat";

    fvector<3> machine_pos[6];
    fvector<3> tire_pos[3][4];
    fmatrix<3, 3> Rr = Jux(2*pi/3), Rl = Jux(-2*pi/3);

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

    float t = 0;

    fvector<6> v;
    fvector<3> x, dx;
    fvector<3> e, pe;

    fvector<3> u;

    fvector<3> omega, tau, i;

    constexpr float max_voltage = 18; // muximum voltage
    constexpr float sat_time = 0.2; // minimum time from 0 V to maximum
    constexpr float de = max_voltage / sat_time * dt;
    constexpr float ke = 0.009;

    constexpr float m_per_pulse = 0.2355e-3;


    constexpr float eps = 0;

    constexpr float ddx = 5;
    constexpr float dx_d_max = 5;
    constexpr float ddx_slow = 3;
    constexpr float ddx_buffer = 0.3;
    constexpr float ddx_rel = ddx * 1.2;

    float ddw_rel = (J(0).inv()* fvector<3>{1, 0, 0} * ddx_rel).abs().max();

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

    std::vector<fvector<4>> xd_vec = {
        {0, 0, 0, 0},
        {2, 1.5, 0, pi/2},
        {2, 2.5, 0, pi/2},
        {1, 4, 0, 3*pi/4}
    };

    // std::vector<fvector<4>> xd_vec = {
    //     {0, 0, 0, 0},
    //     {2, 1, pi/4, pi/2},
    //     {2, 2.5, pi/4, pi/2},
    //     {1, 4, 0, 3*pi/4}
    // };

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

    fvector<3> xd = xdc_arr[0].xd, p_xd;

    int type = xdc_arr[0].type;
    fvector<5> xd_arc = xdc_arr[0].arc;
    float arc_dir = xdc_arr[0].arc[4] > xdc_arr[0].arc[3] ? -1 : 1;
    float vm = xdc_arr[0].v_m, ve = xdc_arr[0].v_e;
    

    int index_xd = 1;

    int enc_u = 0, enc_v = 0;
    int p_enc_u = enc_u, p_enc_v = enc_v;
    fvector<3> del_u;
    fvector<3> x_cal, p_x_cal;

    fvector<3> x_org = upper(v);
    fvector<3> dx_cal, p_dx_cal, dx_d, p_dx_d, duv_d, dw_d, p_dw_d;
    fvector<3> delta_dx_d;

    fvector<2> xd_tmp, p_xd_tmp;
    float arg = xd_arc[3];
    float vel_tmp = 0;
    float vel = 0;
    
    fvector<3> uv, duv;
    float arg_uv;

    if((x_org[0] != xd[0])||(x_org[1] != xd[1])){
        if(type == xd_t::Line){
            arg_uv = atan2(xd[1] - x_org[1], xd[0] - x_org[0]);
        }
        else if(type == xd_t::Arc){
            arg_uv = xd_arc[3];
        }
    }
    else arg_uv = 0;

    if(index_xd < n_xd){
        if(xdc_arr[index_xd].type == xd_t::Arc){
            ve = std::min(sqrtf(xdc_arr[index_xd].arc[2] * ddx), ve);
        }
    }

    bool over = false;

    fvector<3> input;
    fvector<2> x_cal_pol;




    /* simulation for loop starts here */
    for(int n = 0;n < N;n++, t += dt){

        x = upper(v);
        dx = lower(v);

#ifndef PURE_PURSUIT
        out_v[n] = log(t, v, xd, dx_d, x_cal, dx_cal, e, tau, i);
#else
        out_v[n] = log(t, v, fvector<3>{xd_tmp[0], xd_tmp[1], xd[2]}, dx_d, x_cal, dx_cal, e, tau, i);
#endif

        /**** controller ****/
        /* dead reckoning */
        enc_u = u[0] / m_per_pulse;
        enc_v = u[1] / m_per_pulse;

        del_u = fvector<3>(
                (enc_u - p_enc_u)* m_per_pulse, 
                (enc_v - p_enc_v)* m_per_pulse, 
                0
        );

        p_x_cal = x_cal;
        p_xd = xd;

        x_cal[2] = x[2];
        x_cal = x_cal + Jux(x_cal[2]) * del_u;

        p_enc_u = enc_u;
        p_enc_v = enc_v;

        
        p_dx_cal = dx_cal;
        dx_cal = 0.6 * dx_cal + 0.4 * (x_cal - p_x_cal)/ dt;

        /* control input calculation */
        // check convergence and update xd if x_cal converged to xd
        if((xd - x_cal).norm() < eps || over){

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

        if(type == xd_t::Arc){
            fvector<2> xd_arc_org = {xd_arc[0], xd_arc[1]};
            fvector<2> x_cal_tmp = fvector<2>{x_cal[0], x_cal[1]} - xd_arc_org;
            x_cal_pol = {x_cal_tmp.norm(), atan2(x_cal_tmp[1], x_cal_tmp[0])};
            arc_dir = xd_arc[4] > xd_arc[3] ? -1 : 1;
            arg_uv = ((arc_dir*(x_cal_pol[1] - xd_arc[3]))> 0 ? xd_arc[3] : x_cal_pol[1]) - pi / 2 * arc_dir;
            while(arg_uv > pi)arg_uv -= 2*pi;
            while(arg_uv < -pi)arg_uv += 2*pi;
        }
        else arc_dir = - 1;

        if(type == xd_t::Line){
        // || (xd - x_cal).norm() < 0.1){
            uv = Jux(arg_uv).trans()*(x_cal - xd);
            duv = Jux(arg_uv).trans()*(dx_cal - (xd - p_xd)/ dt);
        }
        else {
            float x_t_diff = x_cal_pol[1] - xd_arc[4];
            while(x_t_diff > pi)x_t_diff -= 2*pi;
            while(x_t_diff < -pi)x_t_diff += 2*pi;
            uv = fvector<3>{xd_arc[2]*(x_t_diff), arc_dir*(x_cal_pol[0] - xd_arc[2]), x_cal[2] - xd[2]};
            duv = Jux(arg_uv).trans()*(dx_cal - (xd - p_xd)/ dt);
        }

        float dist;

        if(type == xd_t::Line){
            dist = (Jux(arg_uv).trans()*(merge_v(xd_tmp, 0) - xd)).abs()[0];
        }
        else {
            float x_t_diff = arg - xd_arc[4];
            dist = fabsf(xd_arc[2]*x_t_diff);
        }

        // vel = (Jux(arg_uv).trans()* dx_cal)[0];
        // if(type == xd_t::Arc){
        //     delta_dx_d = - arc_dir * vel * vel * dt * fvector<3>{- sinf(arg_uv), cosf(arg_uv), 0} / xd_arc[2];
        // }
        // else delta_dx_d = {0, 0, 0};

#ifdef POS_PD
        /* PD feedback in uv-coordinate */
        duv += Jux(arg_uv).trans() * delta_dx_d;

        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Jux(arg_uv) * input;

#elif defined VEL_PD

        if(fabsf(uv[0])< ddx_buffer) duv_d[0] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[0])+ ve * ve), uv[0]);
        else duv_d[0] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[0]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer + ve * ve), uv[0]);
        if(fabsf(uv[1])< ddx_buffer) duv_d[1] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[1])), uv[1]);
        else duv_d[1] = - std::copysign(sqrtf(2 * ddx * (fabsf(uv[1]) - ddx_buffer)+ 2 * ddx_slow * ddx_buffer), uv[1]);
        duv_d[2] = - std::copysign(sqrtf(2 * ddx_slow * fabsf(uv[2])), uv[2]);

        duv_d = duv_d.each([=](float _f){
            return fabsf(_f) > vm ? std::copysign(vm, _f) : _f;
        });
        if(duv_d.norm() > vm) duv_d *= vm / duv_d.norm();

        p_dx_d = dx_d;
        dx_d = Jux(arg_uv) * duv_d + delta_dx_d;

        fvector<3> delta_dx = dx_d - p_dx_d;
        float dx_max = delta_dx.abs().max();

        for(int i = 0;i < 3;i++){
            if(dx_max > ddx_rel * dt){
                delta_dx[i] = ddx_rel * dt * delta_dx[i] / dx_max;
                if(fabsf(p_dx_d[i] - dx_d[i])> fabsf(delta_dx[i])) dx_d[i] = p_dx_d[i] + delta_dx[i];
            }
        }
        
        p_dw_d = dw_d;
        dw_d = J(x_cal[2]).inv()* dx_d;
        fvector<3> ddw_d = J(x_cal[2]).inv()* ((dx_d - p_dx_d)/ dt - dJ(x_cal[2], dx_cal[2])* dw_d);
        fvector<3> dw = J(x_cal[2]).inv()* dx_cal;
        fvector<3> ddw = J(x_cal[2]).inv()* ((dx_cal - p_dx_cal)/ dt - dJ(x_cal[2], dx_cal[2])* dw_d);

        input = ke * 14 * dw_d + Kp *(dw_d - dw) + Kd *(ddw_d - ddw);
        e = input;

        if(- uv[0] * arc_dir > 0) over = true;

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

        uv = Jux(arg_uv).trans()*(x_cal - merge_v(xd_tmp, xd[2]));
        duv = Jux(arg_uv).trans()*(dx_cal - dx_d);

        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Jux(arg_uv) * input;


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

        state(v, tau, i, e);
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
            machine = out_v[n].x + Jux(out_v[n].x[2]) * machine_pos[k];
            fprintf(fp, "%6.3f %6.3f ", machine[0], machine[1]);
        }
        for(int l = 0;l < 4*3;l++){
            tire = out_v[n].x + Jux(out_v[n].x[2]) * tire_pos[l / 4][l % 4];
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
