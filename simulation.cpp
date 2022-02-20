#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>

#include "la.h"

using std::array;

constexpr float th_omega = 1e-2; // threshold value for determining stop or not
constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }

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

constexpr int N = 500;

struct output{
    float t;

    fvector<3> x;
    fvector<3> dx;

    fvector<3> xd;

    fvector<3> e;
    fvector<3> tau;
    fvector<3> i;
} out_v[N];

output log(float _t, fvector<6> _x, fvector<3> _xd, fvector<3> _e, fvector<3> _tau, fvector<3> _i){
    output ret;
    ret.t = _t;
    ret.x = upper(_x);
    ret.dx = lower(_x);
    ret.xd = _xd;
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

    constexpr float dt = 0.01;
    float t = 0;

    fvector<6> v;
    fvector<3> x, dx;
    fvector<3> e, pe;

    fvector<3> u;

    fvector<3> omega, tau, i;

    constexpr float max_voltage = 18; // muximum voltage
    constexpr float sat_time = 0.2; // minimum time from 0 V to maximum
    constexpr float de = 18 / sat_time * dt;

    constexpr float m_per_pulse = 0.2355e-3;

    auto Kp = diag<3>(1, 1, 2);
    auto Kd = diag<3>(0.2, 0.2, 0.4);
    constexpr float eps = 0.1;


    constexpr int n_xd = 1;
    // fvector<3> xd_arr[n_xd] = {{2, 1, pi/4}, {2, 3, pi/4}, {1, 4, 0}};
    fvector<3> xd_arr[n_xd] = {{3, 4, pi/4}};
    fvector<3> xd = xd_arr[0];

    int index_xd = 0;

    int enc_u = 0, enc_v = 0;
    int p_enc_u = enc_u, p_enc_v = enc_v;
    fvector<3> del_u;
    fvector<3> x_cal, p_x_cal;

    for(int n = 0;n < N;n++, t += dt){

        x = upper(v);
        dx = lower(v);

        out_v[n] = log(t, v, xd, e, tau, i);

        // controller

        // dead reckoning
        enc_u = u[0] / m_per_pulse;
        enc_v = u[1] / m_per_pulse;

        del_u = fvector<3>(
                (enc_u - p_enc_u)* m_per_pulse, 
                (enc_v - p_enc_v)* m_per_pulse, 
                0
        );

        p_x_cal = x_cal;

        x_cal[2] = x[2];
        x_cal = x_cal + Jux(x_cal[2]) * del_u;

        p_enc_u = enc_u;
        p_enc_v = enc_v;

        // control input

        // checkconvergence
        if((xd - x_cal).norm() < eps){
            if(index_xd < n_xd) xd = xd_arr[index_xd++];
        }

        e = J(x_cal[2]).inv()*(Kp *(xd - x_cal) - Kd * (x_cal - p_x_cal)*(1 / dt));

        for(int i = 0;i < 3;i++){
            if(fabsf(e[i]) > 18) e[i] = std::copysign(18, e[i]);
        }

        fvector<3> diff_e = e - pe;
        float de_max = diff_e.abs().max();

        for(int i = 0;i < 3;i++){
            diff_e[i] = de * diff_e[i] / de_max;
            if(pe[i] - e[i] > - diff_e[i])e[i] = pe[i] + diff_e[i];
            if(e[i] - pe[i] > diff_e[i])e[i] = pe[i] + diff_e[i];
        }
        pe = e;

        // end

        v = rk4<fvector<6>>(v, [=](fvector<6> _v){ return f(_v, e); }, dt);
        u = rk4<fvector<3>>(u, [=](fvector<3> _u){ return du(dx, _u); }, dt);

        state(v, tau, i, e);
    }

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

    constexpr int skip = 5;

    for(int n = 0;n < N;n += skip){
        fprintf(fp, "%5.2f %6.3f %6.3f %6.3f ",
                out_v[n].t, out_v[n].x[0], out_v[n].x[1], out_v[n].x[2]);
        fprintf(fp, "%6.3f %6.3f ", out_v[n].xd[0], out_v[n].xd[1]);
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
    }

    fclose(fp);

    fp = fopen(fname2, "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", fname2);
    fprintf(fp, "# control points\r\n");

    for(int n = 0;n < n_xd;n++){
        fprintf(fp, "%6.3f %6.3f %6.3f ", xd_arr[n][0], xd_arr[n][1], xd_arr[n][2]);
        fprintf(fp, "\r\n");
    }

    fclose(fp);
    printf("file closed.\r\n");

    return 0;
}
