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

constexpr int N = 700;
constexpr float dt = 0.01;
constexpr int skip = 5;

#define POS_PD
// #define VEL_PD
// #define PURE_PURSUIT

#ifdef POS_PD
    auto Kp = diag<3>(1, 5, 2);
    auto Kd = diag<3>(0.3, 1.5, 0.6);
#elif defined VEL_PD
    auto Kp = diag<3>(1, 1, 1);
    auto Kd = diag<3>(0.3, 0.3, 0.3);
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


    constexpr float eps = 0.5;

    constexpr float ddx = 3;
    constexpr float dx_d_max = 5;
    constexpr float ddx_slow = 3;
    constexpr float ddx_buffer = 0.3;

    // fvector<3> xd_arr[] = {{2, 1, pi/4}, {2, 4, pi/4}, {1, 4, 0}};
    // fvector<3> xd_arr[] = {{3, 4, pi/4}};
    fvector<3> xd_arr[] = {
        {0.5, 0.15, pi/4},
        {1, 0.4, pi/4},
        {1.5, 0.75, pi/4},
        {2, 1.2, pi/4},
        {2.2, 1.5, pi/4},
        {2.3, 2, pi/4},
        {2.3, 2.5, pi/4},
        {2.2, 3, pi/4},
        {2, 3.3, pi/4},
        {1.5, 3.7, pi/4},
        {1, 4, 0},
    };
    // int n_xd = sizeof(xd_arr) / sizeof(xd_arr[0]);
    // fvector<3> xd = xd_arr[0], p_xd;

    xd_t xdc_arr[] = {
        {{0.5, 0, 0}, xd_t::Line, {}},
        {{2, 1.5, 0}, xd_t::Arc, {0.5, 1.5, 1.5, -1.57, 0}},
        {{2, 2.5, 0}, xd_t::Line, {}},
        {{1.65, 3.35, 0}, xd_t::Arc, {0.79, 2.5, 1.21, 0, 1.57}},
        {{1, 4, 0}, xd_t::Line, {}},
    };
    int n_xd = sizeof(xdc_arr) / sizeof(xdc_arr[0]);
    fvector<3> xd = xdc_arr[0].xd, p_xd;

    int type = xdc_arr[0].type;
    fvector<5> xd_arc = xdc_arr[0].arc;
    float arc_dir = xdc_arr[0].arc[4] > xdc_arr[0].arc[3] ? -1 : 1;
    

    int index_xd = 1;

    int enc_u = 0, enc_v = 0;
    int p_enc_u = enc_u, p_enc_v = enc_v;
    fvector<3> del_u;
    fvector<3> x_cal, p_x_cal;

    fvector<3> x_org = upper(v);
    fvector<3> dx_cal, p_dx_cal, dx_d, p_dx_d, duv_d, p_duv_d;

    fvector<3> xd_tmp, p_xd_tmp;
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

    fvector<3> input;
    fvector<2> x_cal_pol;

    /* simulation for loop starts here */
    for(int n = 0;n < N;n++, t += dt){

        x = upper(v);
        dx = lower(v);

#ifndef PURE_PURSUIT
        out_v[n] = log(t, v, xd, dx_d, x_cal, dx_cal, e, tau, i);
#else
        out_v[n] = log(t, v, xd_tmp, dx_d, x_cal, dx_cal, e, tau, i);
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

        /* control input calculation */
#ifndef PURE_PURSUIT
        // check convergence and update xd if x_cal converged to xd
        if((xd - x_cal).norm() < eps){
            if(index_xd < n_xd){
                if(xd != xdc_arr[index_xd].xd) x_org = xd;
                
                type = xdc_arr[index_xd].type;
                xd_arc = xdc_arr[index_xd].arc;
                xd = xdc_arr[index_xd++].xd;

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
#endif

        p_dx_cal = dx_cal;
        dx_cal = 0.6 * dx_cal + 0.4 * (x_cal - p_x_cal)/ dt;

#ifdef POS_PD
        /* PD feedback in uv-coordinate */
        if(type == xd_t::Line || (xd - x_cal).norm() < 0.1){
            uv = Jux(arg_uv).trans()*(x_cal - xd);
            duv = Jux(arg_uv).trans()*(dx_cal - (xd - p_xd)/ dt);
        }
        else {
            float x_t_diff = x_cal_pol[1] - xd_arc[4];
            while(x_t_diff > pi)x_t_diff -= 2*pi;
            while(x_t_diff < -pi)x_t_diff += 2*pi;
            uv = fvector<3>{xd_arc[2]*(x_t_diff), arc_dir*(x_cal_pol[0] - xd_arc[2]), x_cal[2] - xd[2]};
            uv[1] += arc_dir *(Jux(arg_uv).trans()* dx_cal * dt)[0]/ xd_arc[2];
            printf("%6.2f \r\n", arg_uv);
            // (Jux(arg_uv).trans()*dx_cal).print();
            duv = Jux(arg_uv).trans()*(dx_cal - (xd - p_xd)/ dt);
            duv.print();
            duv[1] += arc_dir * 5 *((Jux(arg_uv).trans()*dx_cal)[0] * (Jux(arg_uv).trans()*dx_cal)[0] * dt)/ xd_arc[2];

            duv.print();
        }

        // if(fabsf(uv[0]) > 2) uv[0] = std::copysign(1, uv[0]);
        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Jux(arg_uv) * input;

#elif defined VEL_PD
        /* PD feedback in xy-velocity-dimention */
        uv = Jux(arg_uv).trans()*(x_cal - xd);
        // if(fabsf(uv[0]) > dx_d_max * dx_d_max /(2 * ddx)) uv[0] = std::copysign(dx_d_max * dx_d_max /(2 * ddx), uv[0]);
        p_duv_d = duv_d;
        duv_d = uv.each([](float _v){
            if(fabsf(_v)< ddx_buffer) return - std::copysign(sqrtf(fabsf(2 * ddx_slow * _v)), _v);
            else return - std::copysign(sqrtf(fabsf(2 * ddx * (fabsf(_v) - ddx_buffer) + 2 * ddx_slow * ddx_buffer)), _v);
        });

        duv_d = duv_d.each([](float _f){
            return fabsf(_f) > dx_d_max ? std::copysign(dx_d_max, _f) : _f;
        });
        if(duv_d.norm() > dx_d_max) duv_d *= dx_d_max / duv_d.norm();


        p_dx_d = dx_d;
        dx_d = Jux(arg_uv) * duv_d;

        fvector<3> diff_dx = dx_d - p_dx_d;
        float dx_max = diff_dx.abs().max();

        for(int i = 0;i < 3;i++){
            if(diff_dx[i] * dx_d[i] >= 0) diff_dx[i] = ddx * dt * diff_dx[i] / dx_max;

            if(fabsf(p_dx_d[i] - dx_d[i])> fabsf(diff_dx[i])) dx_d[i] = p_dx_d[i] + diff_dx[i];
        }

        input = ke * 14 * J(x_cal[2]).inv() * dx_d
                + Kp * J(x_cal[2]).inv()*(dx_d - dx_cal)
                + Kd * J(x_cal[2]).inv()*((dx_d - p_dx_d)-(dx_cal - p_dx_cal))/ dt;
        e = input;

#elif defined PURE_PURSUIT

        float dist = (xd - xd_tmp).norm();
        for(int k = index_xd;k < n_xd;k++){
            dist += (xd_arr[k] - xd_arr[k - 1]).norm();
        }
        vel = std::min(sqrtf(2 * ddx * fabsf(dist)), vel + ddx * dt);
        if(vel > dx_d_max) vel = dx_d_max;

        float remain;
        fvector<3> p_xd_remain = xd_tmp;

        p_xd_tmp = xd_tmp;
        xd_tmp += (xd - xd_tmp)/(xd - xd_tmp).norm()* vel * dt;
        while((p_xd_remain - xd)*(xd_tmp - xd)< 0){
            if(index_xd >= n_xd) break;

            remain = (xd_tmp - xd).norm();
            p_xd_remain = xd_tmp;
            xd_tmp = xd;

            if(xd != xd_arr[index_xd]) x_org = xd;

            xd = xd_arr[index_xd++];
            xd_tmp += (xd - xd_tmp)/(xd - xd_tmp).norm()* remain;

            if((x_org[0] != xd[0])||(x_org[1] != xd[1])){
                arg_uv = atan2(xd[1] - x_org[1], xd[0] - x_org[0]);
            }
            else arg_uv = 0;
        }
        dx_d = (xd_tmp - p_xd_tmp)/ dt;

        uv = Jux(arg_uv).trans()*(x_cal - xd_tmp);
        duv = Jux(arg_uv).trans()*(dx_cal - dx_d);

        // if(fabsf(uv[0]) > 2) uv[0] = std::copysign(1, uv[0]);
        input = - Kp * uv - Kd * duv;
        e = J(x_cal[2]).inv() * Jux(arg_uv) * input;

        // input = Kp *(xd_tmp - x_cal) + Kd *(dx_d - dx_cal);
        // e = J(x_cal[2]).inv() * input;


#endif

        float e_max = e.abs().max();
        if(e_max > max_voltage){
            for(int i = 0;i < 3;i++) e[i] *= max_voltage / e_max;
        }

        fvector<3> diff_e = e - pe;
        float de_max = diff_e.abs().max();

        for(int i = 0;i < 3;i++){
            diff_e[i] = de * diff_e[i] / de_max;
            if(pe[i] - e[i] > - diff_e[i])e[i] = pe[i] + diff_e[i];
            if(e[i] - pe[i] > diff_e[i])e[i] = pe[i] + diff_e[i];
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

    // i_ave /= cnt;
    // printf("i_ave: %6.2f \r\n", i_ave.mean());
    // i_ave.print();
    // printf("i_peak: \r\n");
    // i_peak.print();

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
