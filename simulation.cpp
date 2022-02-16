#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>

constexpr float th_omega = 1e-2; // threshold value for determining stop or not
constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }
constexpr float pi = M_PI;

template <typename T>
inline int sign(T _v){ return _v >= 0 ? 1 : -1; }

struct vector3{
    float e[3];

    vector3(void){ memset(e, 0, sizeof(float)* 3); }
    vector3(float _v[]){ memcpy(e, _v, sizeof(float)* 3); }
    vector3(float _v1, float _v2, float _v3){
        e[0] = _v1;
        e[1] = _v2;
        e[2] = _v3;
    }
    float &operator [](int i){ return e[i]; }

    vector3 add(vector3 _v){
        vector3 ret(this->e);
        for(int k = 0;k < 3;k++) ret.e[k] += _v.e[k];
        return ret;
    }

    vector3 sub(vector3 _v){
        vector3 ret(this->e);
        for(int k = 0;k < 3;k++) ret.e[k] -= _v.e[k];
        return ret;
    }

    vector3 operator +(vector3 _v){ return add(_v); }
    vector3 operator -(vector3 _v){ return sub(_v); }

    vector3 mul(float _f){
        vector3 ret(this->e);
        for(int k = 0;k < 3;k++) ret.e[k] *= _f;
        return ret;
    }

    vector3 operator *(float _f){ return mul(_f); }

    vector3 abs(void){ return vector3(fabsf(e[0]), fabsf(e[1]), fabsf(e[2])); }
    float max(void){ return std::max(std::max(e[0], e[1]), e[2]); }

    void print(void){
        printf("(%6.2f; %6.2f; %6.2f)\r\n", e[0], e[1], e[2]);
    }
};

vector3 operator *(float _f, vector3 _v){ return _v.mul(_f); }

struct vector6{
    float e[6];

    vector6(void){ memset(e, 0, sizeof(float)* 6); }
    vector6(float _v1, float _v2, float _v3, float _v4, float _v5, float _v6){
        e[0] = _v1;
        e[1] = _v2;
        e[2] = _v3;
        e[3] = _v4;
        e[4] = _v5;
        e[5] = _v6;
    }
    vector6(float _v[]){ memcpy(e, _v, sizeof(float)* 6); }
    vector6(vector3 _v1, vector3 _v2){
        memcpy(e, _v1.e, sizeof(float)* 3);
        memcpy(e + 3, _v2.e, sizeof(float)* 3);
    }
    float &operator [](int i){ return e[i]; }

    vector3 upper(void){
        return vector3(e[0], e[1], e[2]);
    }

    vector3 lower(void){
        return vector3(e[3], e[4], e[5]);
    }

    vector6 add(vector6 _v){
        vector6 ret(this->e);
        for(int k = 0;k < 6;k++) ret.e[k] += _v.e[k];
        return ret;
    }

    vector6 sub(vector6 _v){
        vector6 ret(this->e);
        for(int k = 0;k < 6;k++) ret.e[k] -= _v.e[k];
        return ret;
    }

    vector6 operator +(vector6 _v){ return add(_v); }
    vector6 operator -(vector6 _v){ return sub(_v); }

    vector6 mul(float _f){
        vector6 ret(this->e);
        for(int k = 0;k < 6;k++) ret.e[k] *= _f;
        return ret;
    }

    vector6 operator *(float _f){ return mul(_f); }

    void print(void){
        printf("(%6.2f; %6.2f; %6.2f; %6.2f; %6.2f; %6.2f)\r\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    }
};

vector6 operator *(float _f, vector6 _v){ return _v.mul(_f); }

struct matrix3{
    float e[9];

    matrix3(void){}
    matrix3(vector3 _v1, vector3 _v2, vector3 _v3){
        memcpy(e, _v1.e, sizeof(float) * 3);
        memcpy(e + 3, _v2.e, sizeof(float) * 3);
        memcpy(e + 6, _v3.e, sizeof(float) * 3);
    }
    matrix3(float _e[]){ memcpy(e, _e, sizeof(float)* 9); }
    matrix3(vector3 _v){ 
        for(int k = 0;k < 9;k++) e[k] = 0;
        e[0] = _v.e[0];
        e[4] = _v.e[1];
        e[8] = _v.e[2];
    }

    float &operator ()(int i){ return e[i]; }

    matrix3 trans(void){
        matrix3 ret(*this);
        ret.e[1] = e[3];
        ret.e[2] = e[6];
        ret.e[3] = e[1];
        ret.e[5] = e[7];
        ret.e[6] = e[2];
        ret.e[7] = e[5];

        return ret;
    }

    float det(void){
        return e[0] * e[4] * e[8] + e[3] * e[7] * e[2] + e[6] * e[1] * e[5]
                - e[0] * e[7] * e[5] - e[3] * e[1] * e[8] - e[6] * e[4] * e[2];
    }

    matrix3 inv(void){
        matrix3 ret;
        float d = det();
        ret.e[0] = e[4] * e[8] - e[7] * e[5];
        ret.e[1] = - e[1] * e[8] + e[7] * e[2];
        ret.e[2] = e[1] * e[5] - e[4] * e[2];
        ret.e[3] = - e[3] * e[8] + e[6] * e[5];
        ret.e[4] = e[0] * e[8] - e[6] * e[2];
        ret.e[5] = - e[0] * e[5] + e[3] * e[2];
        ret.e[6] = e[3] * e[7] - e[6] * e[4];
        ret.e[7] = - e[0] * e[7] + e[6] * e[1];
        ret.e[8] = e[0] * e[4] - e[3] * e[1];
        for(int k = 0;k < 9;k++) ret.e[k] /= d;
        return ret;
    }

    matrix3 add(matrix3 _v){
        matrix3 ret(this->e);
        for(int k = 0;k < 9;k++) ret.e[k] += _v.e[k];
        return ret;
    }

    vector3 prod(vector3 _v){
        vector3 ret;
        ret.e[0] = e[0] * _v.e[0] + e[3] * _v.e[1] + e[6] * _v.e[2];
        ret.e[1] = e[1] * _v.e[0] + e[4] * _v.e[1] + e[7] * _v.e[2];
        ret.e[2] = e[2] * _v.e[0] + e[5] * _v.e[1] + e[8] * _v.e[2];
        return ret;
    }

    matrix3 prod(matrix3 _v){
        vector3 v1(_v.e[0], _v.e[1], _v.e[2]),
                v2(_v.e[3], _v.e[4], _v.e[5]),
                v3(_v.e[6], _v.e[7], _v.e[8]);
        matrix3 ret(prod(v1), prod(v2), prod(v3));
        return ret;
    }

    matrix3 operator +(matrix3 _v){ return add(_v); }
    vector3 operator *(vector3 _v){ return prod(_v); }
    matrix3 operator *(matrix3 _v){ return prod(_v); }

    void print(void){
        printf("(%6.2f, %6.2f, %6.2f;\r\n", e[0], e[3], e[6]);
        printf(" %6.2f, %6.2f, %6.2f;\r\n", e[1], e[4], e[7]);
        printf(" %6.2f, %6.2f, %6.2f)\r\n", e[2], e[5], e[8]);
    }
};

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

wheel w1(0.0083 * 14, 0.12, 0.33, 0.8);

matrix3 M(vector3(15, 15, inertia_of_cylinder(15, R)));
matrix3 I(vector3(inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r), inertia_of_cylinder(0.1, r)));

float Juw_elm[] = {
        cosf(2 * pi / 3)/ r, cosf(4 * pi / 3)/ r, cosf(0)/ r,
        sinf(2 * pi / 3)/ r, sinf(4 * pi / 3)/ r, sinf(0)/ r,
        R / r, R / r, R / r
};
matrix3 Juw(Juw_elm);

matrix3 Jux(float _theta){
    float elm[] = {
            cosf(_theta), sinf(_theta), 0, 
            - sinf(_theta), cosf(_theta), 0,
            0, 0, 1
    };
    return matrix3(elm);
}

matrix3 dJux(float _theta, float _dtheta){
    float elm[] = {
            - sinf(_theta) * _dtheta, cosf(_theta) * _dtheta, 0, 
            - cosf(_theta) * _dtheta, - sinf(_theta) * _dtheta, 0,
            0, 0, 0
    };
    return matrix3(elm);
}

matrix3 J(float _theta){
    return Jux(_theta) * Juw.inv();
}

matrix3 dJ(float _theta, float _dtheta){
    return dJux(_theta, _dtheta) * Juw.inv();
}

vector6 f(vector6 _x, vector3 _e){
    vector3 x = _x.upper(), dx = _x.lower();
    float theta = x[2], dtheta = dx[2];

    vector3 omega = J(theta).inv() * dx;

    vector3 tau(
            w1.torque(omega[0], _e[0]), 
            w1.torque(omega[1], _e[1]), 
            w1.torque(omega[2], _e[2])
    );

    vector3 ddx = (M + J(theta).trans().inv() * I * J(theta).inv()).inv() * J(theta).trans().inv() * (tau - I * dJ(theta, dtheta) * dx);

    return vector6(dx, ddx);
}

vector3 du(vector3 _dx, vector3 _u){
    return Jux(_u[2]).inv() * _dx;
}

void state(vector6 &_x, vector3 &_tau, vector3 &_i, vector3 _e){
    vector3 x = _x.upper(), dx = _x.lower();
    float theta = x[2], dtheta = dx[2];

    vector3 omega = J(theta).inv() * dx;

    _tau = vector3(
            w1.torque(omega[0], _e[0]), 
            w1.torque(omega[1], _e[1]), 
            w1.torque(omega[2], _e[2])
    );
    _i = vector3(
            w1.current(omega[0], _e[0]), 
            w1.current(omega[1], _e[1]),
            w1.current(omega[2], _e[2])
    );
}

constexpr int N = 300;

struct output{
    float t;

    vector3 x;
    vector3 dx;

    vector3 e;
    vector3 tau;
    // vector3 i;
} out_v[N];

output log(float _t, vector6 _x, vector3 _e, vector3 _tau, vector3 _i){
    output ret;
    ret.t = _t;
    ret.x = _x.upper();
    ret.dx = _x.lower();
    ret.e = _e;
    ret.tau = _tau;
    // ret.i = _i;
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
    fname[0] = 'a';

    vector3 machine_pos[6];
    matrix3 Rr = Jux(2*pi/3), Rl = Jux(-2*pi/3);

    machine_pos[4] = {-0.08, -0.45, 0};
    machine_pos[5] = {0.08, -0.45, 0};

    machine_pos[0] = Rr * machine_pos[4];
    machine_pos[1] = Rr * machine_pos[5];

    machine_pos[2] = Rl * machine_pos[4];
    machine_pos[3] = Rl * machine_pos[5];


    for(int n = 0;n < N;n++) out_v[n] = {0, vector3(), vector3(), vector3(), vector3()};
    constexpr float dt = 0.01;
    float t = 0;

    vector6 v;
    vector3 x, dx;
    vector3 e, pe;

    vector3 u;

    vector3 omega, tau, i;

    constexpr float max_voltage = 18; // muximum voltage
    constexpr float sat_time = 0.2; // minimum time from 0 V to maximum
    constexpr float de = 18 / sat_time * dt;

    constexpr float m_per_pulse = 0.2355e-3;

    matrix3 Kp(vector3(1, 1, 2));
    matrix3 Kd(vector3(0.1, 0.1, 0.2));

    vector3 xd(2, 1, pi/4);

    int enc_u = 0, enc_v = 0;
    int p_enc_u = enc_u, p_enc_v = enc_v;
    vector3 del_u;
    vector3 x_cal, p_x_cal;

    for(int n = 0;n < N;n++, t += dt){
        // printf("%5.2f  ", t);
        // v.print();
        // printf("       ");
        // e.print();

        x = v.upper();
        dx = v.lower();

        out_v[n] = {t, x, dx, e, tau};

        // controller

        // dead reckoning
        enc_u = u[0] / m_per_pulse;
        enc_v = u[1] / m_per_pulse;

        del_u = vector3(
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
        e = J(x_cal[2]).inv()*(Kp *(xd - x_cal) - Kd * (x_cal - p_x_cal)*(1 / dt));

        for(int i = 0;i < 3;i++){
            if(fabsf(e[i]) > 18) e[i] = std::copysign(18, e[i]);
        }

        vector3 diff_e = e - pe;
        float de_max = diff_e.abs().max();

        for(int i = 0;i < 3;i++){
            diff_e[i] = de * diff_e[i] / de_max;
            if(pe[i] - e[i] > - diff_e[i])e[i] = pe[i] + diff_e[i];
            if(e[i] - pe[i] > diff_e[i])e[i] = pe[i] + diff_e[i];
        }
        pe = e;

        // end

        v = rk4<vector6>(v, [=](vector6 _v){ return f(_v, e); }, dt);
        u = rk4<vector3>(u, [=](vector3 _u){ return du(dx, _u); }, dt);

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
    vector3 machine[6];

    constexpr int skip = 5;

    for(int n = 0;n < N;n += skip){
        fprintf(fp, "%5.2f %6.3f %6.3f %6.1f %5.1f %5.1f %5.1f ",
                out_v[n].t, out_v[n].x[0], out_v[n].x[1], out_v[n].x[2] / pi * 180, out_v[n].e[0], out_v[n].e[1], out_v[n].e[2]);

        for(int k = 0;k < 6;k++){
            machine[k] = out_v[n].x + Jux(out_v[n].x[2]) * machine_pos[k];
            fprintf(fp, "%6.3f %6.3f ", machine[k][0], machine[k][1]);
        }
        fprintf(fp, "\r\n");
    }

    fclose(fp);
    printf("file closed.\r\n");

    return 0;
}
