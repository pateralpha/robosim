#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>

using std::array;

constexpr float th_omega = 1e-2; // threshold value for determining stop or not
constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }
constexpr float pi = M_PI;

template <typename T>
inline int sign(T _v){ return _v >= 0 ? 1 : -1; }
template <typename T>
inline void swap(T &_v1, T &_v2){
    T t = _v1;
    _v1 = _v2;
    _v2 = t;
}

#define for_range(i, n) for(int i = 0;i < n;i++)

template <std::size_t N>
struct fvector : public array<float, N>{

    fvector(void) { for(float &v : *this) v = 0; }
    fvector(const float _f) { for(float &v : *this) v = _f; }
    template <class... A>
    fvector(const float _f, A... _rem){
        init(0, _f, _rem...);
    }

    template <class... A>
    void init(int _i, const float _f, A... _rem){
        if(_i >= N) return;
        (*this)[_i] = _f;
        init(_i + 1, _rem...);
    }

    void init(int _i, const float _f){
        if(_i >= N) return;
        (*this)[_i] = _f;
        for(int i = _i + 1;i < N;i++) (*this)[i] = 0;
    }

    void print(void) const{
        printf("(");
        for(float f : *this) printf("%6.2f; ", f);
        printf("\b\b)\r\n");
    }

    fvector<N> &add_eq(const fvector<N> &_v){
        for_range(i, N) (*this)[i] += _v[i];
        return *this;
    }
    fvector<N> &operator +=(const fvector<N> &_v){ return add_eq(_v); }
    const fvector<N> operator +(const fvector<N> &_v) const{
        fvector<N> v = *this;
        return v += _v;
    }

    fvector<N> &sub_eq(const fvector<N> &_v){
        for_range(i, N) (*this)[i] -= _v[i];
        return *this;
    }
    fvector<N> &operator -=(const fvector<N> &_v){ return sub_eq(_v); }
    const fvector<N> operator -(const fvector<N> &_v) const{
        fvector<N> v = *this;
        return v -= _v;
    }

    fvector<N> &mul_eq(const float _a){
        for(float &v : *this) v *= _a;
        return *this;
    }
    fvector<N> &operator *=(const float _a){ return mul_eq(_a); }
    const fvector<N> operator *(const float _a) const{
        fvector<N> v = *this;
        return v *= _a;
    }

    const float inner_prod(const fvector<N> &_v) const{
        float prod = 0;
        for_range(i, N) prod += (*this)[i] * _v[i];
        return prod;
    }
    const float operator *(const fvector<N> &_v) const{ return inner_prod(_v); }

    fvector<N> &div_eq(const float _a){
        for(float &v : *this) v /= _a;
        return *this;
    }
    fvector<N> &operator /=(const float _a){ return div_eq(_a); }
    const fvector<N> operator /(const float _a) const{
        fvector<N> v = *this;
        return v /= _a;
    }

    const fvector<N> abs(void) const{
        fvector<N> v;
        for_range(i, N) v[i] = fabsf((*this)[i]);
        return v;
    }

    const float max(void) const{
        float max = - std::numeric_limits<float>::infinity();
        for(float v : *this) max = std::max(max, v);
        return max;
    }

    const float min(void) const{
        float min = std::numeric_limits<float>::infinity();
        for(float v : *this) min = std::min(min, v);
        return min;
    }

    const float norm(void) const{
        float quad = 0;
        for(float v : *this) quad += v*v;
        return sqrtf(quad);
    }


};

template <std::size_t N>
const fvector<N> operator *(const float _a, const fvector<N> &_v){
    fvector<N> v = _v;
    return v *= _a;
}


template <std::size_t N, std::size_t M>
struct fmatrix : public array<fvector<M>, N>{

    fmatrix(void) { for(auto &v : *this) v = fvector<M>(0); }
    fmatrix(const int _i) { for(auto &v : *this) v = fvector<M>(_i); }
    fmatrix(const float _f) { for(auto &v : *this) v = fvector<M>(_f); }
    fmatrix(const float _v[]){ for_range(i, N) for_range(j, M) (*this)[i][j] = _v[i * M + j]; }

    template <class... A>
    fmatrix(const float _f, A... _rem){
        init(0, _f, _rem...);
    }

    template <class... A>
    void init(int _i, const float _f, A... _rem){
        if(_i >= N * M) return;
        (*this)[_i / M][_i % M] = _f;
        init(_i + 1, _rem...);
    }

    void init(int _i, const float _f){
        if(_i >= N * M) return;
        (*this)[_i / M][_i % M] = _f;
        for(int i = _i + 1;i < N;i++) (*this)[i] = 0;
    }

    void print(void) const{
        printf("(");
        for_range(i, N){
            for_range(j, M) printf("%6.2f; ", (*this)[i][j]);
            if(i < N - 1) printf("\r\n ");
        }
        printf("\b\b)\r\n");
    }

    const fmatrix<M, N> trans(void) const{
        fmatrix<M, N> mat;
        for_range(i, N){
            for_range(j, M) mat[j][i] = (*this)[i][j];
        }
        return mat;
    }

    const fvector<M> row_vector(int _i) const{
        fvector<M> v((*this)[_i]);
        return v;
    }

    const fvector<N> column_vector(int _j) const{
        return trans().row_vector(_j);
    }

    fmatrix<N, M> &add_eq(const fmatrix<N, M> &_m){
        for_range(i, N){
            for_range(j, M) (*this)[i][j] += _m[i][j];
        }
        return *this;
    }
    fmatrix<N, M> &operator +=(const fmatrix<N, M> &_m){ return add_eq(_m); }
    const fmatrix<N, M> operator +(const fmatrix<N, M> &_m) const{
        fmatrix<N, M> m = *this;
        return m += _m;
    }

    fmatrix<N, M> &sub_eq(const fmatrix<N, M> &_m){
        for_range(i, N){
            for_range(j, M) (*this)[i][j] -= _m[i][j];
        }
        return *this;
    }
    fmatrix<N, M> &operator -=(const fmatrix<N, M> &_m){ return sub_eq(_m); }
    const fmatrix<N, M> operator -(const fmatrix<N, M> &_m) const{
        fmatrix<N, M> m = *this;
        return m -= _m;
    }

    fmatrix<N, M> &mul_eq(const float _a){
        for(auto &v : *this) v *= _a;
        return *this;
    }
    fmatrix<N, M> &operator *=(const float _a){ return mul_eq(_a); }
    const fmatrix<N, M> operator *(const float _a) const{
        fmatrix<N, M> m = *this;
        return m *= _a;
    }

    fmatrix<N, M> &div_eq(const float _a){
        for(auto &v : *this) v /= _a;
        return *this;
    }
    fmatrix<N, M> &operator /=(const float _a){ return div_eq(_a); }
    const fmatrix<N, M> operator /(const float _a) const{
        fmatrix<N, M> m = *this;
        return m /= _a;
    }

    template <std::size_t L>
    const fmatrix<N, L> prod(const fmatrix<M, L> &_m) const{
        fmatrix<N, L> m;
        for_range(i, N){
            for_range(j, L) m[i][j] = row_vector(i) * _m.column_vector(j);
        }
        return m;
    }
    template <std::size_t L>
    const fmatrix<N, L> operator *(const fmatrix<M, L> &_m) const{ return prod(_m); }

    const fvector<N> prod(const fvector<M> &_m) const{
        fvector<N> v;
        for_range(i, N) v[i] = row_vector(i) * _m;
        return v;
    }
    const fvector<N> operator *(const fvector<M> &_v) const{ return prod(_v); }

    const float minor_det(int _i, int _j) const{
        fmatrix<N - 1, N - 1> m;

        int pi = 0, pj = 0;
        for_range(i, N){
            if(i == _i) continue;
            for_range(j, M){
                if(j == _j) continue;
                m[pi][pj++] = (*this)[i][j];
            }
            pj = 0;
            pi++;
        }
        return m.det();
    }

    const float sgn(int _i, int _j) const{ return (_i + _j)% 2 ? - 1 : 1; }

    const float det(void) const{
        if(N != M) return 0;

        fvector<N> v;
        int p;
        auto m = lu(v, p);

        float d = 1;
        for_range(i, N) d *= m[i][i];
        return d *(p % 2 ? - 1 : 1);
    }

    const fmatrix<N, N> adj(void) const{
        fmatrix<N, N> m;
        for_range(i, N){
            for_range(j, N) m[j][i] = sgn(i, j) * minor_det(i, j);
        }
        return m;
    }

    const fmatrix<N, N> inv(void) const{
        fmatrix<N, N> m;
        for_range(i, N){
            fvector<N> v(0);
            v[i] = 1;
            m[i] = solve(v);
        }
        return m.trans();
    }

    const fmatrix<N, N> lu(fvector<N> &_p, int &_pivot) const{
        fmatrix<N, N> m = *this;
        float a_max, a;
        int ip;

        _pivot = 0;
        for_range(i, N){
            a_max = 0;
            ip = 0;
            for(int k = i;k < N;k++){
                if(fabsf(m[k][i]) > a_max){
                    a_max = fabsf(m[k][i]);
                    ip = k;
                }
            }
            if(!a_max) return fmatrix<N, N>(0);

            _p[i] = ip;
            if(i != ip){
                _pivot++;
                for(int j = i;j < N;j++) swap(m[i][j], m[ip][j]);
            }

            for(int j = i + 1;j < N;j++){
                a = - m[j][i] / m[i][i];
                m[j][i] = a;
                for(int k = i + 1;k < N;k++) m[j][k] += a * m[i][k];
            }
        }
        return m;
    }

    const fvector<N> solve(const fvector<N> &_v) const{
        fvector<N> v;
        int p;
        auto m = lu(v, p);

        return solve(_v, m, v);
    }

    const fvector<N> solve(fvector<N> _v, const fmatrix<N, N> &_m, const fvector<N> &_p) const{

        for_range(i, N){
            swap(_v[i], _v[_p[i]]);
            for(int j = i + 1;j < N;j++) _v[j] += _m[j][i] * _v[i];
        }

        float t;
        for(int i = N - 1;i >= 0;i--){
            t = 0;
            for(int j = i + 1;j < N;j++) t += _m[i][j] * _v[j];
            _v[i] = (_v[i] - t)/ _m[i][i];
        }

        return _v;
    }
};

template <>
struct fmatrix<1, 1> : public array<fvector<1>, 1>{

    fmatrix(void){ (*this)[0][0] = 0; }
    fmatrix(const float _f){ (*this)[0][0] = _f; }
    fmatrix(const float _v[]){ (*this)[0][0] = _v[0]; }
    void print(void) const{
        printf("(%6.2f)\r\n", (*this)[0][0]);
    }
    const float minor_dat(int _i, int _j) const{ return 0; }
    const float det(void) const{ return (*this)[0][0]; }
};

template <std::size_t N, std::size_t M>
const fmatrix<N, M> zeros(void) { return fmatrix<N, M>(0); }

template <std::size_t N>
const fmatrix<N, N> eye(void) {
    fmatrix<N, N> m;
    for_range(i, N) m[i][i] = 1;
    return m;
}
template <std::size_t N>
const fmatrix<N, N> diag(const float _f[]){
    fmatrix<N, N> m;
    for_range(i, N) m[i][i] = _f[i];
    return m;
}
template <std::size_t N>
const fmatrix<N, N> diag(const fvector<N> _v){
    fmatrix<N, N> m;
    for_range(i, N) m[i][i] = _v[i];
    return m;
}
template <std::size_t N, class... A>
const fmatrix<N, N> diag(const float _f, A... _rem){
    fmatrix<N, N> m;
    diag(m, 0, _f, _rem...);
    return m;
}
template <std::size_t N, class... A>
void diag(fmatrix<N, N> &_m, int _i, const float _f, A... _rem){
    if(_i >= N) return;
    _m[_i][_i] = _f;
    diag(_m, _i + 1, _rem...);
}
template <std::size_t N>
void diag(fmatrix<N, N> &_m, int _i, const float _f){
    if(_i >= N) return;
    _m[_i][_i] = _f;
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
    constexpr float eps = 1;


    constexpr int n_xd = 3;
    fvector<3> xd_arr[n_xd] = {{2, 1, pi/4}, {2, 3, pi/4}, {1, 4, 0}};
    // fvector<3> xd_arr[n_xd] = {{2, 1, pi/4}};
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
