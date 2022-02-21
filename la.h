#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <array>
#include <limits>


constexpr float pi = M_PI;

using std::array;

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
    const fvector<N> operator -(void) const{ return *this * (- 1); }

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

    const float mean(void) const{
        float sum = 0;
        for(float v : *this) sum += v;
        return sum / N;
    }

    const fvector<N> each(std::function<float(float)> f) const{
        fvector<N> v;
        for_range(i, N) v[i] = f((*this)[i]);
        return v;
    }

};

template <std::size_t N>
const fvector<N> operator *(const float _a, const fvector<N> &_v);


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

    const fmatrix<N, M> operator -(void) const{ return *this * (- 1); }

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


template <std::size_t N>
const fvector<N> operator *(const float _a, const fvector<N> &_v){
    fvector<N> v = _v;
    return v *= _a;
}

template <std::size_t N, std::size_t M>
const fmatrix<N, M> operator *(const float _a, const fmatrix<N, M> &_m){
    fmatrix<N, M> m = _m;
    return m *= _a;
}

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
