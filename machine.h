#pragma once

#include <cstdlib>


constexpr float th_omega = 1e-2; // threshold value for determining stop or not

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

constexpr float inertia_of_cylinder(float _m, float _a){ return (_m * _a * _a)/ 2; }

constexpr float r = 0.05; // radius of the wheel
constexpr float R = 0.5; // radius of the machine / distance between CoG and wheels

constexpr float m_per_pulse = 0.2355e-3;

extern const wheel w1;
extern const wheel w2;
extern const wheel w3;

extern const fmatrix<3, 3> M;
extern const fmatrix<3, 3> I;

extern const fmatrix<3, 3> Juw;

inline fmatrix<3, 3> Rotr(float _theta){
    return fmatrix<3, 3>{
            cosf(_theta), - sinf(_theta), 0, 
            sinf(_theta), cosf(_theta), 0,
            0, 0, 1
    };
}

inline fmatrix<3, 3> dJux(float _theta, float _dtheta){
    return fmatrix<3, 3>{
            - sinf(_theta) * _dtheta, - cosf(_theta) * _dtheta, 0, 
            cosf(_theta) * _dtheta, - sinf(_theta) * _dtheta, 0,
            0, 0, 0
    };
}

inline fmatrix<3, 3> J(float _theta){ return Rotr(_theta) * Juw.inv(); }
inline fmatrix<3, 3> dJ(float _theta, float _dtheta){ return dJux(_theta, _dtheta) * Juw.inv(); }

inline fvector<6> f(fvector<6> _x, fvector<3> _e){
    fvector<3> x = top_v<6, 3>(_x), dx = bottom_v<6, 3>(_x);
    float theta = x[2], dtheta = dx[2];

    fvector<3> w = J(theta).inv() * dx;
    fvector<3> tau(w1.torque(w[0], _e[0]), w2.torque(w[1], _e[1]), w3.torque(w[2], _e[2]));

    fvector<3> ddx = (M + J(theta).trans().inv() * I * J(theta).inv()).inv() * J(theta).trans().inv() * (tau - I * dJ(theta, dtheta) * dx);

    return merge_v(dx, ddx);
}

inline fvector<3> du(fvector<3> _dx, fvector<3> _u){
    return Rotr(_u[2]).trans() * _dx;
}