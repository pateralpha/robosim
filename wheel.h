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