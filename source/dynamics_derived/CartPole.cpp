#ifndef CP_CPP
#define CP_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <vector>
#include "../../include/dynamics.hpp"
#include "../../include/dynamics_derived/CartPole.hpp"
using std::ostream;
using std::endl;

template <int T, typename type>
CartPole<T,type>::CartPole():mass_pole(0.05), mass_cart(1), gravity(9.8), pole_length(2), Dynamics_Abstract<T, 4, 1, type>(std::vector<std::string>{"x","θ","dx/dt","dθ/dt"}){}

template <int T, typename type>
Eigen::Matrix<type, 4, 1> CartPole<T,type>::xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 1, 1>& u_) const{
    // state is [displacement, angle, velocity, angular_velocity]
    type sin_angle = sin(x(1));
    type sin_angle_sqr = pow(sin_angle, 2);
    type cos_angle = cos(x(1));
    type ang_velocity_sqr = pow(x(3),2);
    type u = u_(0);
    return (Eigen::Matrix<type, 4, 1>() << 
        x(2),
        x(3),
        (mass_pole*sin_angle
            * (pole_length * ang_velocity_sqr * gravity * cos_angle)
            + u) 
            / (mass_cart + mass_pole * sin_angle_sqr),
        (-mass_pole * pole_length * ang_velocity_sqr * cos_angle * sin_angle 
            - (mass_cart + mass_pole) * gravity * sin_angle
            - cos_angle * u)
            / (pole_length * (mass_cart + mass_pole * sin_angle_sqr))
        ).finished();
}

template <int T, typename type>
void CartPole<T,type>::differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 1>& U, Dynamics_Jacobians_Struct<T, 4, 1>& djs) const {

    for (int t = 0; t<T-1; t++) {
        type sin_angle = sin(X(t,1));
        type sin_angle_sqr = pow(sin_angle, 2);
        type cos_angle = cos(X(t,1));
        type ang_velocity = X(t,3);
        type ang_velocity_sqr = pow(ang_velocity,2);
        auto u = U(t);
        djs.xdotx.at(t) = (Eigen::Matrix<type, 4, 4>()<< 
        0,0,1,0,
        0,0,0,1,

        0,
        (mass_pole*cos_angle*(pole_length*ang_velocity_sqr+gravity*cos_angle)-gravity*mass_pole*sin_angle_sqr)
        /(mass_pole*sin_angle_sqr+mass_cart)
        -(2*mass_pole*cos_angle*sin_angle*(u+mass_pole*sin_angle*(pole_length*ang_velocity_sqr+gravity*cos_angle)))
        /pow(mass_pole*sin_angle_sqr+mass_cart,2),
        0,
        (2*pole_length*mass_pole*ang_velocity*sin_angle)
        /(mass_pole*sin_angle_sqr+mass_cart),

        0,
        (pole_length*mass_pole*ang_velocity_sqr*(-pow(cos_angle,2)+sin_angle_sqr)
        -gravity*(mass_cart+mass_pole)*cos_angle
        +u*sin_angle)
        /(pole_length*(mass_pole*sin_angle_sqr+mass_cart))
        +(2*mass_pole*cos_angle*sin_angle*(pole_length*mass_pole*cos_angle*sin_angle*ang_velocity_sqr
        +u*cos_angle+gravity*sin_angle*(mass_cart+mass_pole)))
        /(pole_length*pow((mass_pole*sin_angle_sqr+mass_cart),2)),
        0,
        -(2*mass_pole*ang_velocity*cos_angle*sin_angle)
        /(mass_pole*sin_angle_sqr+mass_cart)
        ).finished();

        djs.xdotu.at(t) = ((Eigen::Matrix<type, 4, 1>() <<
        0,
        0,
        1,
        -cos_angle/pole_length
        ).finished())/(mass_pole*sin_angle_sqr+mass_cart);

    }
}

#endif
