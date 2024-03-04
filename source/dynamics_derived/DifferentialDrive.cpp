#ifndef DD_CPP
#define DD_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <vector>
#include "../../include/dynamics.hpp"
#include "../../include/dynamics_derived/DifferentialDrive.hpp"

using std::ostream;
using std::endl;


template <int T, typename type>
DifferentialDrive<T, type>::DifferentialDrive():wheel_radius(0.2),width(0.2),Dynamics_Abstract<T, 3, 2, type>(std::vector<std::string>{"x","y","θ"}) {}

template <int T, typename type>
Eigen::Matrix<type, 3, 1> DifferentialDrive<T, type>::xdot(const Eigen::Matrix<type, 3, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
    return (Eigen::Matrix<type, 3, 1>() << 
    cos(x(2)) * (u(0)+u(1)),
    sin(x(2)) * (u(0)+u(1)),
    (u(0)-u(1)) / width
    ).finished()*wheel_radius/2;
}

template <int T, typename type>
void DifferentialDrive<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 3>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 3, 2>& djs) const {
    // djs is already initialized to zeros when constructed
    for (int t = 0; t<T-1; t++) {
        type sin_angle = sin(X(t,2));
        type cos_angle = cos(X(t,2));

        djs.xdotx.at(t) = (((Eigen::Matrix<type, 3, 3>()<< 
        0,0,-sin_angle,
        0,0,cos_angle,
        0,0,0
        ).finished())
        *wheel_radius*(U(t,0)+U(t,1))/2);

        djs.xdotu.at(t) = (((Eigen::Matrix<type, 3, 2>()<< 
        cos_angle, cos_angle,
        sin_angle, sin_angle,
        1/width, -1/width).finished())
        *wheel_radius/2);
    }
}


#endif
