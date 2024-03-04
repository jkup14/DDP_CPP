#ifndef SI_CPP
#define SI_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <vector>
#include "../../include/dynamics.hpp"
#include "../../include/dynamics_derived/SingleIntegrator.hpp"
using std::ostream;
using std::endl;

template <int T, typename type>
SingleIntegrator<T, type>::SingleIntegrator():Dynamics_Abstract<T, 2, 2, type>(std::vector<std::string>{"x","y"}) {}

template <int T, typename type>
Eigen::Matrix<type, 2, 1> SingleIntegrator<T, type>::xdot(const Eigen::Matrix<type, 2, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
    return Eigen::Matrix<type, 2, 1>(u);
}

template <int T, typename type>
void SingleIntegrator<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const {
    djs.xdotx = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Zero());
    djs.xdotu = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Identity());
}

#endif
