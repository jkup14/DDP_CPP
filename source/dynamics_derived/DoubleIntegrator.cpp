#ifndef DI_CPP
#define DI_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <vector>
#include "../../include/dynamics.hpp"
#include "../../include/dynamics_derived/DoubleIntegrator.hpp"

using std::ostream;
using std::endl;

template <int T, typename type>
DoubleIntegrator<T, type>::DoubleIntegrator():Dynamics_Abstract<T, 4, 2, type>(std::vector<std::string>{"x","y","dx/dt","dy/dt"}) {}

template <int T, typename type>
Eigen::Matrix<type, 4, 1> DoubleIntegrator<T, type>::xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
            return (Eigen::Matrix<type, 4, 1>()<< x(Eigen::seq(2,3)), u).finished();
            return Eigen::Matrix<type, 4, 1>::Zero();
        }

template <int T, typename type>
void DoubleIntegrator<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 4, 2>& djs) const {
    Eigen::Matrix<type, 4, 4> xdotx_t = Eigen::Matrix<type, 4, 4>::Zero();
    xdotx_t(0,2) = 1;
    xdotx_t(1,3) = 1;
    djs.xdotx = std::vector<Eigen::Matrix<type, 4, 4> >(T-1, xdotx_t);
    Eigen::Matrix<type, 4, 2> xdotu_t = Eigen::Matrix<type, 4, 2>::Zero();
    xdotu_t(2,0) = 1;
    xdotu_t(3,1) = 1;
    djs.xdotu = std::vector<Eigen::Matrix<type, 4, 2> >(T-1, xdotu_t);
}

#endif
