#ifndef DI_HPP
#define DI_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../return_type_structs.hpp"
#include "../dynamics.hpp"
using std::ostream;
using std::endl;

template <int T, typename type=float>
class DoubleIntegrator : public Dynamics_Abstract<T, 4, 2, type> {
    public:  
        DoubleIntegrator();
        Eigen::Matrix<type, 4, 1> xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 4, 2>& djs) const override;
};

#include "../../source/dynamics_derived/DoubleIntegrator.cpp"

#endif
