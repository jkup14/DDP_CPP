#ifndef SI_HPP
#define SI_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../return_type_structs.hpp"
#include "../dynamics.hpp"
using std::ostream;
using std::endl;

template <int T, typename type=float>
class SingleIntegrator : public Dynamics_Abstract<T, 2, 2, type> {
    public:  
        SingleIntegrator();
        Eigen::Matrix<type, 2, 1> xdot(const Eigen::Matrix<type, 2, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const override;
};

#include "../../source/dynamics_derived/SingleIntegrator.cpp"

#endif
