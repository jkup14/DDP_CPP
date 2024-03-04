#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../include/return_type_structs.hpp"
using std::ostream;
using std::endl;


template <int T, int nx, int nu, typename type=float>
class Dynamics_Abstract {
    protected:
        Dynamics_Abstract(std::vector<std::string> state_names_);
        std::vector<std::string> state_names;
    public:
        virtual Eigen::Matrix<type, nx, 1> xdot(const Eigen::Matrix<float, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const = 0;
        virtual void differentiate_dynamics(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Dynamics_Jacobians_Struct<T, nx, nu>& dyn_derivs) const = 0;
        constexpr int getNx() const;
        constexpr int getNu() const;
        std::vector<std::string> get_state_names();
};

#include "../source/dynamics.cpp"

#endif
