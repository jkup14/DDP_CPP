#ifndef CP_HPP
#define CP_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../return_type_structs.hpp"
#include "../dynamics.hpp"
using std::ostream;
using std::endl;

template <int T, typename type=float>
class CartPole : public Dynamics_Abstract<T, 4, 1, type> {
    public:  
        CartPole();
        Eigen::Matrix<type, 4, 1> xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 1, 1>& u) const;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 1>& U, Dynamics_Jacobians_Struct<T, 4, 1>& djs) const;
    private:
        const type mass_pole;
        const type mass_cart;
        const type gravity;
        const type pole_length;
};

#include "../../source/dynamics_derived/CartPole.cpp"

#endif
