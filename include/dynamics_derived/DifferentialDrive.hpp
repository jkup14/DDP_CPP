#ifndef DD_HPP
#define DD_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../return_type_structs.hpp"
#include "../dynamics.hpp"
using std::ostream;
using std::endl;

template <int T, typename type=float>
class DifferentialDrive : public Dynamics_Abstract<T, 3, 2, type> {
    public:  
        DifferentialDrive();
        Eigen::Matrix<type, 3, 1> xdot(const Eigen::Matrix<type, 3, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 3>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 3, 2>& djs) const override;
    private:
        const type wheel_radius;
        const type width;
    
};

#include "../../source/dynamics_derived/DifferentialDrive.cpp"

#endif
