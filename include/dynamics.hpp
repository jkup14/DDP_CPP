#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include "../include/return_type_structs.hpp"
using std::ostream;
using std::endl;


template <int T, int nx, int nu, typename type=float>
class Dynamics {
    public:
        virtual Eigen::Matrix<type, nx, 1> xdot(const Eigen::Matrix<float, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const = 0;
        virtual void differentiate_dynamics(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Dynamics_Jacobians_Struct<T, nx, nu>& dyn_derivs) const = 0;
        constexpr int getNx() const;
        constexpr int getNu() const;
};

template <int T, typename type=float>
class SingleIntegrator : public Dynamics<T, 2, 2, type> {
    public:  
        Eigen::Matrix<type, 2, 1> xdot(const Eigen::Matrix<type, 2, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const override;
};

template <int T, typename type=float>
class DoubleIntegrator : public Dynamics<T, 4, 2, type> {
    public:  
        Eigen::Matrix<type, 4, 1> xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;

        void differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 4, 2>& djs) const override;
};

template <int T, typename type=float>
class DifferentialDrive : public Dynamics<T, 3, 2, type> {
    public:  
        DifferentialDrive();
        Eigen::Matrix<type, 3, 1> xdot(const Eigen::Matrix<type, 3, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const override;
        void differentiate_dynamics(const Eigen::Matrix<type, T, 3>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 3, 2>& djs) const override;
    private:
        const type wheel_radius;
        const type width;
    
};
// template <int T, typename type=float>
// class CartPole : public Dynamics<T, 4, 1, type> {
//     public:  
//         CartPole() {
//             mass_pole = 0.05;
//             mass_cart = 1;
//             gravity = 9.8;
//             pole_length = 2;
//         } 
//         Eigen::Matrix<type, 4, 1> xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 1, 1>& u) const{
//             // state is [displacement, angle, velocity, angular_velocity]
//             type sin_angle = sin(x(1));
//             type sin_angle_sqr = pow(sin_angle, 2);
//             type cos_angle = cos(x(1));
//             type ang_velocity_sqr = pow(x(3),2);
//             return (Eigen::Matrix<type, 4, 1>() << 
//                 x(2),
//                 x(3),
//                 (mass_pole*sin_angle
//                     * (pole_length * ang_velocity_sqr * gravity * cos_angle)
//                     + u) 
//                     / (mass_cart + mass_pole * sin_angle_sqr),
//                 (-mass_pole * pole_length * ang_velocity_sqr * cos_angle * sin_angle 
//                     - (mass_cart + mass_pole) * gravity * sin_angle
//                     - cos_angle * u)
//                     / (pole_length * (mass_cart + mass_pole * sin_angle_sqr))
//                 ).finished();
//         }

//         void differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const {
//             djs.xdotx = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Zero());
//             djs.xdotu = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Identity());
//         }

//     private:
//         type mass_pole;
//         type mass_cart;
//         type gravity;
//         type pole_length;
    
// };

#include "../source/dynamics.cpp"

#endif
