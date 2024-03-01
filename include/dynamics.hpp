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
template <int T, typename type=float>
class CartPole : public Dynamics<T, 4, 1, type> {
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

#include "../source/dynamics.cpp"

#endif
