#ifndef COST_HPP
#define COST_HPP

// #include <Eigen/Cholesky>
// #include <iostream>
#include "./return_type_structs.hpp"

template <int T, int nx, int nu, typename type=float>
class Cost {
    public:
        virtual type cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const = 0;
        virtual void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Jacobians_Struct<T, nx, nu, type>& cjs) const = 0;
};

template <int T, int nx, int nu, typename type=float>
class QuadraticCost : public Cost<T, nx, nu, type> {
    public:
        QuadraticCost (const Eigen::Matrix<type, nx, nx>& Q_, const Eigen::Matrix<type, nu, nu>& R_, const Eigen::Matrix<type, nx, nx>& Qf_);

        type cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const override;

        void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Jacobians_Struct<T, nx, nu, type>& cjs) const override;

    private:
        const Eigen::Matrix<type, nx, nx>& Q;
        const Eigen::Matrix<type, nu, nu>& R;
        const Eigen::Matrix<type, nx, nx>& Qf;
};

#include "../source/cost.cpp"

#endif