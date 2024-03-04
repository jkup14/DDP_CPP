#ifndef COST_CPP
#define COST_CPP

#include "../include/cost.hpp"
// #include <iostream>
using std::cout;
using std::endl;


// TODO: maybe store X_diff for performance isntead of calculating twice

template <int T, int nx, int nu, typename type>
QuadraticCost<T, nx, nu, type>::QuadraticCost (const Eigen::Matrix<type, nx, nx> Q_, const Eigen::Matrix<type, nu, nu> R_, const Eigen::Matrix<type, nx, nx> Qf_):
            Q(Q_), R(R_), Qf(Qf_) {}
            // {Q = Q_; R = R_; Qf = Qf_;}

template <int T, int nx, int nu, typename type>
type QuadraticCost<T, nx, nu, type>::cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const {
    type sum = 0;
    // cout << this->Q << this->R << this->Qf << endl;
    const Eigen::Matrix<type, T, nx> X_diff = X-X_track;
    // cout << "Xdiff" << endl << X_diff << endl;
    for (int t = 0; t<T-1; t++) {
        auto x_cost = X_diff.row(t) * Q * X_diff.row(t).transpose();
        auto u_cost = U.row(t) * R * U.row(t).transpose();
        sum += x_cost(0) + u_cost(0);
        // cout << U.row(t) * R * U.row(t).transpose() << endl;
    }
    sum += X_diff.row(T-1) * Qf * X_diff.row(T-1).transpose();
    return sum;
}

template <int T, int nx, int nu, typename type>
void QuadraticCost<T, nx, nu, type>::differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Jacobians_Struct<T, nx, nu, type>& cjs) const {
    const Eigen::Matrix<type, T, nx> X_diff = X-X_track;
    for (int t = 0; t<T-1; t++) {
        cjs.Lx.row(t) = X_diff.row(t) * Q;
        cjs.Lu.row(t) = U.row(t) * R;
        cjs.Lxx.at(t) = Q;
        cjs.Luu.at(t) = R;
        cjs.Lxu.at(t) = Eigen::Matrix<type, nx, nu>::Zero();
    }
    cjs.Lx.row(T-1) = X_diff.row(T-1) * Qf;
    cjs.Lxx.at(T-1) = Qf;
}

#endif