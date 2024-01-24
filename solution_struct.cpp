// #include <iostream>
#include <Eigen/Cholesky>


namespace DDP {
    template <int T, int nx, int nu>
    struct Solution {
        Eigen::Matrix<float, T, nx> X;
        Eigen::Matrix<float, T, nu> U;
        std::vector<Eigen::Matrix<float, nu, nx> > K_u;
        std::vector<float> J;
        Solution(Eigen::Matrix<float, T, nx> X_, 
                 Eigen::Matrix<float, T, nu> U_,
                 std::vector<Eigen::Matrix<float, nu, nx> > K_u_, 
                 std::vector<float> J_)
                 : X(X_), U(U_), K_u(K_u_), J(J_) {}
    };
}

