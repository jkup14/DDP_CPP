#ifndef RETURN_TYPE_STRUCTS_H
#define RETURN_TYPE_STRUCTS_H

#include <iostream>
#include <Eigen/Cholesky>
#include <vector>


namespace DDP {
    template <int T, int nx, int nu, typename type=float>
    struct Solution {
        Eigen::Matrix<type, T, nx> X;
        Eigen::Matrix<type, T-1, nu> U;
        std::vector<Eigen::Matrix<type, nu, nx> > K;
        std::vector<type> J;
        int it;
        long long ms;

        Solution(Eigen::Matrix<type, T, nx> X_, 
                 Eigen::Matrix<type, T-1, nu> U_,
                 std::vector<Eigen::Matrix<type, nu, nx> > K_, 
                 std::vector<type> J_,
                 int it_,
                 long long ms_);

    template <int T_, int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Solution<T_, nx_, nu_, type_>& vjs);
    };

    template <int T, int nx, int nu, typename type>
    std::ostream& operator<<(std::ostream& os, const Solution<T, nx, nu, type>& sol) ;
}

template <int T, int nx, int nu, typename type=float>
struct Dynamics_Jacobians_Struct {

    Dynamics_Jacobians_Struct() ;

    std::vector<Eigen::Matrix<type, nx, nx> > xdotx;
    std::vector<Eigen::Matrix<type, nx, nu> > xdotu;

    template <int T_, int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Dynamics_Jacobians_Struct<T_, nx_, nu_, type_>& djs);
};

template <int T, int nx, int nu, typename type=float>
std::ostream& operator<<(std::ostream& os, const Dynamics_Jacobians_Struct<T, nx, nu, type>& djs) ;

template <int T, int nx, int nu, typename type=float>
struct Integrator_Jacobians_Struct {

    Integrator_Jacobians_Struct() ;

    std::vector<Eigen::Matrix<type, nx, nx> > fx;
    std::vector<Eigen::Matrix<type, nx, nu> > fu;

    template <int T_, int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Integrator_Jacobians_Struct<T_, nx_, nu_, type_>& ijs);
};

template <int T, int nx, int nu, typename type=float>
std::ostream& operator<<(std::ostream& os, const Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) ;

template <int T, int nx, int nu, typename type=float>
struct Cost_Jacobians_Struct {

    Cost_Jacobians_Struct () ;

    Eigen::Matrix<type, T, nx> Lx;
    Eigen::Matrix<type, T-1, nu> Lu;
    std::vector<Eigen::Matrix<type, nx, nx> > Lxx;
    std::vector<Eigen::Matrix<type, nu, nu> > Luu;
    std::vector<Eigen::Matrix<type, nx, nu> > Lxu;

    template <int T_, int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Cost_Jacobians_Struct<T_, nx_, nu_, type_>& cjs);
};

template <int T, int nx, int nu, typename type=float>
std::ostream& operator<<(std::ostream& os, const Cost_Jacobians_Struct<T, nx, nu, type>& cjs) ;

template<int T, int nx, int nu, typename type=float>
struct Feedback_Struct {

    Feedback_Struct () ;

    Eigen::Matrix<type, T-1, nu> k;
    std::vector<Eigen::Matrix<type, nu, nx> > K;

    template <int T_, int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Feedback_Struct<T_, nx_, nu_, type_>& fs);
};

template <int T, int nx, int nu, typename type>
std::ostream& operator<<(std::ostream& os, const Feedback_Struct<T, nx, nu, type>& fs) ;


template<int nx, int nu, typename type=float>
struct Value_Jacobians_Struct {

    Value_Jacobians_Struct(Eigen::Matrix<type, 1, nx> Lx_T, Eigen::Matrix<type, nx, nx> Lxx_T) ;

    Eigen::Matrix<type, nx, 1> Vx;
    Eigen::Matrix<type, nx, nx> Vxx;

    template <int nx_, int nu_, typename type_>
    friend std::ostream& operator<<(std::ostream& os, const Value_Jacobians_Struct<nx_, nu_, type_>& vjs);
};

template <int nx, int nu, typename type>
std::ostream& operator<<(std::ostream& os, const Value_Jacobians_Struct<nx, nu, type>& vjs) ;


#include "../source/return_type_structs.cpp"


#endif