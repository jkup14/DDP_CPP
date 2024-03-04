#include <iostream>
#include "include/ddp.hpp"
#include "include/dynamics_derived/SingleIntegrator.hpp"

//TODO eigen print formatting perhaps

int main() {
    const int T = 10;
    float dt = 0.1;
    const typedef float type;
    const SingleIntegrator<T> Model = SingleIntegrator<T>();
    const int nx = Model.getNx();
    const int nu = Model.getNu();

    Eigen::Matrix<float,nx,1> x0 = Eigen::Matrix<float, nx, 1>::Zero();
    Eigen::Matrix<float,T,nx> X = Eigen::Matrix<float, T, nx>::Zero();
    Eigen::Matrix<float,1,nx> x_goal = (Eigen::Matrix<float,1,nx>()<<1,3).finished();
    Eigen::Matrix<float,T,nx> X_track = x_goal.replicate<T,1>();
    Eigen::Matrix<float,T-1,nu> U = Eigen::Matrix<float, T-1, nu>::Constant(1);

    // Test xdot
    std::cout << Model.xdot(x0, (Eigen::Matrix<float, nu, 1>()<<1,2).finished()) << std::endl;

    // //Test differentiate dynamics and << override
    Dynamics_Jacobians_Struct<T, nx, nu> djs;
    Model.differentiate_dynamics(X, U, djs);
    std::cout << djs << std::endl;

    //Test Integrator
    EulerIntegrator<T, nx, nu> Dyn(Model, dt);
    std::cout << "xnext:" << endl << Dyn.propogate(x0, (Eigen::Matrix<float, nu, 1>()<<1,2).finished()) << std::endl;

    // Test differentiate integrator and << override
    Integrator_Jacobians_Struct<T, nx, nu> ijs;
    Dyn.differentiate_integrator(X, U, ijs);
    std::cout << ijs << std::endl;

    //Test Rollout
    Dyn.rollout_controls(x0, U, X);
    cout << "X:" << endl << X << endl;

    //Test cost
    auto Q = Eigen::Matrix<float, nx, nx>::Identity()*0;
    auto R = Eigen::Matrix<float, nu, nu>::Identity()*0.1;
    auto Qf = Eigen::Matrix<float, nx, nx>::Identity()*1;
    QuadraticCost<T, nx, nu> cost(Q, R, Qf);
    cout << "Q:" << endl << cost.Q << endl 
         << "R:" << endl << cost.R << endl 
         << "Qf:" << endl << cost.Qf << endl;
    cout << "Cost:" << endl;
    cout << cost.cost(X, U, X_track) << endl;
    Cost_Jacobians_Struct<T, nx, nu> cjs;
    cost.differentiate_cost(X, U, X_track, cjs);
    std::cout << cjs << std::endl;

    //Test DDP
    Value_Jacobians_Struct<nx, nu> vjs(cjs.Lx.row(T-1), cjs.Lxx.back());
    Feedback_Struct<T, nx, nu> fs;
    cout << vjs << endl << fs << endl;

    //Initialize Solver
    DDP::DDP_Solver<T, nx, nu> solver(Dyn, cost);
    
    std::pair<float, float> deltaVs = {0,0};
    solver.backward_pass(cjs, ijs, fs, 0, deltaVs);
    cout << fs << endl;
    cout << "deltaV:" << endl << deltaVs.first << "," << deltaVs.second <<endl;

    Eigen::Matrix<float,T,nx> Xnew = Eigen::Matrix<float, T, nx>::Zero();
    Eigen::Matrix<float,T-1,nu> Unew = Eigen::Matrix<float, T-1, nu>::Zero();
    solver.rollout_ddp_gains(x0, X, U, Xnew, Unew, fs, 1);
    cout << "Xnew:" << endl << Xnew << endl;
    cout << "Unew:" << endl << Unew << endl;

    DDP::Solution<T, nx, nu> sol(
        Eigen::Matrix<float, T, nx>::Zero(), 
        Eigen::Matrix<float, T-1, nu>::Zero(), 
        std::vector<Eigen::Matrix<float, nu, nx> >(10-1, Eigen::Matrix<float,nu,nx>::Zero()),
        std::vector<float>(T, 1.5),
        10, 10);
    cout << sol << endl;

}