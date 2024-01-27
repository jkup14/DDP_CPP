#include <iostream>
#include <Eigen/Cholesky>
#include "ddp.cpp"

//TODO eigen print formatting perhaps

int main() {
    const int T = 100;
    float dt = 0.01;
    DifferentialDrive<T> Model = DifferentialDrive<T>();
    const int nx = Model.getNx();
    const int nu = Model.getNu();

    Eigen::Matrix<float,nx,1> x0 = Eigen::Matrix<float, nx, 1>::Zero();
    Eigen::Matrix<float,T,nx> X = Eigen::Matrix<float, T, nx>::Zero();
    Eigen::Matrix<float,1,nx> x_goal = (Eigen::Matrix<float,1,nx>()<<1,2,0).finished();
    Eigen::Matrix<float,T,nx> X_track = x_goal.replicate<T,1>();
    Eigen::Matrix<float,T-1,nu> U = Eigen::Matrix<float, T-1, nu>::Constant(0);

    //Discretize Dynamics
    EulerIntegrator<T, nx, nu> dynamics(Model, dt);

    //Cost
    auto Q = Eigen::Matrix<float, nx, nx>::Identity()*0;
    auto R = Eigen::Matrix<float, nu, nu>::Identity()*0.00001;
    auto Qf = Eigen::Matrix<float, nx, nx>::Identity()*100;
    QuadraticCost<T, nx, nu> cost(Q, R, Qf);
    
    //Initialize Solver
    DDP::solver_args args;
    args.verbose = 1;
    args.conv_threshold = 1e-5;
    DDP::DDP_Solver<T, nx, nu> solver(dynamics, cost, args);
    // Solve!
    cout << "Solving..." << endl;
    DDP::Solution<T, nx, nu> sol = solver.solve(x0, U, X_track);
    cout << "Final error: " << sol.X.row(T-1)-x_goal << endl;
}