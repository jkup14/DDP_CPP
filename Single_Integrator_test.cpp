#include <iostream>
#include <Eigen/Cholesky>
#include "integrator.cpp"
#include "solution_struct.cpp"
#include "cost.cpp"

//TODO eigen print formatting perhaps

int main() {
    const int T = 3;
    const typedef float type;
    const SingleIntegrator<T> Model = SingleIntegrator<T>();
    const int nx = Model.getNx();
    const int nu = Model.getNu();

    auto x0 = Eigen::Matrix<float, nx, 1>::Zero();
    Eigen::MatrixXd u0 = Eigen::MatrixXd::Zero(T, nu);
    // std::cout << x0 << std::endl;
    // std::cout << u0 << std::endl;

    // Test xdot
    std::cout << Model.xdot(x0, (Eigen::Matrix<float, nu, 1>()<<1,2).finished()) << std::endl;

    //Test differentiate dynamics and << override
    Dynamics_Derivatives_Struct<T, nx, nu> dds;
    Model.differentiate_dynamics(Eigen::Matrix<float, T, nx>::Random(), Eigen::Matrix<float, T-1, nu>::Random(), dds);
    std::cout << dds << std::endl;

    //Test Integrator
    EulerIntegrator<T, nx, nu> Dyn(Model, 0.1);
    std::cout << "xnext:" << endl << Dyn.propogate(x0, (Eigen::Matrix<float, nu, 1>()<<1,2).finished()) << std::endl;

    // //Test differentiate integrator and << override
    Integrator_Derivatives_Struct<T, nx, nu> ids;
    Dyn.differentiate_integrator(Eigen::Matrix<float, T, nx>::Random(), Eigen::Matrix<float, T-1, nu>::Random(), ids);
    std::cout << ids << std::endl;

    //Test cost
    auto Q = Eigen::Matrix<float, nx, nx>::Identity()*0;
    auto R = Eigen::Matrix<float, nu, nu>::Identity()*0.001;
    auto Qf = Eigen::Matrix<float, nx, nx>::Identity()*1;
    QuadraticCost<T, nx, nu> cost(Q, R, Qf);
    cout << "Cost:" << endl;
    cout << cost.cost(Eigen::Matrix<float, T, nx>::Constant(0), Eigen::Matrix<float, T-1, nu>::Constant(1), Eigen::Matrix<float, T, nx>::Constant(0)) << endl;
    Cost_Derivatives_Struct<T, nx, nu> cds;
    cost.differentiate_cost(Eigen::Matrix<float, T, nx>::Constant(0), Eigen::Matrix<float, T-1, nu>::Constant(1), Eigen::Matrix<float, T, nx>::Constant(0), cds);
    std::cout << cds << std::endl;

    // DDP:DDP ddp_solver(Dyn, )

    // DDP::Solution<T_, nx, nu> sol(
    //     Eigen::Matrix<float, T_, nx>::Zero(), 
    //     Eigen::Matrix<float, T_, nu>::Zero(), 
    //     std::vector<Eigen::Matrix<float, nu, nx> >(10, Eigen::Matrix<float,nu,nx>::Zero()),
    //     std::vector<float>(T, 1.5));
}