#include <Eigen/Cholesky>
#include <iostream>
using std::ostream;
using std::endl;

template <int T, int nx, int nu, typename type=float>
class Dynamics {
    public:
        virtual Eigen::Matrix<type, nx, 1> xdot(const Eigen::Matrix<float, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const = 0;
        virtual void differentiate_dynamics(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Dynamics_Jacobians_Struct<T, nx, nu>& dyn_derivs) const = 0;
        constexpr int getNx() const {return nx;}
        constexpr int getNu() const {return nu;}
};

template <int T, typename type=float>
class SingleIntegrator : public Dynamics<T, 2, 2, type> {
    public:  
        // SingleIntegrator() {} 
        Eigen::Matrix<type, 2, 1> xdot(const Eigen::Matrix<type, 2, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
            return Eigen::Matrix<type, 2, 1>(u);
        }

        void differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& dyn_derivs) const {
            dyn_derivs.xdotx = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Zero());
            dyn_derivs.xdotu = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Identity());
        }
    
};

// template <typename type, int T>
// class DoubleIntegrator : public Dynamics {
//     public:  
//         DoubleIntegrator(): Dynamics<type, T, 4, 2>() {} 
//         Eigen::Matrix<type, Dynamics.getNx(), 1> xdot(const Eigen::Matrix<float, Dynamics.getNx(), 1>& x, const Eigen::Matrix<float, Dynamics.getNu(), 1>& u) const{
//             return (Eigen::VectorXd(4)<<x(Eigen::seq(2,4)),u).finished();
//         }

//         void differentiate_dynamics(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Dynamics_Jacobians_Struct<T, nx, nu>& dyn_derivs) const {
//             //
//         }
// };



