#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "dynamics.hpp"
using std::ostream;
using std::cout;
using std::endl;

template <int T, int nx, int nu, typename type=float>
class Integrator {
    public: 
        Integrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_);
        virtual Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const = 0;
        virtual void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Jacobians_Struct<T, nx, nu>& ijs) const = 0;
        void rollout_controls(const Eigen::Matrix<float, nx, 1>& x0, const Eigen::Matrix<float, T-1, nu>& U, Eigen::Matrix<float, T, nx>& X) const ;

    protected:
        type dt;
        const Dynamics<T, nx, nu, type>& dynamics; 
};

template <int T, int nx, int nu, typename type=float>
class EulerIntegrator : public Integrator<T, nx, nu, type> {
    public:
        EulerIntegrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_);
        
        Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const override;

        void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) const override;
};

#include "../source/integrator.cpp"

#endif