#include "dynamics.cpp"
using std::ostream;
using std::cout;
using std::endl;

template <int T, int nx, int nu, typename type=float>
class Integrator {
    public: 
        Integrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            :dynamics(dynamics_), dt(dt_) {}

        virtual Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const = 0;
        virtual void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Jacobians_Struct<T, nx, nu>& int_derivs) const = 0;

    protected:
        type dt;
        const Dynamics<T, nx, nu, type>& dynamics; 
};

template <int T, int nx, int nu, typename type=float>
class EulerIntegrator : public Integrator<T, nx, nu, type> {
    public:
        EulerIntegrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            : Integrator<T, nx, nu, type>(dynamics_, dt_) {}
        
        Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const {
            return x + this->dt * this->dynamics.xdot(x,u);
        }

        void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) const {
            Dynamics_Jacobians_Struct<T, nx, nu, type> dds;
            this->dynamics.differentiate_dynamics(X, U, dds);
            for (int t = 0; t<T-1; t++) {
                ijs.fx[t] = Eigen::Matrix<type,nx,nx>::Identity() + this->dt * dds.xdotx[t];
                ijs.fu[t] = this->dt * dds.xdotu[t];
            }
        }

};
