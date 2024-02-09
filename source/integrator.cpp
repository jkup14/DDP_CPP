#ifndef INTEGRATOR_CPP
#define INTEGRATOR_CPP

using std::ostream;
using std::cout;
using std::endl;

template <int T, int nx, int nu, typename type>
Integrator<T, nx, nu, type>::Integrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            :dynamics(dynamics_), dt(dt_) {}

template <int T, int nx, int nu, typename type>
void Integrator<T, nx, nu, type>::rollout_controls(const Eigen::Matrix<float, nx, 1>& x0, const Eigen::Matrix<float, T-1, nu>& U, Eigen::Matrix<float, T, nx>& X) const {
        X.row(0) = x0.transpose();
        for (int t = 0; t<T-1; t++) {
            X.row(t+1) = propogate(X.row(t), U.row(t)).transpose();
        }
}

template <int T, int nx, int nu, typename type>
EulerIntegrator<T, nx, nu, type>::EulerIntegrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            : Integrator<T, nx, nu, type>(dynamics_, dt_) {}

template <int T, int nx, int nu, typename type>
Eigen::Matrix<type, nx, 1> EulerIntegrator<T, nx, nu, type>::propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) const {
    return x + this->dt * this->dynamics.xdot(x,u);
}

template <int T, int nx, int nu, typename type>
void EulerIntegrator<T, nx, nu, type>::differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) const {
    Dynamics_Jacobians_Struct<T, nx, nu, type> dds;
    this->dynamics.differentiate_dynamics(X, U, dds);
    for (int t = 0; t<T-1; t++) {
        ijs.fx[t] = Eigen::Matrix<type,nx,nx>::Identity() + this->dt * dds.xdotx[t];
        ijs.fu[t] = this->dt * dds.xdotu[t];
    }
}

#endif
