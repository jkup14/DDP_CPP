#include "dynamics.cpp"
using std::ostream;
using std::cout;
using std::endl;

template <int T, int nx, int nu, typename type=float>
struct Integrator_Derivatives_Struct {
    Integrator_Derivatives_Struct() {
        fx = std::vector<Eigen::Matrix<type, nx, nx> >(T-1, Eigen::Matrix<type, nx, nx>::Zero());
        fu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
    }
    std::vector<Eigen::Matrix<type, nx, nx> > fx;
    std::vector<Eigen::Matrix<type, nx, nu> > fu;

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Integrator_Derivatives_Struct<T_, nx_, nu_, type_>& ids);
};

template <int T, int nx, int nu, typename type=float>
ostream& operator<<(ostream& os, const Integrator_Derivatives_Struct<T, nx, nu, type>& ids) {
    os << "fx:" << endl;
    // os << ids.fx[0] << endl;
    for (int t = 0; t<T-1; t++) {os << ids.fx.at(t) << endl;}
    os << "fu:" << endl;
    for (int t = 0; t<T-2; t++) {os << ids.fu.at(t) << endl;}
    os << ids.fu.at(T-2);
    return os;
}


template <int T, int nx, int nu, typename type=float>
class Integrator {
    public: 
        Integrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            :dynamics(dynamics_), dt(dt_) {}

        virtual Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) = 0;
        virtual void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Derivatives_Struct<T, nx, nu>& int_derivs) const = 0;

    protected:
        type dt;
        const Dynamics<T, nx, nu, type>& dynamics; 
};

template <int T, int nx, int nu, typename type=float>
class EulerIntegrator : Integrator<T, nx, nu, type> {
    public:
        EulerIntegrator(const Dynamics<T, nx, nu, type>& dynamics_, type dt_)
            : Integrator<T, nx, nu, type>(dynamics_, dt_) {}
        
        Eigen::Matrix<type, nx, 1> propogate(const Eigen::Matrix<type, nx, 1>& x, const Eigen::Matrix<type, nu, 1>& u) {
            return x + this->dt * this->dynamics.xdot(x,u);
        }

        void differentiate_integrator(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, Integrator_Derivatives_Struct<T, nx, nu, type>& ids) const {
            Dynamics_Derivatives_Struct<T, nx, nu, type> dds;
            this->dynamics.differentiate_dynamics(X, U, dds);
            for (int t = 0; t<T-1; t++) {
                ids.fx[t] = Eigen::Matrix<type,nx,nx>::Identity() + this->dt * dds.xdotx[t];
                ids.fu[t] = this->dt * dds.xdotu[t];
            }
        }

};

// class Autodifferentiator {
//     Autodifferentiator (Integrator* integrator_): integrator(integrator_) {}
//     Eigen::MatrixXd jacobian (const Integrator* integrator, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double t) {
//             VectorXreal xv(x);
//             VectorXreal uv(u);
//             real tv(t);
//             VectorXreal F;
//             Eigen::MatrixXd Fx = jacobian(integrator->propogate, wrt(xv), at(xv, uv, tv), F);
//             Eigen::MatrixXd Fu = jacobian(integrator->propogate, wrt(uv), at(xv, uv, tv), F);
//         }
//     private:
//         Integrator* integrator;
// }