#include <vector>
#include "./solution_struct.cpp"
#include "integrator.cpp"
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include "cost.cpp"

template<typename type, int T, int nx, int nu>
struct Feedback_Struct {
    Feedback_Struct () {
        k = Eigen::Matrix<type, T, nu>::Zero();
        K = std::vector<Eigen::Matrix<type, nu, nx> >(T-1, Eigen::Matrix<type, nu, nx>::Zero());
    }
    Eigen::Matrix<type, T, nu> k;
    std::vector<Eigen::Matrix<type, nu, nx> > K;
};

template<typename type, int nx, int nu>
struct Value_Derivatives_Struct {
    Value_Derivatives_Struct() {
        Vx = Eigen::Matrix<type, nx, 1>::Zero();
        Vxx = Eigen::Matrix<type, nx, nx>::Zero();
    }
    Eigen::Matrix<type, nx, 1> Vx;
    Eigen::Matrix<type, nx, nx> Vxx;
};

template <typename type, int T, int nx, int nu>
class DDP_solver {
    public:
        DDP_solver(const Integrator<nx, nu>& dynamics_,
            const Cost<T, nx, nu>& cost_,
            const int max_iters_ = 500,
            const double conv_threshold_ = 1e-3,
            const double reg_mu_min_ = 10e-6,
            const double reg_mu_max_ = 10e10,
            const double reg_delta_0_ = 2,
            const double ls_alpha_0 = 1,
            const double ls_alpha_f = 10e-3,
            const int ls_alpha_len = 11,
            const bool toggle_ls = true,
            const bool toggle_reg = true)
          : dynamics(&dynamics_), cost(cost_),
            max_iters(max_iters_), conv_threshold(conv_threshold_),
            reg_mu_min(reg_mu_min_), reg_mu_max(reg_mu_max_), reg_delta_0(reg_delta_0_) {

            this->differentiateLQR(this->dynamics, this->cost);
            
            // Handling of linesearch and regularization options
            if (!(toggle_ls || toggle_reg)){
                throw std::runtime_error("Must toggle Linesearch and/or Regularization!");
            } 
            if (toggle_ls) {
                if (!(ls_alpha_0 > ls_alpha_f)) {
                    throw std::runtime_error("Alpha0 must be bigger than Alphaf");
                }
                // Building learning rate vector
                double log_a0 = log(ls_alpha_0);
                double log_af = log(ls_alpha_f);
                double step = (log_a0-log_af)/ls_alpha_len;
                this->alphas.reserve(ls_alpha_len);
                for (int i = ls_alpha_len; i--; i >= 0) {
                    this->alphas.push_back(pow(10, log_af+i*step));
                }
            } else {
                // If linesearch is off the learning rate is always 1
                this->alphas.push_back(1);
            }
            
            }

        DDP::Solution solve(const Eigen::Matrix<type, nx, 1>& x0,
              Eigen::Matrix<type, T, nu>& U,
              Eigen::Matrix<type, T, nx>& X_track, double mu = 10e-6, double delta = 0) {

            if (delta == 0) {
                delta = this->reg_delta_0;
            }
            int max_iters = this->max_iters;
            double conv_threshold = this->conv_threshold;
            int N = ubar.size();
            // auto dynamics = this->dyn_integrator.propogate;
            auto cost = this->cost;
            int n = xbar[0].size();
            int m = ubar[0].size();
            std::vector<double> J(max_iters);

            type L = cost(xbar, ubar, xdesired);
            double dV = L;
            int it = 0;
            

            Eigen<type, T, nu> k_u;
            std::vector<Eigen::Matrix<type, nu, nx> > K_u(T);
            std::pair<type, type> deltaV;

            Cost_Derivatives_Struct<type, T, nx, nu> cost_derivative_struct();

            auto Vx = Eigen::Matrix<type, nx, 1>::Zero();
            auto Vxx = Eigen::Matrix<type, nx, nx>::Zero();

            Eigen::Matrix<type, T, nx> X;
            Eigen::Matrix<type, T, nu> U;
            while (it < max_iters && dV > conv_threshold) {
                J[it] = L;

                cost.cost_derivatives(X, U, X_track, cost_derivative_struct)

                

                bool forwardpass = false;
                double L_new = 0;
                x[0] = x0;
                for (int jj = 0; jj++; jj < this->alphas.size()) {
                    double alpha = alphas[jj];
                    Eigen::VectorXd dx = x[0] - xbar[0];
                    for (int k = 0; k++; k<N-1) {
                        Eigen::VectorXd du = alpha*k_u[k] + K_u[k]*dx;
                        u[k] = ubar[k] + du;
                        x[k+1] = this->dyn_integrator->propogate(x[k], u[k]);
                        dx = x[k+1] - xbar[k+1];
                    }
                    L_new = cost(x, u, xdesired);
                    double true_reduction = L-L_new;
                    double expected_reduction = -alpha*(deltaV[0]+alpha*deltaV[1]);
                    double z= true_reduction/expected_reduction;
                    if (expected_reduction < 0) {
                        z = (true_reduction > 0) - (true_reduction < 0);
                        std::cout << "Negative expected reduction in cost!" << std::endl;
                    } if (z > 0) {
                        forwardpass = true;
                        break;
                    }
                }
                if (forwardpass) {
                    this->decreaseReg(mu, delta);
                    xbar = x;
                    ubar = u;
                    dV = L - L_new;
                    L = L_new;
                } else {
                    this->increaseReg(mu, delta);
                    std::cout << "LS failed, increasing mu to " << mu << std::endl;
                    if (mu > reg_mu_max) {
                        x = xbar;
                        u = ubar;
                        std::cout << "Exiting, max regularization reached!" << std::endl;
                        break;
                    }
                }
                it++;
            }

        std::cout << "DDP done in " << it << " iterations." << std::endl;
        DDP::Solution sol = DDP::Solution(x, u, K_u, J);
        return sol;
        }

        void forward_propogate_controls(const Eigen::Matrix<float, T, nx>& x0, const Eigen::Matrix<float, T, nu>& U, const Eigen::Matrix<float, T, nx>& X) {
            X.row(0) = x0;
            for (int t = 0; t<T; t++) {
                X.row(t+1) = Integrator.propogate(X.row(t), U.row(t));
            }
        }

        void backward_pass(const Cost_Derivatives_Struct& cds, const Dynamics_Derivates_Struct& dds, Value_Derivates_Struct& vds, Feedback_Struct& bps, type mu, std::pair<type, type>& deltaV){
            for (int t = T-2; t--; t >= 0) {
                auto fxTVxx = dds.fx.transpose() * vds.Vxx;

                auto Qx = cds.Lx + dds.fx.transpose() * vds.Vx;
                auto Qu = cds.Lu + dds.fu.transpose() * vds.Vx;
                auto Qxx = cds.Lxx + fxTVxx * dds.fx;
                auto Qxu = cds.Lxu + fxTVxx * dds.fu;
                auto Qux = Qxu.transpose();
                auto Quu = cds.Luu + dds.fu.transpose() * vds.Vxx * dds.fu + mu * Eigen::Matrix<type, nu, nu>::Identity();

                auto invQuu = Quu.inverse();
                Backward_Pass_Struct.K.at(t) = -invQuu * Qux;
                Backward_Pass_Struct.k.row(t) = -invQuu * Qu;

                vds.Vx = Qx + Qxu * k.row(t);
                vds.Vxx = Qxx + Qxu * K.at(t);
                vds.Vxx = 0.5 * (Vxx + Vxx.transpose());

                deltaV.first += k.row(t).transpose() * Qu;
                deltaV.second += 0.5 * k.row(t).transpose() * Quu * k.row(t);
                }
        }

        // void forward_pass(const Eigen::Matrix<float, T, nx>& X, const Eigen::Matrix<float, T, nu>& U, )
    private:
        const Integrator<nx, nu>& dynamics;
        const Cost<T, nx, nu>& cost
        const int max_iters;
        const double conv_threshold;
        const double reg_mu_min;
        const double reg_mu_max;
        const double reg_delta_0;
        std::vector<double> alphas; 

        void increaseReg(double &mu, double &delta) {
            delta = std::max(reg_delta_0, delta*reg_delta_0);
            mu = std::max(reg_mu_min, mu*delta);
        }

        void decreaseReg(double &mu, double &delta) {
            delta = std::min(1/reg_delta_0, delta/reg_delta_0);
            double md = mu*delta;
            if (md >= reg_mu_min) {
                mu = md;
            } else {
                mu = 0;
            }
        }



 


};
