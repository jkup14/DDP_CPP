#include <vector>
#include "return_type_structs.cpp"
#include "integrator.cpp"
#include "cost.cpp"
#include <Eigen/LU>
#include "sign_func.cpp"




namespace DDP {
    template <typename type=float>
    struct solver_args {
        int max_iters = 500;
        type conv_threshold = 1e-3;
        type reg_mu_min = 10e-6;
        type reg_mu_max = 10e10;
        type reg_delta_0 = 2;
        type ls_alpha_0 = 1;
        type ls_alpha_f = 10e-3;
        int ls_alpha_len = 11;
        bool toggle_ls = true;
        bool toggle_reg = true;
        int verbose = 0;
    };

    template <int T, int nx, int nu, typename type=float>
    class DDP_Solver {
        public:
            DDP_Solver(const Integrator<T, nx, nu, type>& dynamics_,
                const Cost<T, nx, nu, type>& cost_,
                const solver_args<type> args = solver_args<type>()) 
                : dynamics(dynamics_), 
                  cost(cost_),
                  max_iters(args.max_iters), 
                  conv_threshold(args.conv_threshold),
                  reg_mu_min(args.reg_mu_min),
                  reg_mu_max(args.reg_mu_max),
                  reg_delta_0(args.reg_delta_0),
                  verbose(args.verbose) {

                
                // Handling of linesearch and regularization options
                if (!(args.toggle_ls || args.toggle_reg)){
                    throw std::runtime_error("Must toggle Linesearch and/or Regularization!");
                } 
                if (args.toggle_ls) {
                    if (!(args.ls_alpha_0 > args.ls_alpha_f)) {
                        throw std::runtime_error("Alpha0 must be bigger than Alphaf");
                    }
                    // Building learning rate vector
                    type log_a0 = log(args.ls_alpha_0);
                    type log_af = log(args.ls_alpha_f);
                    type step = (log_a0-log_af)/args.ls_alpha_len;
                    alphas.reserve(args.ls_alpha_len);
                    for (int i = args.ls_alpha_len; i >= 0; i--) {
                        alphas.push_back(pow(10, log_af+i*step));
                    }
                } else {
                    // If linesearch is off the learning rate is always 1
                    alphas.push_back(1);
                }
                
                }

            DDP::Solution<T, nx, nu, type> solve(const Eigen::Matrix<type, nx, 1>& x0,
                  Eigen::Matrix<type, T-1, nu>& U,
                  Eigen::Matrix<type, T, nx>& X_track, double mu = 10e-6, double delta = 0) {

                if (delta == 0) {
                    delta = reg_delta_0;
                }

                Eigen::Matrix<type, T, nx> X;
                dynamics.rollout_controls(x0, U, X);
                // cout << "X_0:" << endl << X << endl << "U_0:" << endl << U << endl;

                std::vector<type> Costs(max_iters);
                type cost_curr = cost.cost(X, U, X_track);
                type cost_change = cost_curr;
                type cost_new;

                Cost_Jacobians_Struct<T, nx, nu, type> cjs;
                Integrator_Jacobians_Struct<T, nx, nu, type> ijs;
                Feedback_Struct<T, nx, nu, type> fs;
                std::pair<type, type> deltaV;
                
                bool ls_success;
                int it = 0;
                while (it < max_iters && cost_change > conv_threshold) {
                    Costs.at(it) = cost_curr;

                    cost.differentiate_cost(X, U, X_track, cjs);
                    dynamics.differentiate_integrator(X, U, ijs);

                    cost_new = 0;
                    backward_pass(cjs, ijs, fs, 0, deltaV);
                    std::pair<bool, type> fp_results = forward_pass(x0, X, U, X_track, fs, deltaV, cost_curr);
                    ls_success = fp_results.first;
                    cost_change = fp_results.second;
                    if (!ls_success) {
                        cout << "Linesearch failed, exiting...";
                        break;
                    } else if (verbose >= 1) {
                        cout << "Cost Reduced!" << endl;
                    }

                    it++;
                }

                cout << "DDP done in " << it << " iterations." << endl;
                DDP::Solution<T, nx, nu, type> sol = DDP::Solution(X, U, fs.K, Costs);
                return sol;
            }

            std::pair<bool, type> forward_pass(const Eigen::Matrix<type, nx, 1>& x0,
                                Eigen::Matrix<type, T, nx>& Xbar, 
                                Eigen::Matrix<type, T-1, nu>& Ubar, 
                                const Eigen::Matrix<type, T, nx>& X_track,                        
                                const Feedback_Struct<T, nx, nu, type>& fs,
                                std::pair<type, type>& deltaV,
                                type& Cost_old) {
                type Cost_new;
                type true_reduction;
                bool ls_success = false;
                Eigen::Matrix<type, T, nx> Xnew;
                Eigen::Matrix<type, T-1, nu> Unew;

                // cout << "Cost_old:" << Cost_old << endl;
                if (verbose >= 1) {
                        cout << "Linesearch: trying alpha=";
                } 

                for (auto alpha: alphas) {
                    if (verbose >= 1) {
                        cout << alpha;
                    }

                    rollout_ddp_gains(x0, Xbar, Ubar, Xnew, Unew, fs, alpha);
                    Cost_new = cost.cost(Xnew, Unew, X_track);
                    // cout << endl << "Cost_new:" << Cost_new << endl;
                    true_reduction = Cost_old - Cost_new;
                    type expected_reduction = - alpha*(deltaV.first + alpha*deltaV.second);
                    // cout << "Expected_reduction:" << expected_reduction << endl;
                    type z = true_reduction / expected_reduction;
                    // cout << "Z:" << z << endl;
                    if (expected_reduction < 0) {
                        z = sgn(true_reduction);
                        cout << "negative expected reduction in cost!!" << endl;
                    }
                    if (z >= 0 && expected_reduction > 0) {
                        Xbar = Xnew;
                        Ubar = Unew;
                        Cost_old = Cost_new;
                        ls_success = true;
                        // cout << "Success!" << endl;
                        break;
                    }
                    
                    if (verbose >= 1) {
                        cout << ",";
                    }
                }

                if (verbose >= 1) {
                            cout << endl;
                }
                return {ls_success, true_reduction};
            }
            

            void rollout_ddp_gains(const Eigen::Matrix<type, nx, 1>& x0,
                                   const Eigen::Matrix<type, T, nx>& Xbar, 
                                   const Eigen::Matrix<type, T-1, nu>& Ubar, 
                                   Eigen::Matrix<type, T, nx>& Xnew,
                                   Eigen::Matrix<type, T-1, nu>& Unew,                                    
                                   const Feedback_Struct<T, nx, nu, type>& fs,
                                   const type alpha) {
                Xnew.row(0) = x0.transpose();
                Eigen::Matrix<type, 1, nx> deltax = Xnew.row(0)-Xbar.row(0);
                Eigen::Matrix<type, 1, nu> deltau;
                for (int t = 0; t<T-1; t++) {
                    deltau = alpha*fs.k.row(t) + (fs.K.at(t) * deltax.transpose()).transpose();
                    Unew.row(t) = Ubar.row(t) + deltau;
                    Xnew.row(t+1) = dynamics.propogate(Xnew.row(t), Unew.row(t)).transpose();
                    deltax = Xnew.row(t+1) - Xbar.row(t+1);
                }
            }

            void backward_pass(const Cost_Jacobians_Struct<T, nx, nu, type>& cjs, 
                               const Integrator_Jacobians_Struct<T, nx, nu, type>& ijs, 
                               Feedback_Struct<T, nx, nu, type>& fs, 
                               type mu, 
                               std::pair<type, type>& deltaV){

                Eigen::Matrix<type, 1, nx> Vx = cjs.Lx.row(T-1);
                Eigen::Matrix<type, nx, nx> Vxx = cjs.Lxx.back();
                for (int t = T-2; t >= 0; t--) {
                    Eigen::Matrix<type, nx, nx> fxTVxx = ijs.fx.at(t).transpose() * Vxx; 
                    Eigen::Matrix<type, nx, 1>   Qx = cjs.Lx.row(t).transpose() + ijs.fx.at(t).transpose() * Vx.transpose(); 
                    Eigen::Matrix<type, nu, 1>   Qu = cjs.Lu.row(t).transpose() + ijs.fu.at(t).transpose() * Vx.transpose(); 
                    Eigen::Matrix<type, nx, nx> Qxx = cjs.Lxx.at(t) + fxTVxx * ijs.fx.at(t); 
                    Eigen::Matrix<type, nx, nu> Qxu = cjs.Lxu.at(t) + fxTVxx * ijs.fu.at(t); 
                    Eigen::Matrix<type, nu, nu> Quu = cjs.Luu.at(t) + ijs.fu.at(t).transpose() * Vxx * ijs.fu.at(t) + mu * Eigen::Matrix<type, nu, nu>::Identity();
                    Eigen::Matrix<type, nu, nu> invQuu = Quu.inverse(); //[nu,nu]

                    fs.k.row(t) = (-invQuu * Qu).transpose(); // [1,nu]
                    fs.K.at(t) = -invQuu * Qxu.transpose(); //[nu,nx]

                    Vx = Qx + Qxu * fs.k.row(t).transpose(); //[nx,1]
                    Vxx = Qxx + Qxu * fs.K.at(t); //[nx,nx]
                    Vxx = 0.5 * (Vxx + Vxx.transpose()); //[nx,nx]

                    deltaV.first += fs.k.row(t) * Qu;
                    deltaV.second += 0.5 * fs.k.row(t) * Quu * fs.k.row(t).transpose();
                    }
            }

        private:
            void increaseReg(type &mu, type &delta) {
                delta = std::max(reg_delta_0, delta*reg_delta_0);
                mu = std::max(reg_mu_min, mu*delta);
            }

            void decreaseReg(type &mu, type &delta) {
                delta = std::min(1/reg_delta_0, delta/reg_delta_0);
                type md = mu*delta;
                if (md >= reg_mu_min) {
                    mu = md;
                } else {
                    mu = 0;
                }
            }

            const Integrator<T, nx, nu>& dynamics;
            const Cost<T, nx, nu>& cost;
            const int max_iters;
            const type conv_threshold;
            const type reg_mu_min;
            const type reg_mu_max;
            const type reg_delta_0;
            const int verbose;
            std::vector<type> alphas; 
    };
}




