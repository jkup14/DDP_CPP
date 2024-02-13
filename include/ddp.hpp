#ifndef DDP_HPP
#define DDP_HPP

#include <vector>
#include "../include/return_type_structs.hpp"
#include "../include/integrator.hpp"
#include "../include/cost.hpp"
#include <Eigen/LU>
#include "../include/sign_func.hpp"


namespace DDP {
    template <typename type=float>
    struct solver_args {
        int max_iters = 500;
        type conv_threshold = 1e-3;
        type reg_mu_min = 10e-6;
        type reg_mu_max = 1e10;
        type reg_delta_factor = 2;
        type ls_alpha_0 = 1;
        type ls_alpha_f = 10e-3;
        int ls_alpha_len = 11;
        bool toggle_ls = true;
        bool toggle_reg = true;
        int timing_level = 0;
        int verbose = 0;
    };

    template <int T, int nx, int nu, typename type=float>
    class DDP_Solver {
        public:
            DDP_Solver(const Integrator<T, nx, nu, type>& dynamics_,
                const Cost<T, nx, nu, type>& cost_,
                const solver_args<type> args = solver_args<type>());

            DDP::Solution<T, nx, nu, type> solve(const Eigen::Matrix<type, nx, 1>& x0,
                Eigen::Matrix<type, T-1, nu>& U,
                Eigen::Matrix<type, T, nx>& X_track, type mu = 10e-6, type delta = 0) ;

            void printIteration (const int& it, const type& cost_curr, const type& mu);

            std::pair<bool, type> forward_pass(const Eigen::Matrix<type, nx, 1>& x0,
                Eigen::Matrix<type, T, nx>& Xbar, 
                Eigen::Matrix<type, T-1, nu>& Ubar, 
                const Eigen::Matrix<type, T, nx>& X_track,                        
                const Feedback_Struct<T, nx, nu, type>& fs,
                std::pair<type, type>& deltaV,
                type& Cost_old) ;
            

            void rollout_ddp_gains(const Eigen::Matrix<type, nx, 1>& x0,
                                   const Eigen::Matrix<type, T, nx>& Xbar, 
                                   const Eigen::Matrix<type, T-1, nu>& Ubar, 
                                   Eigen::Matrix<type, T, nx>& Xnew,
                                   Eigen::Matrix<type, T-1, nu>& Unew,                                    
                                   const Feedback_Struct<T, nx, nu, type>& fs,
                                   const type alpha) ;

            void backward_pass(const Cost_Jacobians_Struct<T, nx, nu, type>& cjs, 
                               const Integrator_Jacobians_Struct<T, nx, nu, type>& ijs, 
                               Feedback_Struct<T, nx, nu, type>& fs, 
                               type mu, 
                               std::pair<type, type>& deltaV);

        private:
            bool increaseReg(type &mu, type &delta) ;

            void decreaseReg(type &mu, type &delta) ;

            const Integrator<T, nx, nu>& dynamics;
            const Cost<T, nx, nu>& cost;
            const int max_iters;
            const type conv_threshold;
            const type reg_mu_min;
            const type reg_mu_max;
            const type reg_delta_factor;
            const int verbose;
            std::vector<type> alphas; 
            const bool toggle_ls;
            const bool toggle_reg;
            const int timing_level;
    };
}

#include "../source/ddp.cpp"

#endif


