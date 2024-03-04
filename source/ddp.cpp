#ifndef DDP_CPP
#define DDP_CPP

#include <vector>
#include "../include/return_type_structs.hpp"
#include "../include/integrator.hpp"
#include "../include/cost.hpp"
#include <Eigen/LU>
#include "../include/helper_funcs.hpp"
#include "../include/timer.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration;

template <int T, int nx, int nu, typename type>
DDP::DDP_Solver<T, nx, nu, type>::DDP_Solver(const Integrator<T, nx, nu, type>& dynamics_,
    const Cost<T, nx, nu, type>& cost_,
    const solver_args<type> args) 
    :   dynamics(dynamics_), 
        cost(cost_),
        max_iters(args.max_iters), 
        conv_threshold(args.conv_threshold),
        reg_mu_min(args.reg_mu_min),
        reg_mu_max(args.reg_mu_max),
        reg_delta_factor(args.reg_delta_factor),
        verbose(args.verbose),
        toggle_ls(args.toggle_ls),
        toggle_reg(args.toggle_reg),
        timing_level(args.timing_level) {

    
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
template <int T, int nx, int nu, typename type>
DDP::Solution<T, nx, nu, type> DDP::DDP_Solver<T, nx, nu, type>::solve(const Eigen::Matrix<type, nx, 1>& x0,
        Eigen::Matrix<type, T-1, nu>& U,
        Eigen::Matrix<type, T, nx>& X_track, type mu, type delta) {
    
    Timer total_Timer;
    Timer fp_Timer;
    Timer bp_Timer;
    Timer dCost_Timer;
    Timer dIntegrator_Timer;
    long long total_Time;
    if (timing_level >= 1) {total_Timer.Start();}


    if (delta == 0) {
        delta = reg_delta_factor;
    }


    Eigen::Matrix<type, T, nx> X;
    dynamics.rollout_controls(x0, U, X);

    std::vector<type> Costs(max_iters);
    type cost_curr = cost.cost(X, U, X_track);
    type cost_change = cost_curr;

    Cost_Jacobians_Struct<T, nx, nu, type> cjs;
    Integrator_Jacobians_Struct<T, nx, nu, type> ijs;
    Feedback_Struct<T, nx, nu, type> fs;
    std::pair<type, type> deltaV;
    
    bool forwardpass_success;
    bool cost_continuation_criteria = true;
    bool regularization_possible;
    int it = 0;
    while (it < max_iters && cost_continuation_criteria) {
        Costs.at(it) = cost_curr;

        if (timing_level >= 2) {dCost_Timer.Start();}
        cost.differentiate_cost(X, U, X_track, cjs);
        if (timing_level >= 2) {
            dCost_Timer.Stop();
            dIntegrator_Timer.Start();
        }
        dynamics.differentiate_integrator(X, U, ijs);
        if (timing_level >= 2) {
            dIntegrator_Timer.Stop();
            bp_Timer.Start();
        }

        backward_pass(cjs, ijs, fs, mu, deltaV);
        if (timing_level >= 2) {
            bp_Timer.Stop();
            fp_Timer.Start();
        }
        std::pair<bool, type> fp_results = forward_pass(x0, X, U, X_track, fs, deltaV, cost_curr);
        if (timing_level >= 2) {fp_Timer.Stop();}

        forwardpass_success = fp_results.first;
        cost_change = fp_results.second;

        // If change in cost is <= 0, we keep going
        cost_continuation_criteria = forwardpass_success ? cost_change >= conv_threshold : true;

        // Decide on forward pass success by continuing, applying regularization, or exiting
        it++;
        
        if (forwardpass_success) {
            if (toggle_reg) {decreaseReg(mu, delta);}
            printIteration(it, cost_curr, mu);
            continue; // go to next iteration
        } 
        // go to next iteration if regularization is successful
        if (toggle_reg && increaseReg(mu, delta)) {
            continue;
        } 
        if (verbose >= 1) {cout << "Solution not found, exiting..." << endl;}
        break;
    }
    if (timing_level >= 1) {total_Time = total_Timer.Stop();}

    if (it >= max_iters && verbose >= 1) { cout << "Max iterations reached" << endl;}

    if (timing_level >= 2) {
        cout << "FP Average and Total Time: " << fp_Timer.avg_Time() << "µs, " << fp_Timer.total_Time() << "µs" << endl;
        cout << "BP Average and Total Time: " << bp_Timer.avg_Time() << "µs, " << bp_Timer.total_Time() << "µs" << endl;
        cout << "Cost Differentiation and Total Time: " << dCost_Timer.avg_Time() << "µs, " << dCost_Timer.total_Time()<< "µs" << endl;
        cout << "Integrator Differentiation and Total Time: " << dIntegrator_Timer.avg_Time() << "µs, " << dIntegrator_Timer.total_Time()<< "µs" << endl;
    }

    DDP::Solution<T, nx, nu, type> sol = DDP::Solution<T, nx, nu, type>(X, U, fs.K, std::vector<type>(Costs.begin(), Costs.begin()+it), it-1, total_Time*(0.001));
    return sol;
}

template <int T, int nx, int nu, typename type>
void DDP::DDP_Solver<T, nx, nu, type>::printIteration (const int& it, const type& cost_curr, const type& mu) {
    if (verbose >= 2) {
        cout << "Cost Reduced! ";
        if (toggle_reg && mu > 0) {cout << "Mu is now " << mu << endl;}
        else { cout << endl;}
    }
    if (verbose >= 1) {cout << "Iteration " << it << " Cost: " <<  cost_curr << endl;}
}

template <int T, int nx, int nu, typename type>
std::pair<bool, type> DDP::DDP_Solver<T, nx, nu, type>::forward_pass(const Eigen::Matrix<type, nx, 1>& x0,
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

    if (verbose >= 2) {
            cout << "Linesearch: trying alpha=";
    } 

    for (auto alpha: alphas) {
        if (verbose >= 2) {
            cout << alpha;
        }

        rollout_ddp_gains(x0, Xbar, Ubar, Xnew, Unew, fs, alpha);
        Cost_new = cost.cost(Xnew, Unew, X_track);
        true_reduction = Cost_old - Cost_new;
        type expected_reduction = - alpha*(deltaV.first + alpha*deltaV.second);
        type z = true_reduction / expected_reduction;
        if (expected_reduction < 0) {
            z = sgn(true_reduction);
            cout << "negative expected reduction in cost!!" << endl;
        }
        if (z >= 0 && expected_reduction > 0) {
            Xbar = Xnew;
            Ubar = Unew;
            Cost_old = Cost_new;
            ls_success = true;
            break;
        }
        
        if (verbose >= 2) {
            cout << ",";
        }
    }

    if (verbose >= 2) {
                cout << " ...";
    }
    return {ls_success, true_reduction};
}

template <int T, int nx, int nu, typename type>
void DDP::DDP_Solver<T, nx, nu, type>::rollout_ddp_gains(const Eigen::Matrix<type, nx, 1>& x0,
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

template <int T, int nx, int nu, typename type>
void DDP::DDP_Solver<T, nx, nu, type>::backward_pass(const Cost_Jacobians_Struct<T, nx, nu, type>& cjs, 
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
        auto deltaV_second = 0.5 * fs.k.row(t) * Quu * fs.k.row(t).transpose();
        deltaV.second += deltaV_second(0);
        }
}

template <int T, int nx, int nu, typename type>
bool DDP::DDP_Solver<T, nx, nu, type>::increaseReg(type &mu, type &delta) {
    delta = std::max(reg_delta_factor, delta*reg_delta_factor);
    mu = std::max(reg_mu_min, mu*delta);
    if (verbose >= 2) {
        cout << "Failed, Mu is now " << mu << endl;
    }
    return mu <= reg_mu_max;
}

template <int T, int nx, int nu, typename type>
void DDP::DDP_Solver<T, nx, nu, type>::decreaseReg(type &mu, type &delta) {
    delta = std::min(1/reg_delta_factor, delta/reg_delta_factor);
    type md = mu*delta;
    mu = (md >= reg_mu_min) ? md : 0;
}

#endif
    





