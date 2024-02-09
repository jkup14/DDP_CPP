#include "../include/return_type_structs.hpp"
using std::ostream;
using std::cout;
using std::endl;

namespace DDP {
    template <int T, int nx, int nu, typename type>
    Solution<T, nx, nu, type>::Solution(Eigen::Matrix<type, T, nx> X_, 
        Eigen::Matrix<type, T-1, nu> U_,
        std::vector<Eigen::Matrix<type, nu, nx> > K_, 
        std::vector<type> J_,
        int it_)
        : X(X_), U(U_), J(J_), K(K_), it(it_) {
        if (K_.size() != T-1) {
            throw std::runtime_error("K must be size T-1");
        } 
    }
        
    template <int T, int nx, int nu, typename type>
    ostream& operator<<(ostream& os, const Solution<T, nx, nu, type>& sol) {
        os << "X:" << endl << sol.X << endl;
        os << "U:" << endl << sol.U << endl;
        os << "K:" << endl;
        for (int t = 0; t<sol.K.size(); t++) {os << sol.K.at(t) << endl;}
        os << "J:" << endl;
        for (int i = 0; i < sol.J.size(); i++) {os << sol.J.at(i) << endl;}
        os << "Iterations: " << sol.it;
        return os;
    }
}


template <int T, int nx, int nu, typename type>
Dynamics_Jacobians_Struct<T, nx, nu, type>::Dynamics_Jacobians_Struct() {
    xdotx = std::vector<Eigen::Matrix<type, nx, nx> >(T-1, Eigen::Matrix<type, nx, nx>::Zero());
    xdotu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
} 

template <int T, int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Dynamics_Jacobians_Struct<T, nx, nu, type>& djs) {
    os << "xdotx:" << endl;
    for (int t = 0; t<djs.xdotx.size(); t++) {os << djs.xdotx.at(t) << endl;}
    os << "xdotu:" << endl;
    for (int t = 0; t<djs.xdotu.size()-1; t++) {os << djs.xdotu.at(t) << endl;}
    os << djs.xdotu.at(T-2);
    return os;
}


template <int T, int nx, int nu, typename type>
Integrator_Jacobians_Struct<T, nx, nu, type>::Integrator_Jacobians_Struct() {
    fx = std::vector<Eigen::Matrix<type, nx, nx> >(T-1, Eigen::Matrix<type, nx, nx>::Zero());
    fu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
}

template <int T, int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) {
    os << "fx:" << endl;
    // os << ijs.fx[0] << endl;
    for (int t = 0; t<T-1; t++) {os << ijs.fx.at(t) << endl;}
    os << "fu:" << endl;
    for (int t = 0; t<T-2; t++) {os << ijs.fu.at(t) << endl;}
    os << ijs.fu.at(T-2);
    return os;
}


template <int T, int nx, int nu, typename type>
Cost_Jacobians_Struct<T, nx, nu, type>::Cost_Jacobians_Struct () {
    Lx = Eigen::Matrix<type, T, nx>::Zero();
    Lu = Eigen::Matrix<type, T-1, nu>::Zero();
    Lxx = std::vector<Eigen::Matrix<type, nx, nx> >(T, Eigen::Matrix<type, nx, nx>::Zero());
    Luu =  std::vector<Eigen::Matrix<type, nu, nu> >(T-1, Eigen::Matrix<type, nu, nu>::Zero());
    Lxu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
}

template <int T, int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Cost_Jacobians_Struct<T, nx, nu, type>& cjs) {
    os << "Lx:" << endl << cjs.Lx << endl;
    os << "Lu:" << endl << cjs.Lu << endl;
    os << "Lxx:" << endl;
    for (int t = 0; t<T; t++) {os << cjs.Lxx.at(t) << endl;}
    os << "Luu:" << endl;
    for (int t = 0; t<T-1; t++) {os << cjs.Luu.at(t) << endl;}
    os << "Lxu:" << endl;
    for (int t = 0; t<T-2; t++) {os << cjs.Lxu.at(t) << endl;}
    os << cjs.Lxu.at(T-2);
    return os;
}


template<int T, int nx, int nu, typename type>
Feedback_Struct<T, nx, nu, type>::Feedback_Struct () {
    k = Eigen::Matrix<type, T-1, nu>::Zero();
    K = std::vector<Eigen::Matrix<type, nu, nx> >(T-1, Eigen::Matrix<type, nu, nx>::Zero());
}

template <int T, int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Feedback_Struct<T, nx, nu, type>& fs) {
    os << "k:" << endl << fs.k << endl;
    os << "K:" << endl;
    for (int t = 0; t<fs.K.size()-1; t++) {os << fs.K.at(t) << endl;}
    os << fs.K.at(fs.K.size()-1);
    return os;
}


template<int nx, int nu, typename type>
Value_Jacobians_Struct<nx, nu, type>::Value_Jacobians_Struct(Eigen::Matrix<type, 1, nx> Lx_T, Eigen::Matrix<type, nx, nx> Lxx_T) {
    Vx = Lx_T;
    Vxx = Lxx_T;
}

template <int nx, int nu, typename type>
ostream& operator<<(std::ostream& os, const Value_Jacobians_Struct<nx, nu, type>& vjs) {
    os << "Vx:" << endl << vjs.Vx << endl;
    os << "Vxx:" << endl << vjs.Vxx;
    return os;
}

