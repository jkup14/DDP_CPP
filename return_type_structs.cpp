#include <iostream>
#include <Eigen/Cholesky>
using std::ostream;
using std::cout;
using std::endl;


namespace DDP {
    template <int T, int nx, int nu, typename=float>
    struct Solution {
        Eigen::Matrix<float, T, nx> X;
        Eigen::Matrix<float, T-1, nu> U;
        std::vector<Eigen::Matrix<float, nu, nx> > K;
        std::vector<float> J;
        int it;
        Solution(Eigen::Matrix<float, T, nx> X_, 
                 Eigen::Matrix<float, T-1, nu> U_,
                 std::vector<Eigen::Matrix<float, nu, nx> > K_, 
                 std::vector<float> J_,
                 int it_)
                 : X(X_), U(U_), J(J_), it(it_) {
                    if (K_.size() != T-1) {
                        throw std::runtime_error("K must be size T-1");
                    } 
                    K = K_;
                 }

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Solution<T_, nx_, nu_, type_>& vjs);
    };

    template <int T, int nx, int nu, typename type>
    ostream& operator<<(ostream& os, const Solution<T, nx, nu, type>& sol) {
        os << "X:" << endl << sol.X << endl;
        os << "U:" << endl << sol.U << endl;
        os << "K:" << endl;
        for (int t = 0; t<sol.K.size(); t++) {os << sol.K.at(t) << endl;}
        os << "J:" << endl;
        for (int i = 0; i < sol.J.size()-1; i++) {os << sol.J.at(i) << endl;}
        os << sol.J.back();
        return os;
    }
}

template <int T, int nx, int nu, typename type=float>
struct Dynamics_Jacobians_Struct {
    Dynamics_Jacobians_Struct() {
        xdotx = std::vector<Eigen::Matrix<type, nx, nx> >(T-1, Eigen::Matrix<type, nx, nx>::Zero());
        xdotu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
    } 
    std::vector<Eigen::Matrix<type, nx, nx> > xdotx;
    std::vector<Eigen::Matrix<type, nx, nu> > xdotu;

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Dynamics_Jacobians_Struct<T_, nx_, nu_, type_>& djs);
};

template <int T, int nx, int nu, typename type=float>
ostream& operator<<(ostream& os, const Dynamics_Jacobians_Struct<T, nx, nu, type>& djs) {
    os << "xdotx:" << endl;
    for (int t = 0; t<djs.xdotx.size(); t++) {os << djs.xdotx.at(t) << endl;}
    os << "xdotu:" << endl;
    for (int t = 0; t<djs.xdotu.size()-1; t++) {os << djs.xdotu.at(t) << endl;}
    os << djs.xdotu.at(T-2);
    return os;
}

template <int T, int nx, int nu, typename type=float>
struct Integrator_Jacobians_Struct {
    Integrator_Jacobians_Struct() {
        fx = std::vector<Eigen::Matrix<type, nx, nx> >(T-1, Eigen::Matrix<type, nx, nx>::Zero());
        fu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
    }
    std::vector<Eigen::Matrix<type, nx, nx> > fx;
    std::vector<Eigen::Matrix<type, nx, nu> > fu;

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Integrator_Jacobians_Struct<T_, nx_, nu_, type_>& ijs);
};

template <int T, int nx, int nu, typename type=float>
ostream& operator<<(ostream& os, const Integrator_Jacobians_Struct<T, nx, nu, type>& ijs) {
    os << "fx:" << endl;
    // os << ijs.fx[0] << endl;
    for (int t = 0; t<T-1; t++) {os << ijs.fx.at(t) << endl;}
    os << "fu:" << endl;
    for (int t = 0; t<T-2; t++) {os << ijs.fu.at(t) << endl;}
    os << ijs.fu.at(T-2);
    return os;
}

template <int T, int nx, int nu, typename type=float>
struct Cost_Jacobians_Struct {
    Cost_Jacobians_Struct () {
        Lx = Eigen::Matrix<type, T, nx>::Zero();
        Lu = Eigen::Matrix<type, T-1, nu>::Zero();
        Lxx = std::vector<Eigen::Matrix<type, nx, nx> >(T, Eigen::Matrix<type, nx, nx>::Zero());
        Luu =  std::vector<Eigen::Matrix<type, nu, nu> >(T-1, Eigen::Matrix<type, nu, nu>::Zero());
        Lxu = std::vector<Eigen::Matrix<type, nx, nu> >(T-1, Eigen::Matrix<type, nx, nu>::Zero());
    }
    Eigen::Matrix<type, T, nx> Lx;
    Eigen::Matrix<type, T-1, nu> Lu;
    std::vector<Eigen::Matrix<type, nx, nx> > Lxx;
    std::vector<Eigen::Matrix<type, nu, nu> > Luu;
    std::vector<Eigen::Matrix<type, nx, nu> > Lxu;

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Cost_Jacobians_Struct<T_, nx_, nu_, type_>& cjs);
};

template <int T, int nx, int nu, typename type=float>
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

template<int T, int nx, int nu, typename type=float>
struct Feedback_Struct {
    Feedback_Struct () {
        k = Eigen::Matrix<type, T-1, nu>::Zero();
        K = std::vector<Eigen::Matrix<type, nu, nx> >(T-1, Eigen::Matrix<type, nu, nx>::Zero());
    }
    Eigen::Matrix<type, T-1, nu> k;
    std::vector<Eigen::Matrix<type, nu, nx> > K;

    template <int T_, int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Feedback_Struct<T_, nx_, nu_, type_>& fs);
};

template <int T, int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Feedback_Struct<T, nx, nu, type>& fs) {
    os << "k:" << endl << fs.k << endl;
    os << "K:" << endl;
    for (int t = 0; t<fs.K.size()-1; t++) {os << fs.K.at(t) << endl;}
    os << fs.K.at(fs.K.size()-1);
    return os;
}

template<int nx, int nu, typename type=float>
struct Value_Jacobians_Struct {
    Value_Jacobians_Struct(Eigen::Matrix<type, 1, nx> Lx_T, Eigen::Matrix<type, nx, nx> Lxx_T) {
        Vx = Lx_T;
        Vxx = Lxx_T;
    }
    Eigen::Matrix<type, nx, 1> Vx;
    Eigen::Matrix<type, nx, nx> Vxx;

    template <int nx_, int nu_, typename type_>
    friend ostream& operator<<(ostream& os, const Value_Jacobians_Struct<nx_, nu_, type_>& vjs);
};

template <int nx, int nu, typename type>
ostream& operator<<(ostream& os, const Value_Jacobians_Struct<nx, nu, type>& vjs) {
    os << "Vx:" << endl << vjs.Vx << endl;
    os << "Vxx:" << endl << vjs.Vxx;
    return os;
}

