#include <Eigen/Cholesky>
using std::cout;
using std::endl;


// TODO: maybe store X_diff for performance isntead of calculating twice

template <int T, int nx, int nu, typename type=float>
struct Cost_Derivatives_Struct {
    Cost_Derivatives_Struct () {
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
    friend ostream& operator<<(ostream& os, const Cost_Derivatives_Struct<T_, nx_, nu_, type_>& cds);
};

template <int T, int nx, int nu, typename type=float>
ostream& operator<<(ostream& os, const Cost_Derivatives_Struct<T, nx, nu, type>& cds) {
    os << "Lx:" << endl << cds.Lx << endl;
    os << "Lu:" << endl << cds.Lu << endl;
    os << "Lxx:" << endl;
    for (int t = 0; t<T; t++) {os << cds.Lxx.at(t) << endl;}
    os << "Luu:" << endl;
    for (int t = 0; t<T-1; t++) {os << cds.Luu.at(t) << endl;}
    os << "Lxu:" << endl;
    for (int t = 0; t<T-2; t++) {os << cds.Lxu.at(t) << endl;}
    os << cds.Lxu.at(T-2);
    return os;
}

template <int T, int nx, int nu, typename type=float>
class Cost {
    virtual type cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const = 0;

    virtual void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Derivatives_Struct<T, nx, nu, type>& cds) const = 0;
};

template <int T, int nx, int nu, typename type=float>
class QuadraticCost : public Cost<T, nx, nu, type> {
    public:
        QuadraticCost (const Eigen::Matrix<type, nx, nx>& Q_, const Eigen::Matrix<type, nu, nu>& R_, const Eigen::Matrix<type, nx, nx>& Qf_):
            Q(Q_), R(R_), Qf(Qf_) {}

        type cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const {
            type sum = 0;
            const Eigen::Matrix<type, T, nx> X_diff = X-X_track;
            for (int t = 0; t<T-1; t++) {
                sum += X_diff.row(t) * Q * X_diff.row(t).transpose();
                sum += U.row(t) * R * U.row(t).transpose();
            }
            sum += X_diff.row(T-1) * Qf * X_diff.row(T-1).transpose();
            return sum;
        }

        void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Derivatives_Struct<T, nx, nu, type>& cds) const {
            const Eigen::Matrix<type, T, nx> X_diff = X-X_track;
            for (int t = 0; t<T-1; t++) {
                cds.Lx.row(t) = X_diff.row(t) * Q;
                cds.Lu.row(t) = U.row(t) * R;
                cds.Lxx.at(t) = Q;
                cds.Luu.at(t) = R;
                cds.Lxu.at(t) = Eigen::Matrix<type, nx, nu>::Zero();
            }
            cds.Lx.row(T-1) = X_diff.row(T-1) * Qf;
            cds.Lxx.at(T-1) = Qf;
        }

    private:
        const Eigen::Matrix<type, nx, nx>& Q;
        const Eigen::Matrix<type, nu, nu>& R;
        const Eigen::Matrix<type, nx, nx>& Qf;
};

