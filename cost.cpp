#include <Eigen/Cholesky>
using std::cout;
using std::endl;


// TODO: maybe store X_diff for performance isntead of calculating twice

template <int T, int nx, int nu, typename type=float>
class Cost {
    public:
        virtual type cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track) const = 0;

        virtual void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Jacobians_Struct<T, nx, nu, type>& cjs) const = 0;
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

        void differentiate_cost(const Eigen::Matrix<type, T, nx>& X, const Eigen::Matrix<type, T-1, nu>& U, const Eigen::Matrix<type, T, nx> X_track, Cost_Jacobians_Struct<T, nx, nu, type>& cjs) const {
            const Eigen::Matrix<type, T, nx> X_diff = X-X_track;
            for (int t = 0; t<T-1; t++) {
                cjs.Lx.row(t) = X_diff.row(t) * Q;
                cjs.Lu.row(t) = U.row(t) * R;
                cjs.Lxx.at(t) = Q;
                cjs.Luu.at(t) = R;
                cjs.Lxu.at(t) = Eigen::Matrix<type, nx, nu>::Zero();
            }
            // cout << "X last " << X.row(T-1) << endl;
            // cout << "X_track last" << X_track.row(T-1) << endl;
            // cout << "X_diff last " << X_diff.row(T-1) << endl;
            // cout << "Qf" << Qf << endl;
            cjs.Lx.row(T-1) = X_diff.row(T-1) * Qf;
            cjs.Lxx.at(T-1) = Qf;
        }

    private:
        const Eigen::Matrix<type, nx, nx>& Q;
        const Eigen::Matrix<type, nu, nu>& R;
        const Eigen::Matrix<type, nx, nx>& Qf;
};

