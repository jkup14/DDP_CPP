#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
using std::ostream;
using std::endl;

template <int T, int nx, int nu, typename type>
constexpr int Dynamics<T, nx, nu, type>::getNx() const {return nx;}

template <int T, int nx, int nu, typename type>
constexpr int Dynamics<T, nx, nu, type>::getNu() const {return nu;}


template <int T, typename type>
Eigen::Matrix<type, 2, 1> SingleIntegrator<T, type>::xdot(const Eigen::Matrix<type, 2, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
    return Eigen::Matrix<type, 2, 1>(u);
}

template <int T, typename type>
void SingleIntegrator<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const {
    djs.xdotx = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Zero());
    djs.xdotu = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Identity());
}

template <int T, typename type>
Eigen::Matrix<type, 4, 1> DoubleIntegrator<T, type>::xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
            return (Eigen::Matrix<type, 4, 1>()<< x(Eigen::seq(2,3)), u).finished();
            return Eigen::Matrix<type, 4, 1>::Zero();
        }

template <int T, typename type>
void DoubleIntegrator<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 4>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 4, 2>& djs) const {
    Eigen::Matrix<type, 4, 4> xdotx_t = Eigen::Matrix<type, 4, 4>::Zero();
    xdotx_t(0,2) = 1;
    xdotx_t(1,3) = 1;
    djs.xdotx = std::vector<Eigen::Matrix<type, 4, 4> >(T-1, xdotx_t);
    Eigen::Matrix<type, 4, 2> xdotu_t = Eigen::Matrix<type, 4, 2>::Zero();
    xdotu_t(2,0) = 1;
    xdotu_t(3,1) = 1;
    djs.xdotu = std::vector<Eigen::Matrix<type, 4, 2> >(T-1, xdotu_t);
}

template <int T, typename type>
DifferentialDrive<T, type>::DifferentialDrive():wheel_radius(0.2),width(0.2) {}

template <int T, typename type>
Eigen::Matrix<type, 3, 1> DifferentialDrive<T, type>::xdot(const Eigen::Matrix<type, 3, 1>& x, const Eigen::Matrix<type, 2, 1>& u) const{
    return (Eigen::Matrix<type, 3, 1>() << 
    cos(x(2)) * (u(0)+u(1)),
    sin(x(2)) * (u(0)+u(1)),
    (u(0)-u(1)) / width
    ).finished()*wheel_radius/2;
}

template <int T, typename type>
void DifferentialDrive<T, type>::differentiate_dynamics(const Eigen::Matrix<type, T, 3>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 3, 2>& djs) const {
    // djs is already initialized to zeros when constructed
    for (int t = 0; t<T-1; t++) {
        type sin_angle = sin(X(t,2));
        type cos_angle = cos(X(t,2));

        djs.xdotx.at(t) = (((Eigen::Matrix<type, 3, 3>()<< 
        0,0,-sin_angle,
        0,0,cos_angle,
        0,0,0
        ).finished())
        *wheel_radius*(U(t,0)+U(t,1))/2);

        djs.xdotu.at(t) = (((Eigen::Matrix<type, 3, 2>()<< 
        cos_angle, cos_angle,
        sin_angle, sin_angle,
        1/width, -1/width).finished())
        *wheel_radius/2);
    }
}

// template <int T, typename type=float>
// class CartPole : public Dynamics<T, 4, 1, type> {
//     public:  
//         CartPole() {
//             mass_pole = 0.05;
//             mass_cart = 1;
//             gravity = 9.8;
//             pole_length = 2;
//         } 
//         Eigen::Matrix<type, 4, 1> xdot(const Eigen::Matrix<type, 4, 1>& x, const Eigen::Matrix<type, 1, 1>& u) const{
//             // state is [displacement, angle, velocity, angular_velocity]
//             type sin_angle = sin(x(1));
//             type sin_angle_sqr = pow(sin_angle, 2);
//             type cos_angle = cos(x(1));
//             type ang_velocity_sqr = pow(x(3),2);
//             return (Eigen::Matrix<type, 4, 1>() << 
//                 x(2),
//                 x(3),
//                 (mass_pole*sin_angle
//                     * (pole_length * ang_velocity_sqr * gravity * cos_angle)
//                     + u) 
//                     / (mass_cart + mass_pole * sin_angle_sqr),
//                 (-mass_pole * pole_length * ang_velocity_sqr * cos_angle * sin_angle 
//                     - (mass_cart + mass_pole) * gravity * sin_angle
//                     - cos_angle * u)
//                     / (pole_length * (mass_cart + mass_pole * sin_angle_sqr))
//                 ).finished();
//         }

//         void differentiate_dynamics(const Eigen::Matrix<type, T, 2>& X, const Eigen::Matrix<type, T-1, 2>& U, Dynamics_Jacobians_Struct<T, 2, 2>& djs) const {
//             djs.xdotx = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Zero());
//             djs.xdotu = std::vector<Eigen::Matrix<type, 2, 2> >(T-1, Eigen::Matrix<type, 2, 2>::Identity());
//         }

//     private:
//         type mass_pole;
//         type mass_cart;
//         type gravity;
//         type pole_length;
    
// };


#endif
