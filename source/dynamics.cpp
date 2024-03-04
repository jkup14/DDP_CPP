#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <vector>
using std::ostream;
using std::endl;

template <int T, int nx, int nu, typename type>
Dynamics_Abstract<T, nx, nu, type>::Dynamics_Abstract(std::vector<std::string> state_names_):state_names(state_names_) {}

template <int T, int nx, int nu, typename type>
constexpr int Dynamics_Abstract<T, nx, nu, type>::getNx() const {return nx;}

template <int T, int nx, int nu, typename type>
constexpr int Dynamics_Abstract<T, nx, nu, type>::getNu() const {return nu;}

template <int T, int nx, int nu, typename type>
std::vector<std::string> Dynamics_Abstract<T, nx, nu, type>::get_state_names() {return state_names;}


#endif
