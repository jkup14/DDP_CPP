#ifndef HELPER_FUNCS_CPP
#define HELPER_FUNCS_CPP

#include <Eigen/Cholesky>


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T> 
T random_range(T min, T max) {
    rand();
    return (T)rand()/RAND_MAX * (max-min) + min;
}

template<typename T, int r, int c>
Eigen::Matrix<T, r, c> random_matrix(Eigen::Matrix<T, r, c> min, Eigen::Matrix<T, r, c> max) {
    rand(); // To encourage divergence for similar seeds
    auto random_matrix = Eigen::Matrix<T, r, c>::Random().cwiseAbs();
    auto range_matrix = max-min;
    return (random_matrix.array()*range_matrix.array()) + min.array();
}
 

#endif