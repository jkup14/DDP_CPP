#ifndef HELPER_FUNCS_HPP
#define HELPER_FUNCS_HPP

#include <ctime>
#include <Eigen/Cholesky>

template <typename T> 
int sgn(T val);

template <typename T>
T random_range(T min, T max);

template<typename T, int r, int c>
Eigen::Matrix<T, r, c> random_matrix(Eigen::Matrix<T, r, c> min, Eigen::Matrix<T, r, c> max);

#include "../source/helper_funcs.cpp"

#endif