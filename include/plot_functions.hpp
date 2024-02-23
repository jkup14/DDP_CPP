#ifndef PLOT_FUNCTIONS_HPP
#define PLOT_FUNCTIONS_HPP

#include "matplot/matplot.h"
#include "Eigen/Cholesky"
#include <iterator>
#include <thread>
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

void multi_set_hold(std::vector<matplot::axes_handle>& axes, bool is_hold);
void multi_set_aspect(std::vector<matplot::axes_handle>& axes, double aspect);
void animate_2d_solution(matplot::figure_handle &f, const Eigen::MatrixXf &X, const int T, const float dt);

#endif