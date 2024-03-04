#ifndef PLOT_FUNCTIONS_HPP
#define PLOT_FUNCTIONS_HPP

#include "matplot/matplot.h"
#include "Eigen/Cholesky"
#include "dynamics.hpp"
#include <iterator>
#include <thread>
#include <chrono>
#include <algorithm>
#include <math.h>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

enum Visualization_Type {two_d, cart_pole};

class Visualizer {
    public:
        Visualizer(const Visualization_Type vis_type_); 
        void animate(matplot::figure_handle &f, const Eigen::MatrixXf &X, const Eigen::MatrixXf& X_track, const int T, const float dt, std::vector<std::string> state_names);
    private:
        const Visualization_Type vis_type;
        void animate_2d_solution(matplot::figure_handle &f, const Eigen::MatrixXf &X, const Eigen::MatrixXf& X_track, const int T, const float dt, std::vector<std::string> state_names);
        // void animate_cart_pole(matplot::figure_handle &f, const Eigen::MatrixXf &X, const int T, const float dt);
        void multi_set_hold(std::vector<matplot::axes_handle>& axes, bool is_hold);
        void multi_set_aspect(std::vector<matplot::axes_handle>& axes, double aspect);
        std::vector<std::vector<double> > eigen2vec(const Eigen::MatrixXf& matrix);
        void set_axes_limits_and_names(std::vector<matplot::axes_handle> axes, std::vector<std::vector<double>>& min_and_maxes, float end_time, std::vector<std::string> state_names);
        std::vector<std::vector<double>> get_min_and_max(std::vector<std::vector<double> >& vecs);
        std::vector<double> get_time_vec(int T, float dt);
        void plot_state_trajectories_and_targets(std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> >& vecs, std::vector<std::vector<double> >& vecs_track, std::vector<double>& times);
        void get_time_lines(std::vector<matplot::line_handle>& time_lines, std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> > axis_limits);
        void set_up_state_axes_and_get_time_lines(std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> >& vecs, std::vector<std::vector<double> >& vecs_track, std::vector<matplot::line_handle>& time_lines, int T, float dt, std::vector<std::string> state_names);
};



#endif