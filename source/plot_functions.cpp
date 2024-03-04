#include "../include/plot_functions.hpp"

Visualizer::Visualizer(Visualization_Type vis_type_):vis_type(vis_type_){}

void Visualizer::animate(matplot::figure_handle &f, const Eigen::MatrixXf &X, const Eigen::MatrixXf& X_track, const int T, const float dt, std::vector<std::string> state_names) {
    switch (vis_type) {
        case two_d:
            animate_2d_solution(f, X, X_track, T, dt, state_names); break;
        // case cart_pole:
        //     animate_cart_pole(f, X, T, dt); break;
    }
}

void Visualizer::animate_2d_solution(matplot::figure_handle &f, const Eigen::MatrixXf &X, const Eigen::MatrixXf& X_track, const int T, const float dt, std::vector<std::string> state_names) {
    std::vector<std::vector<double> > vecs = eigen2vec(X);
    std::vector<std::vector<double> > vecs_track = eigen2vec(X_track);

    int w = 2;
    int h = ceil((vecs.size()+1)/2.0);
    matplot::tiledlayout(w,h);
    std::vector<matplot::axes_handle> axes, state_axes;
    std::generate_n(std::back_inserter(axes), vecs.size()+1, []() {return matplot::nexttile();});
    
    auto traj_ax = axes[0];
    state_axes = std::vector<matplot::axes_handle>(axes.begin()+1, axes.end());

    multi_set_hold(axes, true);
    multi_set_aspect(axes, 1);

    traj_ax->title("Trajectory");
    traj_ax->xlabel("x");
    traj_ax->ylabel("y");
    

    std::vector<matplot::line_handle> time_lines;
    set_up_state_axes_and_get_time_lines(state_axes, vecs, vecs_track, time_lines, T, dt, state_names);

    matplot::xlim(traj_ax, {*std::min_element(vecs[0].begin(), vecs[0].end())-0.5, *std::max_element(vecs[0].begin(), vecs[0].end())+0.5});
    matplot::ylim(traj_ax, {*std::min_element(vecs[1].begin(), vecs[1].end())-0.5, *std::max_element(vecs[1].begin(), vecs[1].end())+0.5});
    matplot::line_handle trajectory = traj_ax->plot(std::vector<double>((vecs[0]).begin(), vecs[0].begin()+1), std::vector<double>(vecs[1].begin(), vecs[1].begin()+1));

    traj_ax->draw();

    using namespace std::chrono_literals;
    double ani_time = 5;
    double ani_dt = ani_time/T;
    for (int t=1; t < T; t++) {

        trajectory->x_data(std::vector<double>(vecs[0].begin(), vecs[0].begin()+t));
        trajectory->y_data(std::vector<double>(vecs[1].begin(), vecs[1].begin()+t));

        for (auto& time_line: time_lines) {
            time_line->x_data(std::vector<double>{dt*t,dt*t});
        }
        
        traj_ax->draw();
        std::this_thread::sleep_for(std::chrono::duration<float>(ani_dt));
    }
}

void Visualizer::multi_set_hold(std::vector<matplot::axes_handle>& axes, bool is_hold) {
    for (auto& ax: axes) {
        ax->hold(is_hold);
    }
}
void Visualizer::multi_set_aspect(std::vector<matplot::axes_handle>& axes, double aspect) {
    for (auto& ax: axes) {
        ax->axes_aspect_ratio(aspect);
    }
}

std::vector<std::vector<double> > Visualizer::eigen2vec(const Eigen::MatrixXf& matrix) {
    std::vector<std::vector<double> > vecs;
    for (int i =0 ; i < matrix.cols(); i++) {
        vecs.push_back(std::vector<double>(matrix.col(i).data(), matrix.col(i).data()+matrix.rows()));
    }
    return vecs;
}

void Visualizer::set_axes_limits_and_names(std::vector<matplot::axes_handle> state_axes, std::vector<std::vector<double>>& min_and_maxes, float end_time, std::vector<std::string> state_names) {
    for (int i = 0; i < state_axes.size(); i++) {
        matplot::xlim(state_axes[i], {0, end_time});
        matplot::ylim(state_axes[i], {min_and_maxes[i][0], min_and_maxes[i][1]});
        state_axes[i]->xlabel("Time (s)");
        state_axes[i]->ylabel(state_names[i]);
    }
}

std::vector<std::vector<double>> Visualizer::get_min_and_max(std::vector<std::vector<double> >& vecs) {
    std::vector<std::vector<double>> min_and_maxes;
    for (int i = 0; i < vecs.size(); i++) {
        std::vector<double> min_and_max = {*std::min_element(vecs[i].begin(), vecs[i].end())-0.1, *std::max_element(vecs[i].begin(), vecs[i].end())+0.1};
        min_and_maxes.push_back(min_and_max);
    }
    return min_and_maxes;
}

std::vector<double> Visualizer::get_time_vec(int T, float dt) {
    std::vector<double> times;
    std::generate_n(
        std::back_inserter(times),
        T,
        [t = -dt, dt]() mutable {t+=dt; return t; }
    );
    return times;
}

void Visualizer::plot_state_trajectories_and_targets(std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> >& vecs, std::vector<std::vector<double> >& vecs_track, std::vector<double>& time_vec) {
    for (int i = 0; i< state_axes.size(); i++) {
        state_axes[i]->plot(time_vec, vecs[i]);
        state_axes[i]->plot(time_vec, vecs_track[i], "r--");
    }
}

void Visualizer::get_time_lines(std::vector<matplot::line_handle>& time_lines, std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> > axis_limits) {
    std::generate_n(
        std::back_inserter(time_lines),
        state_axes.size(),
        [i=-1, state_axes, axis_limits]() mutable { i+=1; return state_axes[i]->plot(std::vector<double>{0,0}, std::vector<double>{axis_limits[i][0],axis_limits[i][1]}); }
    );
}

void Visualizer::set_up_state_axes_and_get_time_lines(std::vector<matplot::axes_handle>& state_axes, std::vector<std::vector<double> >& vecs, std::vector<std::vector<double> >& vecs_track, std::vector<matplot::line_handle>& time_lines, int T, float dt, std::vector<std::string> state_names) {
    std::vector<std::vector<double>> axis_limits = get_min_and_max(vecs);
    std::vector<double> time_vec = get_time_vec(T, dt);
    set_axes_limits_and_names(state_axes, axis_limits, T*dt, state_names);
    plot_state_trajectories_and_targets(state_axes, vecs, vecs_track, time_vec);
    get_time_lines(time_lines, state_axes, axis_limits);
}