#include "../include/plot_functions.hpp"

void multi_set_hold(std::vector<matplot::axes_handle>& axes, bool is_hold) {
    for (auto& ax: axes) {
        ax->hold(is_hold);
    }
}
void multi_set_aspect(std::vector<matplot::axes_handle>& axes, double aspect) {
    for (auto& ax: axes) {
        ax->axes_aspect_ratio(aspect);
    }
}

void animate_2d_solution(matplot::figure_handle &f, const Eigen::MatrixXf &X, const int T, const float dt) {
    std::vector<double> x(X.col(0).data(), X.col(0).data() + X.col(0).rows() * X.col(0).cols());
    std::vector<double> y(X.col(1).data(), X.col(1).data() + X.col(1).rows() * X.col(1).cols());

    matplot::tiledlayout(2,2);
    std::vector<matplot::axes_handle> axes;
    std::generate_n(std::back_inserter(axes), 3, []() {return matplot::nexttile();});
    
    auto ax1 = axes[0];
    auto ax_X = axes[1];
    auto ax_Y = axes[2];

    multi_set_hold(axes, true);
    multi_set_aspect(axes, 1);

    ax1->title("Trajectory");
    ax1->xlabel("x");
    ax1->ylabel("y");
    ax_X->xlabel("Time (s)");
    ax_X->ylabel("x");
    ax_Y->xlabel("Time (s)");
    ax_Y->ylabel("y");

    matplot::xlim(ax1, {-0.5, 1.5});
    matplot::ylim(ax1, {-0.5, 2.5});
    matplot::xlim(ax_X, {0, T*dt});
    matplot::ylim(ax_X, {0, 1});
    matplot::xlim(ax_Y, {0, T*dt});
    matplot::ylim(ax_Y, {0, 2});

    std::vector<double> times;
    std::generate_n(
        std::back_insert_iterator<std::vector<double>>(times),
        T,
        [t = -dt, dt]() mutable { t+=dt; return t; }
    );

    matplot::line_handle lh = ax1->plot(std::vector<double>(x.begin(), x.begin()+1), std::vector<double>(y.begin(), y.begin()+1));
    ax_X->plot(times,x);
    ax_Y->plot(times,y);
    auto time_line_X = ax_X->plot(std::vector<double>{0,0}, std::vector<double>{0,1});
    auto time_line_Y = ax_Y->plot(std::vector<double>{0,0}, std::vector<double>{0,2});

    ax1->draw();
    // multi_draw(axes);


    using namespace std::chrono_literals;
    double ani_time = 5;
    double ani_dt = ani_time/T;
    for (int t=1; t < T; t++) {
        lh->x_data(std::vector<double>(x.begin(), x.begin()+t));
        lh->y_data(std::vector<double>(y.begin(), y.begin()+t));
        time_line_X->x_data(std::vector<double>{dt*t,dt*t});
        time_line_Y->x_data(std::vector<double>{dt*t,dt*t});

        ax1->draw();
        std::this_thread::sleep_for(std::chrono::duration<float>(ani_dt));
    }
    
}

