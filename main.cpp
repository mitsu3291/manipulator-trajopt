#include "matplotlibcpp.h"
#include "traj_opt.hpp"

namespace plt = matplotlibcpp;

int main(){
    // TrajOpt Params
    double dim_x = 3;
    double dim_u = 2;
    double T = 7.5;
    double dt = 1e-2;
    int N = T / dt;
    double alpha = 1e-2;

    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(dim_x, N+1);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(dim_u, N);
    Eigen::MatrixXd lam = Eigen::MatrixXd::Zero(dim_x, N+1);
    Eigen::MatrixXd s = Eigen::MatrixXd::Zero(dim_u, N);
    Eigen::VectorXd x_ref(3);
    x_ref << 3, 3, M_PI/4 ;

    // std::vector<double> time_traj(N), x_ref(N), y_ref(N), th_ref(N);
    // // for (int i = 0; i < N; i++){
    // //     time_traj[i] = i*dt;
    // //     x_ref[i] = 2*M_PI*(i*dt / T);
    // //     y_ref[i] = std::sin(x_ref[i]);
    // //     th_ref[i] = std::atan2(y_ref[i], x_ref[i]);
    // // }

    // // plt::named_plot("reference", x_ref, y_ref, "--");
    // plt::named_plot("reference", time_traj, th_ref, "--");
    // plt::legend();
    // plt::show();

    // std::cout << th_ref[0] << std::endl;
    // std::cout << th_ref[1] << std::endl;

    TrajOpt traj_opt(alpha,
                     dt,
                     N, 
                     x,
                     u,
                     lam,
                     s,
                     x_ref);

    traj_opt.solve_opt_traj();

    std::vector<double> time_traj(N), x_traj(N), y_traj(N), th_traj(N), v_traj(N), omega_traj(N), vel_th(N);
    for (int i = 0; i < N; i++){
        time_traj[i] = i*dt;
        x_traj[i] = traj_opt.m_x(0, i);
        y_traj[i] = traj_opt.m_x(1, i);
        th_traj[i] = traj_opt.m_x(2, i);
        v_traj[i] = traj_opt.m_u(0, i);
        vel_th[i] = 1.8;
        omega_traj[i] = traj_opt.m_u(1, i);
    }
    
    // plt::named_plot("position", x_traj, y_traj);
    plt::named_plot("velocity", time_traj, v_traj);
    plt::named_plot("omega", time_traj, omega_traj);
    plt::named_plot("theta", time_traj, th_traj);
    plt::named_plot("vel_th", time_traj, vel_th, "--");
    plt::legend();
    plt::show();
}