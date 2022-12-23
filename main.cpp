#include "matplotlibcpp.h"
#include "traj_opt.hpp"

namespace plt = matplotlibcpp;

int main(){
    // TrajOpt Params
    double dim_x = 3;
    double dim_u = 2;
    double T = 10;
    double dt = 1e-2;
    int N = T / dt;
    double alpha = 1e-2;

    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(dim_x, N+1);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(dim_u, N);
    Eigen::MatrixXd lam = Eigen::MatrixXd::Zero(dim_x, N+1);
    Eigen::MatrixXd s = Eigen::MatrixXd::Zero(dim_u, N);
    Eigen::VectorXd x_ref(3);
    x_ref << 3, 3, 0;

    TrajOpt traj_opt(alpha,
                     dt,
                     N, 
                     x,
                     u,
                     lam,
                     s,
                     x_ref);

    traj_opt.solve_opt_traj();

    std::vector<double> time_traj(N), x_traj(N), y_traj(N), th_traj(N), v_traj(N), omega_traj(N);
    for (int i = 0; i < N; i++){
        time_traj[i] = i*dt;
        x_traj[i] = traj_opt.m_x(0, i);
        y_traj[i] = traj_opt.m_x(1, i);
        th_traj[i] = traj_opt.m_x(2, i);
        v_traj[i] = traj_opt.m_u(0, i);
        omega_traj[i] = traj_opt.m_u(1, i);
    }
    std::cout << th_traj[N-1] << std::endl;
    plt::plot(x_traj, y_traj, {{"label", "position"}});
    plt::plot(time_traj, v_traj, {{"label", "velocity"}});
    plt::plot(time_traj, omega_traj, {{"label", "omega"}});
    plt::legend();
    plt::show();
}