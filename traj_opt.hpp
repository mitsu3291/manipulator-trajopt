#include <bits/stdc++.h>
#include <Eigen/Core>
#include "robot.hpp"

#pragma once

class TrajOpt{
public:
    Eigen::MatrixXd m_x, m_u, m_lam, m_s;
    Eigen::VectorXd m_x_ref;
    double m_alpha, m_dt;
    int m_N;
    TrajOpt(double alpha,
            double dt,
            int N,
            Eigen::MatrixXd x,
            Eigen::MatrixXd u,
            Eigen::MatrixXd lam,
            Eigen::MatrixXd s,
            Eigen::VectorXd x_ref);
    Eigen::VectorXd f(Eigen::VectorXd x, Eigen::VectorXd u);
    void solve_opt_traj();
    double calc_stage_cost(Eigen::VectorXd x, Eigen::VectorXd u);
    double calc_terminal_cost(Eigen::VectorXd x);
    Eigen::VectorXd calc_dphidx(Eigen::VectorXd x);
    Eigen::VectorXd calc_dhdx(Eigen::VectorXd x, Eigen::VectorXd u, Eigen::VectorXd lam);
    Eigen::VectorXd calc_dhdu(Eigen::VectorXd x, Eigen::VectorXd u, Eigen::VectorXd lam);
};