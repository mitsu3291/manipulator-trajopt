#include "traj_opt.hpp"

TrajOpt::TrajOpt(double alpha,
                 double dt,
                 int N,
                 Eigen::MatrixXd x,
                 Eigen::MatrixXd u,
                 Eigen::MatrixXd lam,
                 Eigen::MatrixXd s,
                 Eigen::VectorXd x_ref){
    m_alpha = alpha;
    m_dt = dt;
    m_N = N;
    m_x = x;
    m_u = u;
    m_lam = lam;
    m_s = s;
    m_x_ref = x_ref;
}

Eigen::VectorXd TrajOpt::f(Eigen::VectorXd x, Eigen::VectorXd u){
    // control info
    Eigen::VectorXd f(3);
    f[0] = u[0]*std::cos(x[2]);
    f[1] = u[0]*std::sin(x[2]);
    f[2] = u[1];

    return f;
}

double TrajOpt::calc_stage_cost(Eigen::VectorXd x, Eigen::VectorXd u){
    double L = (1/2)*((x-m_x_ref).squaredNorm() + u.squaredNorm());
    if (u[0] > 1.8){
        L += 1000*std::pow((u[0] - 1.8), 2);
    }

    return L;
}

double TrajOpt::calc_terminal_cost(Eigen::VectorXd x){
    double phi = (1/2)*(x-m_x_ref).squaredNorm();

    return phi;
}

Eigen::VectorXd TrajOpt::calc_dphidx(Eigen::VectorXd x){
    Eigen::VectorXd dphidx = x - m_x_ref;

    return dphidx;
}

Eigen::VectorXd TrajOpt::calc_dhdx(Eigen::VectorXd x, Eigen::VectorXd u, Eigen::VectorXd lam){
    Eigen::VectorXd dhdx(3);
    dhdx[0] = x[0] - m_x_ref[0];
    dhdx[1] = x[1] - m_x_ref[1];
    dhdx[2] = x[2] - m_x_ref[2] - lam[0]*u[0]*std::sin(x[2]) + lam[1]*u[0]*std::cos(x[2]);

    return dhdx;
}

Eigen::VectorXd TrajOpt::calc_dhdu(Eigen::VectorXd x, Eigen::VectorXd u, Eigen::VectorXd lam){
    Eigen::VectorXd dhdu(2);
    dhdu[0] = u[0] + lam[0]*std::cos(x[2]) + lam[1]*std::sin(x[2]);
    dhdu[1] = u[1] + lam[2];
    if (u[0] > 1.8){
        dhdu[0] += 200*(u[0] - 1.8);
    }

    return dhdu;
}

void TrajOpt::solve_opt_traj(){
    for (int i = 0; i < 10000; i++){
        // advance simulation and get x
        for (int j = 0; j < m_N; j++){
            m_x.col(j+1) = m_x.col(j) + this->f(m_x.col(j), m_u.col(j))*m_dt;
        }

        // caluculate lambda
        m_lam.col(m_N) = this->calc_dphidx(m_x.col(m_N));
        for (int j = m_N-1; j > 0; j--){
            m_lam.col(j) = m_lam.col(j+1) - this->calc_dhdx(m_x.col(j), m_u.col(j), m_lam.col(j+1))*(-m_dt);
        }

        // Update control input or end iteration
        for (int j = 0; j < m_N; j++){
            m_s.col(j) = - this->calc_dhdu(m_x.col(j), m_u.col(j), m_lam.col(j+1));
        }
        double sol_err = std::pow(m_s.squaredNorm(), 0.5);
        if (sol_err < 1e-2){
            std::cout << "converge" << std::endl;
            break;
        }

        for (int j = 0; j < m_N; j++){
            m_u.col(j) += m_alpha*m_s.col(j); 
        }
    }
    std::cout << "error did not converge" << std::endl;
}