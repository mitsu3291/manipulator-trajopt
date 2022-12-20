#include "robot.hpp"

Robot::Robot(double m1,
             double m2,
             double l1,
             double l2,
             double I1,
             double I2){
            m_m1 = m1;
            m_m2 = m2;
            m_l1 = l1; 
            m_l2 = l2;
            m_lg1 = m_l1 / 2;
            m_lg2 = m_l2 / 2;
            m_I1 = I1;
            m_I2 = I2;

            m_th1 = 0;
            m_th2 = 0;
            m_th1dot = 0;
            m_th2dot = 0;
            m_th1ddot = 0;
            m_th2ddot = 0;

            m_g = 0.98;
            m_dt = 0.001;
}

void Robot::forward_dynamics(double tau1, double tau2){
    Eigen::MatrixXd M(2,2);
    Eigen::MatrixXd h(2,1);
    Eigen::MatrixXd g(2,1);
    Eigen::MatrixXd tau(2,1);
    Eigen::MatrixXd thddot(2,1);

    M(0,0) = m_m1*std::pow(m_lg1,2) + m_I1 + m_m2*(std::pow(m_l1,2) + std::pow(m_lg2,2) + 2*m_l1*m_lg2*std::cos(m_th2)) + m_I2;
    M(0,1) = m_m2*(std::pow(m_lg2,2) + m_l1*m_lg2*std::cos(m_th2)) + m_I2;
    M(1,0) = m_m2*(std::pow(m_lg2,2) + m_l1*m_lg2*std::cos(m_th2)) + m_I2;
    M(1,1) = m_m2*std::pow(m_lg2,2) + m_I2;

    h(0,0) = -m_m2*m_l1*m_lg2*std::sin(m_th2)*(2*m_th1dot*m_th2dot+std::pow(m_th2dot, 2));
    h(1,0) = m_m2*m_l1*m_lg2*std::sin(m_th2)*std::pow(m_th1dot, 2);

    g(0,0) = m_m1*m_g*m_lg1*std::cos(m_th1) + m_m2*m_g*(m_l1*std::cos(m_th1) + m_lg2*std::cos(m_th1 + m_th2));
    g(1,0) = m_m2*m_g*m_lg2*std::cos(m_th1 + m_th2);

    tau(0,0) = tau1;
    tau(1,0) = tau2;

    thddot = M.inverse()*(tau - h - g);

    // Update Joint Angle Params
    m_th1ddot = thddot(0,0);
    m_th2ddot = thddot(1,0);

    m_th1dot += m_dt*m_th1ddot;
    m_th2dot += m_dt*m_th2ddot;

    m_th1 += m_dt*m_th1dot;
    m_th2 += m_dt*m_th2dot;
}

std::vector<double> Robot::get_joint_pos(int joint_id){
    double x, y;
    if (joint_id == 1){
        x = m_l1*std::cos(m_th1);
        y = m_l1*std::sin(m_th1);
    }else if (joint_id == 2){
        x = m_l1*std::cos(m_th1) + m_l2*std::cos(m_th1 + m_th2);
        y = m_l1*std::sin(m_th1) + m_l2*std::sin(m_th1 + m_th2);
    }else{
        std::cout << "Joint Id " << joint_id << " is not defined" << std::endl;
        std::exit(0);
    }

    std::vector<double> joint_pos = {x, y};

    return joint_pos;
}