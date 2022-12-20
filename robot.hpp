#include <bits/stdc++.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

class Robot{
public:
    double m_m1, m_m2, m_l1, m_l2, m_lg1, m_lg2, m_I1, m_I2, m_g, m_dt;
    double m_th1, m_th2, m_th1dot, m_th2dot, m_th1ddot, m_th2ddot;
    double x1, y1, x2, y2;

    Robot(double m1,
          double m2,
          double l1,
          double l2,
          double I1,
          double I2);

    void forward_dynamics(double tau1, double tau2);
    std::vector<double> get_joint_pos(int joint_id);
};