#include "robot.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(){
    double m1 = 1;
    double m2 = 1;
    double l1 = 1;
    double l2 = 1;
    double I1 = 1;
    double I2 = 1;

    Robot robot(m1,
                m2,
                l1,
                l2,
                I1,
                I2);

    int n = 50000;
    std::vector<double> origin = {0,0};
    std::vector<double> joint_pos1, joint_pos2;
    std::vector<double> x1_traj(n), y1_traj(n), x2_traj(n), y2_traj(n);
    for (int i = 0; i < n; i++){
        robot.forward_dynamics(0,0);
        joint_pos1 = robot.get_joint_pos(1);
        joint_pos2 = robot.get_joint_pos(2);
        if (i % 50 == 0){
            // figure settings
            plt::clf();
            plt::xlim(-2.5, 2.5);
            plt::ylim(-2.5, 2.5);

            // plot
            plt::plot({origin[0], joint_pos1[0]}, {origin[1], joint_pos1[1]}, "bo-");
            plt::plot({joint_pos1[0], joint_pos2[0]}, {joint_pos1[1], joint_pos2[1]}, "bo-");
            plt::pause(0.01);
        }
    }
    std::cout << robot.m_th1 << " " << robot.m_th2 << std::endl;
    // plt::show();
}