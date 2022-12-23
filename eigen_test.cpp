#include <Eigen/Core>
#include <bits/stdc++.h>

int main(){
    Eigen::VectorXd a(3);
    a << 1, 1, 2;
    Eigen::VectorXd b(3);
    b << 1, 0, 0;

    std::cout << (a-b).squaredNorm() << std::endl;

    return 0;
}