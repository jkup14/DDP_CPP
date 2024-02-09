// my_base.hpp

#ifndef MY_BASE_HPP
#define MY_BASE_HPP

#include <vector>
#include <Eigen/Cholesky>

class MyBase {
public:
    virtual ~MyBase() = default;  // Virtual destructor
    virtual void someVirtualFunction(std::vector<Eigen::VectorXd> x) = 0;  // Pure virtual function
};

#endif  // MY_BASE_HPP