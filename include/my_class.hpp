// my_class.hpp

#ifndef MY_CLASS_HPP
#define MY_CLASS_HPP

#include "./my_base.hpp"  // Include the virtual parent class header

class MyClass : public MyBase {
public:
    MyClass();          // Constructor declaration
    void doSomething();  // Member function declaration

    // Override the virtual function from MyBase
    void someVirtualFunction(std::vector<Eigen::VectorXd> x) override;
private:
    int data;
};

#endif  // MY_CLASS_HPP
