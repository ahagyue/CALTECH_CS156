#include<iostream>
#include "Matrix.h"

int main() {
    Matrix<double> a(3, 3, 0);
    std::cout << (a > 0);
}