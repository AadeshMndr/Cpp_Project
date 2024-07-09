#include<iostream>
#include<iomanip>
#include"extendedDual.h"
#include<limits>

ExtendedDualNum someFunc(ExtendedDualNum x){
    return (ExtendedDual::log(x, x/9)) * (x^0) / ((x ^ (-1)) / 99);
    // return ((x^20) + 5);
}

int main() {

    std::cout << std::setprecision(10) << ExtendedDual::partialDerivative(someFunc, 0.5, 6)[5];

    // double h = std::numeric_limits<double>::min();
    double h = 0.000000000000001;

    ExtendedDualNum x(0.5, 3);

    // std::cout << "\n->" << h;

    std::cout << "\n " << ((someFunc(x + h) - someFunc(x)) / h).getNum();

    return 0;
}