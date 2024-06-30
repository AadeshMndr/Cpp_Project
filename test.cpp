#include<iostream>
#include<iomanip>
#include"extendedDual.h"

ExtendedDualNum someFunc(ExtendedDualNum x){
    return (ExtendedDual::log(x, x/9)) / ((x ^ (-1)) / 99);
}

int main() {

    std::cout << std::setprecision(10) << ExtendedDual::partialDerivative(someFunc, 0.5, 3)[1];

    return 0;
}