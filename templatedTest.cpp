#include<iostream>
#include"dual.h"

DualNum<DualNum<double>> myFunc(DualNum<DualNum<double>> x){
    return ((x ^ DualNum(3.0)));
}

int main(){
    std::cout << evaluatePartialDerivative<DualNum<double>>(myFunc, DualNum<double>(1.0, 1.0)).getDual().getExpression();

    return 0;
}