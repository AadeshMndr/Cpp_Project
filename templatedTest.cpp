#include<iostream>
#include"dual.h"

DualNum<double> myFunc(DualNum<double> x){
    return ((x ^ 3));
}

int main(){
    std::cout << evaluatePartialDerivative<double>(myFunc, 1.0).getExpression();

    return 0;
}