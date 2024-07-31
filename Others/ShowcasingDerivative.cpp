#include<iostream>
#include<iomanip>
#include"../Dual_Library/extendedDual.h"
#include<limits>

ExtendedDualNum someFunc(ExtendedDualNum x){
    return (ExtendedDual::log(x, x/9)) * (x^0) / ((x ^ (-1)) / 99);
}

int main() {

    std::cout << std::setprecision(10);

    int orderOfExtendedDualNumber = 9; //upto which derivative

    vector<double> answers = ExtendedDual::partialDerivative(someFunc, 0.5, orderOfExtendedDualNumber);

    for (int i = 0; i < orderOfExtendedDualNumber; i++){
        std::cout << "\nThe " << i << " order derivative is: " << answers[i] << std::endl;
    }

    // double h = std::numeric_limits<double>::min();
    double h = 0.000000000000001;

    ExtendedDualNum x(0.5, 3);

    // std::cout << "\n->" << h << std::endl;

    std::cout << "\n\nThe derivative by numerical method is: " << ((someFunc(x + h) - someFunc(x)) / h).getNum() << std::endl;
    std::cout << "\nWhile the one obtained by automatic differentiation is: " << answers[1] << std::endl;

    return 0;
}