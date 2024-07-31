#include<iostream>
#include "../Dual_Library/extendedDual.h"

ExtendedDualNum potential(vector<ExtendedDualNum> params){
    return (((params[0])^2) + params[1] + params[2]);
}

ExtendedDualNum fx(vector<ExtendedDualNum> params){
    return (2 * params[0] + params[1]);
}

ExtendedDualNum fy(vector<ExtendedDualNum> params){
    return (3 * params[0] + params[1]);
}

ExtendedDualNum someFunction(ExtendedDualNum x){
    return ((x ^ 2) + 5 * x - 50);
}

int main(){
    vector<double> gradientValue = ExtendedDual::gradient(potential, {1, 1, 1});

    std::cout << "\n\nf(x, y, z) = x^2 + y + z   at  (1, 1, 1)" << std::endl;
    std::cout << "\nThe gradient is: " << "(" << gradientValue[0] << ", " << gradientValue[1] << ", " << gradientValue[2] << ")" << std::endl;

    std::cout << "\nThe laplacian is: " << ExtendedDual::laplacian(potential, {ExtendedDualNum(1, 3), ExtendedDualNum(1, 3), ExtendedDualNum(1, 3)}) << std::endl;

    vector<vector<double>> jacobianMatrix = ExtendedDual::jacobian({fx, fy}, {1, 1});

    std::cout << "\n\nf(x, y) = (2x + y, 3x + y)   at  (1, 1)" << std::endl;
    std::cout << "\nThe Jacobian Matrix is: \n\n";
    for (int i = 0; i < jacobianMatrix.size(); i++){
        for (int j = 0; j < jacobianMatrix[0].size(); j++){
            std::cout << jacobianMatrix[i][j] << "\t";
        }
        std::cout << "\n\n";
    }

    double answer = ExtendedDual::solveUsingNewtonRaphson(someFunction);


    std::cout << "\nequation: x^2 + 5x - 50 = 0" << std::endl;
    std::cout << "\nThus the solution of the equation is " << answer << "\n\n" << std::endl;
    
    return 0;
}   