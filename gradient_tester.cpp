#include<iostream>
#include"rough.h"

DualNum potential(vector<DualNum> params){
    return (((params[0])^2) + params[1] + params[2]);
}

int main(){
    vector<DualNum> gradientValue = gradient(potential, {1, 1, 1}, coordinate_system::cartesian);

    std::cout << "(" << gradientValue[0].getReal() << ", " << gradientValue[1].getReal() << ", " << gradientValue[2].getReal() << " )" << std::endl;

    return 0;
}   