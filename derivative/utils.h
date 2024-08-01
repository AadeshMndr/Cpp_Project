#pragma once
#include<iostream>
#include<string.h>
#include<chrono>
#include<thread>
#include<fstream>

using namespace std;
string ccode=R"(
#include<iostream>
#include<iomanip>
#include<limits>
#include "../Dual_Library/extendedDual.h"
using namespace ExtendedDual;
ExtendedDualNum someFunc(ExtendedDualNum);
double value_of_x();
int main() {
    std::cout << std::setprecision(10);

    int orderOfExtendedDualNumber = 9;

    vector<double> answers = ExtendedDual::partialDerivative(someFunc, value_of_x(), orderOfExtendedDualNumber);
    for (int i = 0; i < orderOfExtendedDualNumber; i++){
        std::cout << answers[i] <<std::endl;
    }

    return 0;
})";


void get_derivative(string s, double x) {
    system("rm -fv runtime_compile.cpp runtime_compile.out");
    string code=ccode+"ExtendedDualNum someFunc(ExtendedDualNum x){return "+s+";}"+"double value_of_x(){return "+to_string(x)+";}";

    // Writing 'code' to runtime_compile.cpp
    ofstream outfile("runtime_compile.cpp");
    outfile << code;
    outfile.close();

    // Compiling the code using g++
    system("g++ runtime_compile.cpp ../Dual_Library/extendedDual.cpp -o runtime_compile.out -std=c++11");

    // Executing the compiled program and redirecting output to output.txt
    system("./runtime_compile.out > output.txt 2>&1");
    // Deleting the generated cpp and out files
    system("rm runtime_compile.cpp runtime_compile.out");
    cout << "noo";
}

