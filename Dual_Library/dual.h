#ifndef DUALNUMBER
#define DUALNUMBER

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include <thread>

using std::vector;
using std::string;

class DualNum {
private:
    long double real;
    long double dual;

public:
    DualNum(long double r = 0, long double e = 0);
    DualNum(const DualNum& num);

    long double getReal() const;
    long double getDual() const;

    void setReal(const long double& r);
    void setDual(const long double& e);

    string getExpression() const;

    void operator = (const long double& realVal);
    DualNum operator * (const DualNum& num) const;
    void operator *= (const DualNum& num);
    DualNum operator * (const long double& k) const;
    void operator *= (const long double& k);
    DualNum operator / (const DualNum& num) const;
    void operator /= (const DualNum& num);
    DualNum operator / (const long double& k) const;
    void operator /= (const long double& k);
    DualNum operator + (const DualNum& num) const;
    void operator += (const DualNum& num);
    DualNum operator + (const long double& k) const;
    void operator += (const long double& k);
    DualNum operator - (const DualNum& num) const;
    void operator -= (const DualNum& num);
    DualNum operator - (const long double& k) const;
    void operator -= (const long double& k);
    DualNum operator ^ (const DualNum& num) const;
    void operator ^= (const DualNum& num);
    DualNum operator ^ (const long double& k) const;
    void operator ^= (const long double& k);
};

// Normal Numbers Operator Overloading
DualNum operator * (const long double& k, const DualNum& num);
DualNum operator / (const long double& k, const DualNum& num);
DualNum operator + (const long double& k, const DualNum& num);
DualNum operator - (const long double& k, const DualNum& num);
DualNum operator ^ (const long double& k, const DualNum& num);

// Some standard functions
namespace Dual {
    DualNum pow(const DualNum& x, const DualNum& y);
    DualNum exp(const DualNum& x);
    DualNum log(const DualNum& x);
    DualNum relu(const DualNum& x);
    DualNum tanh(const DualNum& x);
    DualNum sigmoid(const DualNum& x);
    DualNum softmax(const vector<DualNum>& X, int index, DualNum sum = DualNum(0, 0));
    DualNum mse(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);
    DualNum binary_crossentropy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);
    DualNum categorical_crossentropy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);
    DualNum accuracy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat, const double threshold = 0.5);
    vector<vector<DualNum>> originalMatMul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B);
    void partialMatMul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B, vector<vector<DualNum>>& result, int start, int end);
    vector<vector<DualNum>> matmul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B);
    void displayMatrix(const vector<vector<DualNum>>& mat);
}

// Some extra functions
long double partialDerivative(DualNum(*func)(vector<DualNum>), vector<DualNum> params, int paramIndex = 0, long double at = 1);

template<typename Func>
long double partialDerivative(Func func, vector<DualNum> params, int paramIndex = 0, long double at = 1) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return (func(params)).getDual();
}

long double partialDerivative(DualNum(*func)(DualNum), long double at = 1);
DualNum evaluatePartialDerivative(DualNum(*func)(DualNum), long double at = 1);

enum class coordinate_system { cartesian, cylindrical, spherical };

vector<DualNum> gradient(DualNum(*func)(vector<DualNum>), vector<DualNum> at, coordinate_system system = coordinate_system::cartesian);

double solveUsingNewtonRaphson(DualNum(*func)(DualNum), double initialGuess = 1, int max_no_of_iterations = 10000);

#endif