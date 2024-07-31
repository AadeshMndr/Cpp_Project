#ifndef EXTENDED_DUAL 
#define EXTENDED_DUAL

#include<vector>
#include<iostream>
#include<stdexcept>

using std::vector;

class ExtendedDualNum;

namespace ExtendedDual {
    //General
    long int factorial(int);
    ExtendedDualNum exp(const ExtendedDualNum& num);
    ExtendedDualNum log(const ExtendedDualNum& num);
    ExtendedDualNum log(const ExtendedDualNum& num, const ExtendedDualNum& base);
    ExtendedDualNum log(const ExtendedDualNum& num, const double base);
    ExtendedDualNum sin(const ExtendedDualNum& num);
    ExtendedDualNum cos(const ExtendedDualNum& num);
    ExtendedDualNum tan(const ExtendedDualNum& num);
    ExtendedDualNum sinh(const ExtendedDualNum& num);
    ExtendedDualNum cosh(const ExtendedDualNum& num);
    ExtendedDualNum tanh(const ExtendedDualNum& num);
    ExtendedDualNum asin(const ExtendedDualNum& num);
    ExtendedDualNum acos(const ExtendedDualNum& num);
    ExtendedDualNum atan(const ExtendedDualNum& num);

    //Utility
    ExtendedDualNum derivativeOfASin(ExtendedDualNum num);
    ExtendedDualNum derivativeOfATan(ExtendedDualNum num);

    //Matrix related
    vector<vector<double>> matrixFiller(vector<double> elems);
    vector<vector<double>> minorMatrix(vector<vector<double>> mat, int row, int col);
    double determinant(vector<vector<double>> mat);
    double originalDeterminant(vector<vector<double>> mat);
    vector<vector<double>> originalInverseMatrix(vector<vector<double>> mat);
    vector<vector<double>> inverseMatrix(vector<vector<double>> mat);
    vector<vector<double>> originalMatMul(vector<vector<double>> A, vector<vector<double>> B);
    vector<vector<double>> matMul(vector<vector<double>> A, vector<vector<double>> B);

    //Extra
    vector<double> partialDerivative(ExtendedDualNum (*func)(vector<ExtendedDualNum>), vector<ExtendedDualNum> params, int paramIndex = 0);
    vector<double> partialDerivative(ExtendedDualNum (*func)(ExtendedDualNum), double at = 1, int order = 3);
    ExtendedDualNum evaluatePartialDerivative(ExtendedDualNum(*func)(ExtendedDualNum), long double at = 1);


    enum class coordinate_system { cartesian, cylindrical, spherical };
    using VectorFunctionPointer = ExtendedDualNum(*)(vector<ExtendedDualNum>);
    vector<double> gradient(VectorFunctionPointer, vector<ExtendedDualNum> at, coordinate_system system = coordinate_system::cartesian);
    double laplacian(VectorFunctionPointer, vector<ExtendedDualNum> at, coordinate_system system = coordinate_system::cartesian);
    vector<vector<double>> jacobian(vector<ExtendedDual::VectorFunctionPointer>, vector<ExtendedDualNum> at);
    double solveUsingNewtonRaphson(ExtendedDualNum(*func)(ExtendedDualNum), double initialGuess = 1, int max_no_of_iterations = 10000);
};

class ExtendedDualNum {
private:
    vector<vector<double>> data;
    int order = 2;

public:
    ExtendedDualNum(double num, int orderOfMat = 2);

    ExtendedDualNum(vector<vector<double>> dataMatrix);

    ExtendedDualNum(const ExtendedDualNum& another);

    vector<vector<double>> getData() const;

    double getNum() const;

    int getSize() const;

    vector<vector<double>> getDualMatrix() const;

    void setData(vector<vector<double>> newData);

    void setDual(double dual = 1);

    void displayData() const;

    ExtendedDualNum operator * (const ExtendedDualNum& another) const;
    void operator *= (const ExtendedDualNum& another);
    ExtendedDualNum operator * (const double num) const;
    void operator *= (const double num);
    friend ExtendedDualNum operator * (const double num, const ExtendedDualNum& another);

    ExtendedDualNum operator + (const ExtendedDualNum& another) const;
    void operator += (const ExtendedDualNum& another);
    ExtendedDualNum operator + (const double num) const;
    void operator += (const double num);
    friend ExtendedDualNum operator + (const double num, const ExtendedDualNum& another);

    ExtendedDualNum operator - (const ExtendedDualNum& another) const;
    void operator -= (const ExtendedDualNum& another);
    ExtendedDualNum operator - (const double num) const;
    void operator -= (const double num);
    friend ExtendedDualNum operator - (const double num, const ExtendedDualNum& another);

    ExtendedDualNum operator - ();
    ExtendedDualNum operator + ();

    ExtendedDualNum operator / (const ExtendedDualNum& another) const;
    void operator /= (const ExtendedDualNum& another);
    ExtendedDualNum operator / (const double num) const;
    void operator /= (const double num);
    friend ExtendedDualNum operator / (const double num, const ExtendedDualNum& another);

    ExtendedDualNum operator ^ (const int n) const;
    void operator ^= (const int n);
    ExtendedDualNum operator ^ (const double n) const;
    void operator ^= (const double n);
    ExtendedDualNum operator ^ (const ExtendedDualNum& another) const;
    void operator ^= (const ExtendedDualNum& another);
    friend ExtendedDualNum operator ^ (const double n, const ExtendedDualNum& num); //This is not actually required to be a friend 

    // operator double();

    bool operator == (const ExtendedDualNum& num) const;
    bool operator > (const ExtendedDualNum& num) const;
    bool operator < (const ExtendedDualNum& num) const;
    bool operator >= (const ExtendedDualNum& num) const;
    bool operator <= (const ExtendedDualNum& num) const;

    double getPartialDerivative(int i = 1) const;
};

std::ostream& operator << (std::ostream &os, const ExtendedDualNum& num);

#endif