#include<iostream>
#include<vector>
#include<stdexcept>
#include<cmath>
#include"extendedDual.h"

using std::vector;

long int ExtendedDual::factorial(int n) {
    if (n < 0) {
        throw std::runtime_error("The factorial of a Negative Number is trying to be taken !");
    }
    else if (n == 0) {
        return 1;
    }
    else if (n <= 2) {
        return n;
    }
    else if (n == 3) {
        return 6;
    }
    else if (n == 4) {
        return 24;
    }
    else if (n == 5) {
        return 120;
    }
    else {
        long int fact = 1;
        for (int i = 2; i <= n; i++) {
            fact *= i;
        }

        return fact;
    }
}

ExtendedDualNum ExtendedDual::exp(const ExtendedDualNum& num) {
    int order = num.getSize();

    ExtendedDualNum result(1, order);

    ExtendedDualNum dualMat = num.getDualMatrix();

    //This is correct...(2)
    for (int i = 1; i < order; i++) {
        result += ((dualMat ^ i) / ExtendedDual::factorial(i));
    }

    return (result * std::exp(num.getNum()));
}

ExtendedDualNum ExtendedDual::log(const ExtendedDualNum& num) {
    int order = num.getSize();

    double real = num.getNum();

    ExtendedDualNum result(std::log(real), order);

    ExtendedDualNum dualMat = num.getDualMatrix();

    for (int i = 1; i < order; i++) {
        result += (std::pow(-1, i + 1) * ((dualMat / real) ^ i) / i);
    }

    return result;
}

ExtendedDualNum ExtendedDual::log(const ExtendedDualNum& num, const ExtendedDualNum& base) {
    return (ExtendedDual::log(num) / ExtendedDual::log(base));
}

ExtendedDualNum ExtendedDual::log(const ExtendedDualNum& num, const double base) {
    return (ExtendedDual::log(num) / std::log(base));
}

ExtendedDualNum ExtendedDual::sin(const ExtendedDualNum& num) {
    int order = num.getSize();

    ExtendedDualNum sinMat(std::sin(num.getNum()), order);
    ExtendedDualNum cosMat(std::cos(num.getNum()), order);

    ExtendedDualNum dualMat = num.getDualMatrix();

    ExtendedDualNum result = sinMat;

    for (int i = 1; i < num.getSize(); i++) {
        ExtendedDualNum intermediate = (dualMat ^ i) / ExtendedDual::factorial(i);

        if (i % 2 == 0) {
            result += std::pow(-1, (i / 2)) * sinMat * intermediate;
        }
        else {
            result += std::pow(-1, ((i - 1) / 2)) * cosMat * intermediate;
        }
    }

    return result;
}

ExtendedDualNum ExtendedDual::cos(const ExtendedDualNum& num) {
    int order = num.getSize();

    ExtendedDualNum sinMat(std::sin(num.getNum()), order);
    ExtendedDualNum cosMat(std::cos(num.getNum()), order);

    ExtendedDualNum dualMat = num.getDualMatrix();

    ExtendedDualNum result = cosMat;

    for (int i = 1; i < num.getSize(); i++) {
        ExtendedDualNum intermediate = (dualMat ^ i) / ExtendedDual::factorial(i);

        if (i % 2 == 0) {
            result += std::pow(-1, (i / 2)) * cosMat * intermediate;
        }
        else {
            result += std::pow(-1, ((i - 1) / 2) + 1) * sinMat * intermediate;
        }
    }

    return result;
}

ExtendedDualNum ExtendedDual::tan(const ExtendedDualNum& num) {
    return (ExtendedDual::sin(num) / ExtendedDual::cos(num));
}

ExtendedDualNum ExtendedDual::sinh(const ExtendedDualNum& num) {
    return ((ExtendedDual::exp(num) - ExtendedDual::exp(-1 * num)) / 2);
}

ExtendedDualNum ExtendedDual::cosh(const ExtendedDualNum& num) {
    return ((ExtendedDual::exp(num) + ExtendedDual::exp(-1 * num)) / 2);
}

ExtendedDualNum ExtendedDual::tanh(const ExtendedDualNum& num) {
    return ((ExtendedDual::exp(num) - ExtendedDual::exp(-1 * num)) / (ExtendedDual::exp(num) + ExtendedDual::exp(-1 * num)));
}

ExtendedDualNum ExtendedDual::asin(const ExtendedDualNum& num){
    
    int order = num.getSize();

    double real = num.getNum();

    ExtendedDualNum dualMat = num.getDualMatrix();

    vector<double> derivatives = ExtendedDual::partialDerivative(derivativeOfASin, real, order-1);

    ExtendedDualNum result(std::asin(real), order);

    for (int i = 1; i < order; i++){
        result += derivatives[i-1] * (dualMat ^ i) / ExtendedDual::factorial(i);
    }

    return result;
}

ExtendedDualNum ExtendedDual::acos(const ExtendedDualNum& num){
    
    int order = num.getSize();

    double real = num.getNum();

    ExtendedDualNum dualMat = num.getDualMatrix();

    vector<double> derivatives = ExtendedDual::partialDerivative(derivativeOfASin, real, order-1);

    ExtendedDualNum result(std::acos(real), order);

    for (int i = 1; i < order; i++){
        result += (-1) * derivatives[i-1] * (dualMat ^ i) / ExtendedDual::factorial(i);
    }

    return result;
}

ExtendedDualNum ExtendedDual::atan(const ExtendedDualNum& num){
    
    int order = num.getSize();

    double real = num.getNum();

    ExtendedDualNum dualMat = num.getDualMatrix();

    vector<double> derivatives = ExtendedDual::partialDerivative(derivativeOfATan, real, order-1);

    ExtendedDualNum result(std::atan(real), order);

    for (int i = 1; i < order; i++){
        result += derivatives[i-1] * (dualMat ^ i) / ExtendedDual::factorial(i);
    }

    return result;
}

ExtendedDualNum ExtendedDual::derivativeOfASin(ExtendedDualNum num){
    return (1 / ((1 - (num ^ 2)) ^ (0.5)));
}

ExtendedDualNum ExtendedDual::derivativeOfATan(ExtendedDualNum num){
    return (1 / (1 + (num ^ 2)));
}

vector<vector<double>> ExtendedDual::matrixFiller(vector<double> elems) {
    int rows = elems.size();

    vector<vector<double>> result(rows, vector<double>(rows, 0));

    for (int i = 0; i < rows; i++) {
        int count = 0;
        for (int j = i; j < rows; j++) {
            result[i][j] = elems[count];
            count++;
        }
    }

    return result;
}

vector<vector<double>> ExtendedDual::minorMatrix(vector<vector<double>> mat, int row, int col) {
    vector<vector<double>> minor;

    for (int i = 0; i < mat.size(); i++) {
        if (i == row) {
            continue;
        }

        vector<double> eachRow;

        for (int j = 0; j < mat[0].size(); j++) {
            if (j == col) {
                continue;
            }

            eachRow.push_back(mat[i][j]);
        }

        minor.push_back(eachRow);
    }

    return minor;
}

double ExtendedDual::determinant(vector<vector<double>> mat) {

    return std::pow(mat[0][0], mat.size());
}

double ExtendedDual::originalDeterminant(vector<vector<double>> mat) {
    int rows = mat.size();

    int cols = mat[0].size();

    if (rows != cols) {
        std::cout << "\nThe no of rows != no of cols !!" << std::endl;
        return 0;
    }

    if (rows == 2) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }
    else if (rows == 1) {
        return mat[0][0];
    }
    else {
        double result = 0;

        for (int i = 0; i < cols; i++) {
            result += pow(-1, i) * originalDeterminant(minorMatrix(mat, 0, i)) * mat[0][i];
        }


        return result;
    }
}

vector<vector<double>> ExtendedDual::originalInverseMatrix(vector<vector<double>> mat) {
    int rows = mat.size();
    int cols = mat[0].size();

    if (rows != cols) {
        throw std::runtime_error("\nThe no of rows != no of cols !!");
    }

    vector<vector<double>> cofTransMat(rows, vector<double>(cols, 0));

    double det = determinant(mat);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cofTransMat[j][i] = pow(-1, i + j) * originalDeterminant(minorMatrix(mat, i, j)) / det;
        }
    }

    return cofTransMat;
}

vector<vector<double>> ExtendedDual::inverseMatrix(vector<vector<double>> mat) {
    int rows = mat.size();

    vector<double> elems;

    double det = determinant(mat);

    elems.push_back((1 / mat[0][0]));

    for (int i = 1; i < rows; i++) {
        elems.push_back(pow(-1, i) * originalDeterminant(minorMatrix(mat, i, 0)) / det);
    }

    return matrixFiller(elems);
}

vector<vector<double>> ExtendedDual::originalMatMul(vector<vector<double>> A, vector<vector<double>> B) {
    int Arows = A.size();
    int Acols = A[0].size();
    int Brows = B.size();
    int Bcols = B[0].size();

    if (Acols != Brows) {
        throw std::runtime_error("The shapes of the Matrices don't match !");
    }

    vector<vector<double>> result(Arows, vector<double>(Bcols, 0));

    for (int i = 0; i < Arows; i++) {
        for (int j = 0; j < Bcols; j++) {
            for (int k = 0; k < Brows; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

vector<vector<double>> ExtendedDual::matMul(vector<vector<double>> A, vector<vector<double>> B) {
    int rows = A.size();

    vector<double> elems(rows, 0);

    for (int j = 0; j < rows; j++) {
        for (int k = 0; k < rows; k++) {
            elems[j] += A[0][k] * B[k][j];
        }
    }

    return matrixFiller(elems);
}


//default int paramIndex = 0;
vector<double> ExtendedDual::partialDerivative(ExtendedDualNum(*func)(vector<ExtendedDualNum>), vector<ExtendedDualNum> params, int paramIndex) {
    params[paramIndex].setDual(1);

    vector<double> result = func(params).getData()[0];

    for (int i = 0; i < result.size(); i++) {
        result[i] *= ExtendedDual::factorial(i);
    }

    return result;
}

//default double at = 1 and int order = 3;
vector<double> ExtendedDual::partialDerivative(ExtendedDualNum(*func)(ExtendedDualNum), double at, int order) {
    ExtendedDualNum num(at, order);
    num.setDual(1);

    vector<double> result = func(num).getData()[0];

    for (int i = 0; i < result.size(); i++) {
        result[i] *= ExtendedDual::factorial(i);
    }

    return result;
}




















//Class related definitions only..........................................................


//orderOfMat is a default argument with value = 2, it is omitted here because it is declared in the header file
ExtendedDualNum::ExtendedDualNum(double num, int orderOfMat) {
    data = vector<vector<double>>(orderOfMat, vector<double>(orderOfMat, 0));
    for (int i = 0; i < orderOfMat; i++) {
        data[i][i] = num;
    }

    order = orderOfMat;
}

ExtendedDualNum::ExtendedDualNum(vector<vector<double>> dataMatrix) {
    data = dataMatrix;
    order = dataMatrix.size();
}

ExtendedDualNum::ExtendedDualNum(const ExtendedDualNum& another) {
    data = another.data;
    order = another.order;
}

vector<vector<double>> ExtendedDualNum::getData() const {
    return data;
}

double ExtendedDualNum::getNum() const {
    return data[0][0];
}

int ExtendedDualNum::getSize() const {
    return order;
}

vector<vector<double>> ExtendedDualNum::getDualMatrix() const {
    vector<vector<double>> temp = data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] = 0;
    }

    return temp;
}

void ExtendedDualNum::setData(vector<vector<double>> newData) {
    data = newData;
}


//default double dual = 1;
void ExtendedDualNum::setDual(double dual) {
    for (int i = 1; i < order; i++) {
        data[i - 1][i] = dual;
    }
}

void ExtendedDualNum::displayData() const {

    for (int i = 0; i < order; i++) {
        std::cout << "\n\n";
        for (int j = 0; j < order; j++) {
            std::cout << data[i][j] << "\t";
        }
    }
}

ExtendedDualNum ExtendedDualNum::operator * (const ExtendedDualNum& another) const {
    return ExtendedDual::matMul(data, another.data);
}

void ExtendedDualNum::operator *= (const ExtendedDualNum& another) {
    data = ExtendedDual::matMul(data, another.data);
}

ExtendedDualNum ExtendedDualNum::operator * (const double num) const {

    vector<double> elems = data[0];

    for (int i = 0; i < elems.size(); i++) {
        elems[i] *= num;
    }

    return ExtendedDual::matrixFiller(elems);
}

void ExtendedDualNum::operator *= (const double num) {

    vector<double> elems = data[0];

    for (int i = 0; i < elems.size(); i++) {
        elems[i] *= num;
    }

    data = ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum operator * (const double num, const ExtendedDualNum& another) {

    vector<double> elems = another.data[0];

    for (int i = 0; i < elems.size(); i++) {
        elems[i] *= num;
    }

    return ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum ExtendedDualNum::operator + (const ExtendedDualNum& another) const {

    vector<double> elems(order, 0);

    for (int i = 0; i < order; i++) {
        elems[i] = data[0][i] + another.data[0][i];
    }

    return ExtendedDual::matrixFiller(elems);
}

void ExtendedDualNum::operator += (const ExtendedDualNum& another) {

    vector<double> elems(order, 0);

    for (int i = 0; i < order; i++) {
        elems[i] = data[0][i] + another.data[0][i];
    }

    data = ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum ExtendedDualNum::operator + (const double num) const {

    vector<vector<double>> temp = data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] += num;
    }

    return temp;
}

void ExtendedDualNum::operator += (const double num) {

    vector<vector<double>> temp = data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] += num;
    }

    data = temp;
}

ExtendedDualNum operator + (const double num, const ExtendedDualNum& another) {

    vector<vector<double>> temp = another.data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] += num;
    }

    return temp;
}

ExtendedDualNum ExtendedDualNum::operator - (const ExtendedDualNum& another) const {

    vector<double> elems(order, 0);

    for (int i = 0; i < order; i++) {
        elems[i] = data[0][i] - another.data[0][i];
    }

    return ExtendedDual::matrixFiller(elems);
}

void ExtendedDualNum::operator -= (const ExtendedDualNum& another) {

    vector<double> elems(order, 0);

    for (int i = 0; i < order; i++) {
        elems[i] = data[0][i] - another.data[0][i];
    }

    data = ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum ExtendedDualNum::operator - (const double num) const {

    vector<vector<double>> temp = data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] -= num;
    }

    return temp;
}

void ExtendedDualNum::operator -= (const double num) {

    vector<vector<double>> temp = data;

    for (int i = 0; i < temp.size(); i++) {
        temp[i][i] -= num;
    }

    data = temp;
}

ExtendedDualNum operator - (const double num, const ExtendedDualNum& another) {

    return (ExtendedDualNum(num, another.getSize()) -  another);
}

ExtendedDualNum ExtendedDualNum::operator - (){
    vector<double> elems(order, 0);

    for (int i = 0; i < order; i++) {
        elems[i] = (-1) * data[0][i];
    }

    return ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum ExtendedDualNum::operator + (){
    return *this;
}

ExtendedDualNum ExtendedDualNum::operator / (const ExtendedDualNum& another) const {
    vector<vector<double>> inverse = ExtendedDual::inverseMatrix(another.data);

    return ExtendedDual::matMul(data, inverse);
}

void ExtendedDualNum::operator /= (const ExtendedDualNum& another) {
    vector<vector<double>> inverse = ExtendedDual::inverseMatrix(another.data);

    data = ExtendedDual::matMul(data, inverse);
}

ExtendedDualNum ExtendedDualNum::operator / (const double num) const {

    vector<double> elems = data[0];

    for (int i = 0; i < elems.size(); i++) {
        elems[i] /= num;
    }

    return ExtendedDual::matrixFiller(elems);
}

void ExtendedDualNum::operator /= (const double num) {

    vector<double> elems = data[0];

    for (int i = 0; i < elems.size(); i++) {
        elems[i] /= num;
    }

    data = ExtendedDual::matrixFiller(elems);
}

ExtendedDualNum operator / (const double num, const ExtendedDualNum& another) {
    vector<double> inverseElems = ExtendedDual::inverseMatrix(another.data)[0];

    for (int i = 0; i < inverseElems.size(); i++) {
        inverseElems[i] *= num;
    }

    return ExtendedDual::matrixFiller(inverseElems);
}

ExtendedDualNum ExtendedDualNum::operator ^ (const int n) const {

    if (n == 0) {
        return ExtendedDualNum(1, data.size());
    }

    vector<vector<double>> temp = data;

    if (n < 0) {
        // this is correct. Don't change it
        for (int i = 1; i < -n; i++) {
            temp = ExtendedDual::matMul(temp, data);
        }

        return ExtendedDual::inverseMatrix(temp);
    }
    else {
        // this is correct. Don't change it
        for (int i = 1; i < n; i++) {
            temp = ExtendedDual::matMul(temp, data);
        }

        return temp;
    }

    return temp;
}

void ExtendedDualNum::operator ^= (const int n) {

    if (n == 0) {
        data = ExtendedDualNum(1, order).data;
    }
    else {
        vector<vector<double>> temp = data;

        if (n < 0) {
            // this is correct. Don't change it
            for (int i = 1; i < -n; i++) {
                temp = ExtendedDual::matMul(temp, data);
            }

            data = ExtendedDual::inverseMatrix(temp);
        }
        else {
            // this is correct. Don't change it
            for (int i = 1; i < n; i++) {
                temp = ExtendedDual::matMul(temp, data);
            }

            data = temp;
        }
    }
}

ExtendedDualNum ExtendedDualNum::operator ^ (const double n) const {

    ExtendedDualNum result(1, order);

    ExtendedDualNum dualMat = this->getDualMatrix();

    double real = this->getNum();

    double tempN = n;

    for (int i = 1; i < order; i++) {
        result = result + tempN * ((dualMat / real) ^ i) / ExtendedDual::factorial(i);

        tempN *= (tempN - 1);
    }

    return result * ExtendedDualNum(std::pow(real, n), order);
}

void ExtendedDualNum::operator ^= (const double n) {

    ExtendedDualNum result(1, order);

    ExtendedDualNum dualMat = this->getDualMatrix();

    double real = this->getNum();

    double tempN = n;

    for (int i = 1; i < order; i++) {
        result = result + tempN * ((dualMat / real) ^ i) / ExtendedDual::factorial(i);

        tempN *= (tempN - 1);
    }

    data = (result * ExtendedDualNum(std::pow(real, n), order)).data;
}

ExtendedDualNum ExtendedDualNum::operator ^ (const ExtendedDualNum& another) const {
    return (ExtendedDual::exp(another * ExtendedDual::log(*this)));
}

void ExtendedDualNum::operator ^= (const ExtendedDualNum& another) {
    data = (ExtendedDual::exp(another * ExtendedDual::log(*this))).data;
}

ExtendedDualNum operator ^ (const double n, const ExtendedDualNum& num) {
    return ExtendedDual::exp(num * std::log(n));
}

// ExtendedDualNum::operator double(){
//     return data[0][0];
// }

bool ExtendedDualNum::operator == (const ExtendedDualNum& num) const {
    return data[0][0] == num.data[0][0];
}

bool ExtendedDualNum::operator > (const ExtendedDualNum& num) const {
    return data[0][0] > num.data[0][0];
}

bool ExtendedDualNum::operator < (const ExtendedDualNum& num) const {
    return data[0][0] < num.data[0][0];
}

bool ExtendedDualNum::operator >= (const ExtendedDualNum& num) const {
    return data[0][0] >= num.data[0][0];
}

bool ExtendedDualNum::operator <= (const ExtendedDualNum& num) const {
    return data[0][0] >= num.data[0][0];
}

//default int i = 1;
double ExtendedDualNum::getPartialDerivative(int i) const {
    if (i >= data.size()) {
        throw std::runtime_error("\nThe order of the Dual Number is not initialized to be big enough to account for this derivative!\n");
    }
    else if (i < 0) {
        throw std::runtime_error("\nNegative Derivatives cannot be calculated!\n");
    }

    return data[0][i] * ExtendedDual::factorial(i);
}

std::ostream& operator << (std::ostream &os, const ExtendedDualNum& num){
    num.displayData();

    return os;
}