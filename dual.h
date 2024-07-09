#ifndef DUALNUMBER
#define DUALNUMBER

#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<stdexcept>
#include<thread>

using std::vector;
using std::string;

class DualNum {
private:
    long double real;
    long double dual;

public:
    DualNum(long double r = 0, long double e = 0) {
        real = r;
        dual = e;
    }

    DualNum(const DualNum& num) {
        real = num.real;
        dual = num.dual;
    }

    long double getReal() {
        return real;
    }

    long double getDual() {
        return dual;
    }

    void setReal(long double r) {
        real = r;
    }

    void setDual(long double e) {
        dual = e;
    }

    string getExpression() {
        return std::to_string(real) + " + " + std::to_string(dual) + "E";
    }

    void operator = (long double realVal) {
        real = realVal;
    }

    DualNum operator * (DualNum num) {
        return DualNum(real * num.real, (real * num.dual + dual * num.real));
    }

    void operator *= (DualNum num) {

        long double tempReal = real * num.real;

        long double tempDual = (real * num.dual + dual * num.real);

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator * (const long double k) {
        return DualNum(real * k, dual * k);
    }

    void operator *= (const long double k) {

        long double tempReal = real * k;

        long double tempDual = dual * k;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator / (DualNum num) {
        return DualNum(real / num.real, (dual * num.real - real * num.dual) / (num.real * num.real));
    }

    void operator /= (DualNum num) {

        long double tempReal = real / num.real;

        long double tempDual = (dual * num.real - real * num.dual) / (num.real * num.real);

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator / (const long double k) {
        return DualNum(real / k, dual / k);
    }

    void operator /= (const long double k) {

        long double tempReal = real / k;

        long double tempDual = dual / k;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator + (DualNum num) {
        return DualNum(real + num.real, dual + num.dual);
    }

    void operator += (DualNum num) {

        long double tempReal = real + num.real;

        long double tempDual = dual + num.dual;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator + (const long double k) {
        return DualNum(real + k, dual);
    }

    void operator += (const long double k) {

        long double tempReal = real + k;

        long double tempDual = dual;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator - (DualNum num) {
        return DualNum(real - num.real, dual - num.dual);
    }

    void operator -= (DualNum num) {

        long double tempReal = real - num.real;

        long double tempDual = dual - num.dual;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator - (const long double k) {
        return DualNum(real - k, dual);
    }

    void operator -= (const long double k) {

        long double tempReal = real - k;

        long double tempDual = dual;

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator ^ (DualNum num) {
        long double realTerm = std::pow(real, num.real);

        return DualNum(realTerm, realTerm * ((dual * num.real / real) + num.dual * std::log(real)));
    }

    void operator ^= (DualNum num) {
        long double realTerm = std::pow(real, num.real);

        long double tempReal = realTerm;

        long double tempDual = realTerm * ((dual * num.real / real) + num.dual * std::log(real));

        real = tempReal;
        dual = tempDual;
    }

    DualNum operator ^ (const long double k) {
        long double realTerm = std::pow(real, k);

        return DualNum(realTerm, realTerm * ((dual * k / real)));
    }

    void operator ^= (const long double k) {
        long double realTerm = std::pow(real, k);

        long double tempReal = realTerm;

        long double tempDual = realTerm * ((dual * k / real));

        real = tempReal;
        dual = tempDual;
    }
};

//Normal Numbers Operator Overloading

DualNum operator * (long double k, DualNum num) {
    return DualNum(k * num.getReal(), k * num.getDual());
}

DualNum operator / (long double k, DualNum num) {
    return DualNum(k / num.getReal(), -1 * (k * num.getDual()) / (num.getReal() * num.getReal()));
}

DualNum operator + (long double k, DualNum num) {
    return DualNum(k + num.getReal(), num.getDual());
}

DualNum operator - (long double k, DualNum num) {
    return DualNum(k - num.getReal(), -1 * num.getDual());
}

DualNum operator ^ (long double k, DualNum num) {
    long double realTerm = std::pow(k, num.getReal());

    return DualNum(realTerm, realTerm * num.getDual() * std::log(k));
}


//Some standard functions

namespace Dual {

    DualNum pow(DualNum x, DualNum y) {
        long double realTerm = std::pow(x.getReal(), y.getReal());

        return DualNum(realTerm, realTerm * ((x.getDual() * y.getReal() / x.getReal()) + y.getDual() * std::log(x.getReal())));
    }

    DualNum exp(DualNum x) {
        long double real = std::exp(x.getReal());

        return DualNum(real, real * x.getDual());
    }

    DualNum log(DualNum x) {
        return DualNum(std::log(x.getReal()), x.getDual() / x.getReal());
    }

    DualNum relu(DualNum x) {
        if (x.getReal() <= 0) {
            x.setReal(0);
            x.setDual(0);
        }

        return x;
    }

    DualNum tanh(DualNum x) {

        DualNum posexp = Dual::exp(x);

        DualNum negexp = Dual::exp(-1 * x);

        return ((posexp - negexp) / ((posexp + negexp)));
    }

    DualNum sigmoid(DualNum x) {
        return (1 / (1 + Dual::exp(-1 * x)));
    }

    DualNum softmax(vector<DualNum> X, int index, DualNum sum = DualNum(0, 0)) {

        if (sum.getReal() == 0) {
            vector<DualNum>::iterator iter;

            for (iter = X.begin(); iter != X.end(); iter++) {
                sum += Dual::exp(*iter);
            }
        }

        return (Dual::exp(X[index]) / sum);
    }

    DualNum mse(vector<DualNum> y_train, vector<vector<DualNum>> yhat) {
        DualNum sum(0, 0);

        for (int i = 0; i < yhat.size(); i++) {
            sum += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
        }

        return (sum / (double)(yhat.size()));
    }

    DualNum binary_crossentropy(std::vector<DualNum>& y_train, std::vector<std::vector<DualNum>>& yhat) {
        DualNum sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum = sum - (y_train[i] * Dual::log(yhat[i][0]) + (1.0 - y_train[i]) * Dual::log(1.0 - yhat[i][0]));
        }
        return (sum / (yhat.size()));
    }

    DualNum categorical_crossentropy(std::vector<DualNum>& y_train, std::vector<std::vector<DualNum>>& yhat) {
        DualNum sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum = sum - Dual::log(yhat[i][static_cast<int>(y_train[i].getReal())]);
        }
        return (sum / (yhat.size()));
    }

    DualNum accuracy(vector<DualNum> y_train, vector<vector<DualNum>> yhat, double threshold = 0.5) {
        DualNum correct(0, 0);

        for (int i = 0; i < yhat.size(); i++) {
            if (y_train[i].getReal() == 1 && (yhat[i][0].getReal() >= threshold)) {
                correct += 1;
            }
            else if (y_train[i].getReal() == 0 && (yhat[i][0].getReal() < threshold)) {
                correct += 1;
            }
        }

        return correct / yhat.size();
    }

    vector<vector<DualNum>> originalMatMul(vector<vector<DualNum>> A, vector<vector<DualNum>> B) {
        int Arows = A.size();
        int Acols = A[0].size();
        int Brows = B.size();
        int Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match first=(" + std::to_string(Arows) + ", " + std::to_string(Acols) + ") and second = (" + std::to_string(Brows) + ", " + std::to_string(Bcols) + ")");
        }

        vector<vector<DualNum>> result(Arows, vector<DualNum>(Bcols, DualNum(0)));

        for (int i = 0; i < Arows; i++) {
            for (int j = 0; j < Bcols; j++) {
                for (int k = 0; k < Brows; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return result;
    }

    void partialMatMul(vector<vector<DualNum>>& A, vector<vector<DualNum>>& B, vector<vector<DualNum>>& result, int start, int end) {
        int Arows = A.size();
        int Acols = A[0].size();
        int Brows = B.size();
        int Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match first=(" + std::to_string(Arows) + ", " + std::to_string(Acols) + ") and second = (" + std::to_string(Brows) + ", " + std::to_string(Bcols) + ")");
        }

        if (Arows >= Bcols) {
            for (int i = start; i < end; i++) {
                for (int j = 0; j < Bcols; j++) {
                    for (int k = 0; k < Brows; k++) {
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
        else {
            for (int i = 0; i < Arows; i++) {
                for (int j = start; j < end; j++) {
                    for (int k = 0; k < Brows; k++) {
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
    }

    vector<vector<DualNum>> matmul(vector<vector<DualNum>>& A, vector<vector<DualNum>>& B) {
        int Arows = A.size();
        int Acols = A[0].size();
        int Brows = B.size();
        int Bcols = B[0].size();

        int no_of_threads = std::thread::hardware_concurrency();

        int LargerDimension = (Arows >= Bcols) ? Arows : Bcols;

        if (LargerDimension < no_of_threads * 5) {
            return Dual::originalMatMul(A, B);
        }

        vector<std::thread> threads;

        vector<vector<DualNum>> result(Arows, vector<DualNum>(Bcols, DualNum(0)));


        int portion_for_each_thread = LargerDimension / no_of_threads;

        int start = 0;
        for (int i = 0; i < no_of_threads; i++) {
            int end = (i == (no_of_threads - 1) ? LargerDimension : start + portion_for_each_thread);
            threads.emplace_back(partialMatMul, std::ref(A), std::ref(B), std::ref(result), start, end);
            start = end;
        }

        for (auto& eachThread : threads) {
            eachThread.join();
        }

        return result;
    }

    void displayMatrix(vector<vector<DualNum>> mat) {
        int rows = mat.size();
        int cols = mat[0].size();

        for (int i = 0; i < rows; i++) {
            std::cout << "\n\n";
            for (int j = 0; j < cols; j++) {
                std::cout << mat[i][j].getExpression() << "\t";
            }
        }
    }

}



//some extra functions 

long double partialDerivative(DualNum(*func)(vector<DualNum>), vector<DualNum> params, int paramIndex = 0, long double at = 1) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return ((*func)(params)).getDual();
}

template<typename Func>
long double partialDerivative(Func func, vector<DualNum> params, int paramIndex = 0, long double at = 1) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return (func(params)).getDual();
}

long double partialDerivative(DualNum(*func)(DualNum), long double at = 1) {
    DualNum extX(at, 1);
    DualNum x(at, 0);

    return ((*func)(extX) - (*func)(x)).getDual();
}

DualNum evaluatePartialDerivative(DualNum(*func)(DualNum), long double at = 1) {
    DualNum x(at, 1);

    return ((*func)(x));
}

double solveUsingNewtonRaphson(DualNum(*func)(DualNum), double initialGuess = 1, int max_no_of_iterations = 10000) {

    long double threshold = 0.00001;

    double x = initialGuess;

    DualNum result;

    int iterations = 0;

    do {
        result = evaluatePartialDerivative(func, x);

        x = x - (result.getReal() / result.getDual());

        iterations++;

        if (iterations > max_no_of_iterations) {
            std::cout << "\nCouldn't converge within given no of iterations." << std::endl;
            break;
        }

    } while (result.getReal() > threshold || result.getReal() < -threshold);

    return x;
}

#endif