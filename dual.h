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

template <typename T = double>
class DualNum {
private:
    double real;
    double dual;

public:
    DualNum(double r = 0, double e = 0) {
        real = r;
        dual = e;
    }

    DualNum(const T& num) {
        real = num.real;
        dual = num.dual;
    }

    double getReal() {
        return real;
    }

    double getDual() {
        return dual;
    }

    void setReal(double r) {
        real = r;
    }

    void setDual(double e) {
        dual = e;
    }

    string getExpression() {
        return std::to_string(real) + " + " + std::to_string(dual) + "E";
    }

    void operator = (double realVal) {
        real = realVal;
    }

    T operator * (T num) {
        return T(real * num.real, (real * num.dual + dual * num.real));
    }

    void operator *= (T num) {

        double tempReal = real * num.real;

        double tempDual = (real * num.dual + dual * num.real);

        real = tempReal;
        dual = tempDual;
    }

    T operator * (const double k) {
        return T(real * k, dual * k);
    }

    void operator *= (const double k) {

        double tempReal = real * k;

        double tempDual = dual * k;

        real = tempReal;
        dual = tempDual;
    }

    T operator / (T num) {
        return T(real / num.real, (dual * num.real - real * num.dual) / (num.real * num.real));
    }

    void operator /= (T num) {

        double tempReal = real / num.real;

        double tempDual = (dual * num.real - real * num.dual) / (num.real * num.real);

        real = tempReal;
        dual = tempDual;
    }

    T operator / (const double k) {
        return T(real / k, dual / k);
    }

    void operator /= (const double k) {

        double tempReal = real / k;

        double tempDual = dual / k;

        real = tempReal;
        dual = tempDual;
    }

    T operator + (T num) {
        return T(real + num.real, dual + num.dual);
    }

    void operator += (T num) {

        double tempReal = real + num.real;

        double tempDual = dual + num.dual;

        real = tempReal;
        dual = tempDual;
    }

    T operator + (const double k) {
        return T(real + k, dual);
    }

    void operator += (const double k) {

        double tempReal = real + k;

        double tempDual = dual;

        real = tempReal;
        dual = tempDual;
    }

    T operator - (T num) {
        return T(real - num.real, dual - num.dual);
    }

    void operator -= (T num) {

        double tempReal = real - num.real;

        double tempDual = dual - num.dual;

        real = tempReal;
        dual = tempDual;
    }

    T operator - (const double k) {
        return T(real - k, dual);
    }

    void operator -= (const double k) {

        double tempReal = real - k;

        double tempDual = dual;

        real = tempReal;
        dual = tempDual;
    }

    T operator ^ (T num) {
        double realTerm = std::pow(real, num.real);

        return T(realTerm, realTerm * ((dual * num.real / real) + num.dual * std::log(real)));
    }

    void operator ^= (T num) {
        double realTerm = std::pow(real, num.real);

        double tempReal = realTerm;

        double tempDual = realTerm * ((dual * num.real / real) + num.dual * std::log(real));

        real = tempReal;
        dual = tempDual;
    }

    T operator ^ (const double k) {
        double realTerm = std::pow(real, k);

        return T(realTerm, realTerm * ((dual * k / real)));
    }

    void operator ^= (const double k) {
        double realTerm = std::pow(real, k);

        double tempReal = realTerm;

        double tempDual = realTerm * ((dual * k / real));

        real = tempReal;
        dual = tempDual;
    }
};

//Normal Numbers Operator Overloading

// Normal Numbers Operator Overloading
template<typename T = double>
DualNum<T> operator * (double k, DualNum<T> num) {
    return DualNum<T>(k * num.getReal(), k * num.getDual());
}

template<typename T = double>
DualNum<T> operator / (double k, DualNum<T> num) {
    return DualNum<T>(k / num.getReal(), -1 * (k * num.getDual()) / (num.getReal() * num.getReal()));
}

template<typename T = double>
DualNum<T> operator + (double k, DualNum<T> num) {
    return DualNum<T>(k + num.getReal(), num.getDual());
}

template<typename T = double>
DualNum<T> operator - (double k, DualNum<T> num) {
    return DualNum<T>(k - num.getReal(), -1 * num.getDual());
}

template<typename T = double>
DualNum<T> operator ^ (double k, DualNum<T> num) {
    T realTerm = std::pow(k, num.getReal());
    return DualNum<T>(realTerm, realTerm * num.getDual() * std::log(k));
}



//Some standard functions

namespace Dual {

    template<typename T = double>
    DualNum<T> pow(DualNum<T> x, DualNum<T> y) {
        T realTerm = std::pow(x.getReal(), y.getReal());
        return DualNum<T>(realTerm, realTerm * ((x.getDual() * y.getReal() / x.getReal()) + y.getDual() * std::log(x.getReal())));
    }

    template<typename T = double>
    DualNum<T> exp(DualNum<T> x) {
        T real = std::exp(x.getReal());
        return DualNum<T>(real, real * x.getDual());
    }

    template<typename T = double>
    DualNum<T> log(DualNum<T> x) {
        return DualNum<T>(std::log(x.getReal()), x.getDual() / x.getReal());
    }

    template<typename T = double>
    DualNum<T> relu(DualNum<T> x) {
        if (x.getReal() <= 0) {
            x.setReal(0);
            x.setDual(0);
        }
        return x;
    }

    template<typename T = double>
    DualNum<T> tanh(DualNum<T> x) {
        DualNum<T> posexp = Dual::exp(x);
        DualNum<T> negexp = Dual::exp(-1 * x);
        return ((posexp - negexp) / (posexp + negexp));
    }

    template<typename T = double>
    DualNum<T> sigmoid(DualNum<T> x) {
        return (1 / (1 + Dual::exp(-1 * x)));
    }

    template<typename T = double>
    DualNum<T> softmax(std::vector<DualNum<T>> X, int index, DualNum<T> sum = DualNum<T>(0, 0)) {
        if (sum.getReal() == 0) {
            for (auto& elem : X) {
                sum += Dual::exp(elem);
            }
        }
        return (Dual::exp(X[index]) / sum);
    }

    template<typename T = double>
    DualNum<T> mse(std::vector<DualNum<T>> y_train, std::vector<std::vector<DualNum<T>>> yhat) {
        DualNum<T> sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
        }
        return (sum / static_cast<T>(yhat.size()));
    }

    template<typename T = double>
    DualNum<T> accuracy(std::vector<DualNum<T>> y_train, std::vector<std::vector<DualNum<T>>> yhat, double threshold = 0.5) {
        DualNum<T> correct(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            if ((y_train[i].getReal() == 1 && yhat[i][0].getReal() >= threshold) ||
                (y_train[i].getReal() == 0 && yhat[i][0].getReal() < threshold)) {
                correct += 1;
            }
        }
        return correct / static_cast<T>(yhat.size());
    }

    template<typename T = double>
    std::vector<std::vector<DualNum<T>>> originalMatMul(std::vector<std::vector<DualNum<T>>> A, std::vector<std::vector<DualNum<T>>> B) {
        size_t Arows = A.size();
        size_t Acols = A[0].size();
        size_t Brows = B.size();
        size_t Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match.");
        }

        std::vector<std::vector<DualNum<T>>> result(Arows, std::vector<DualNum<T>>(Bcols, DualNum<T>(0)));

        for (size_t i = 0; i < Arows; i++) {
            for (size_t j = 0; j < Bcols; j++) {
                for (size_t k = 0; k < Brows; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return result;
    }

    template<typename T = double>
    void partialMatMul(std::vector<std::vector<DualNum<T>>>& A, std::vector<std::vector<DualNum<T>>>& B, std::vector<std::vector<DualNum<T>>>& result, size_t start, size_t end) {
        size_t Arows = A.size();
        size_t Acols = A[0].size();
        size_t Brows = B.size();
        size_t Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match.");
        }

        if (Arows >= Bcols) {
            for (size_t i = start; i < end; i++) {
                for (size_t j = 0; j < Bcols; j++) {
                    for (size_t k = 0; k < Brows; k++) {
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
        else {
            for (size_t i = 0; i < Arows; i++) {
                for (size_t j = start; j < end; j++) {
                    for (size_t k = 0; k < Brows; k++) {
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
    }

    template<typename T = double>
    std::vector<std::vector<DualNum<T>>> matmul(std::vector<std::vector<DualNum<T>>>& A, std::vector<std::vector<DualNum<T>>>& B) {
        size_t Arows = A.size();
        size_t Acols = A[0].size();
        size_t Brows = B.size();
        size_t Bcols = B[0].size();

        size_t no_of_threads = std::thread::hardware_concurrency();
        size_t LargerDimension = (Arows >= Bcols) ? Arows : Bcols;

        if (LargerDimension < no_of_threads * 5) {
            return Dual::originalMatMul(A, B);
        }

        std::vector<std::thread> threads;
        std::vector<std::vector<DualNum<T>>> result(Arows, std::vector<DualNum<T>>(Bcols, DualNum<T>(0)));

        size_t portion_for_each_thread = LargerDimension / no_of_threads;
        size_t start = 0;

        for (size_t i = 0; i < no_of_threads; i++) {
            size_t end = (i == (no_of_threads - 1)) ? LargerDimension : start + portion_for_each_thread;
            threads.emplace_back(partialMatMul<T>, std::ref(A), std::ref(B), std::ref(result), start, end);
            start = end;
        }

        for (auto& eachThread : threads) {
            eachThread.join();
        }

        return result;
    }

    template<typename T = double>
    void displayMatrix(const std::vector<std::vector<DualNum<T>>>& mat) {
        size_t rows = mat.size();
        size_t cols = mat[0].size();

        for (size_t i = 0; i < rows; i++) {
            std::cout << "\n\n";
            for (size_t j = 0; j < cols; j++) {
                std::cout << mat[i][j] << "\t";
            }
        }
    }

}




//some extra functions 

template<typename T = double>
double partialDerivative(DualNum<T>(*func)(std::vector<DualNum<T>>), std::vector<DualNum<T>> params, int paramIndex = 0, double at = 1) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return ((*func)(params)).getDual();
}

template<typename Func, typename T = double>
double partialDerivative(Func func, std::vector<DualNum<T>> params, int paramIndex = 0, double at = 1) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return (func(params)).getDual();
}

template<typename T = double>
double partialDerivative(DualNum<T>(*func)(DualNum<T>), double at = 1) {
    DualNum<T> extX(at, 1);
    DualNum<T> x(at, 0);

    return ((*func)(extX) - (*func)(x)).getDual();
}

template<typename T = double>
DualNum<T> evaluatePartialDerivative(DualNum<T>(*func)(DualNum<T>), double at = 1) {
    DualNum<T> x(at, 1);

    return ((*func)(x));
}


template<typename T = double>
double solveUsingNewtonRaphson(DualNum<T>(*func)(DualNum<T>), double initialGuess = 1, int max_no_of_iterations = 10000) {
    double threshold = 0.00001;
    double x = initialGuess;
    DualNum<T> result;
    int iterations = 0;

    do {
        result = evaluatePartialDerivative(func, x);
        x = x - (result.getReal() / result.getDual());
        iterations++;

        if (iterations > max_no_of_iterations) {
            std::cout << "\nCouldn't converge within given number of iterations." << std::endl;
            break;
        }

    } while (result.getReal() > threshold || result.getReal() < -threshold);

    return x;
}


#endif