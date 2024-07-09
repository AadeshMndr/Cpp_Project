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

template <typename T = double>
class DualNum {
private:
    T real;
    T dual;

public:
    DualNum(T r = 0, T e = 0) : real(r), dual(e) {}

    DualNum(const DualNum& num) : real(num.real), dual(num.dual) {}

    T getReal() const {
        return real;
    }

    T getDual() const {
        return dual;
    }

    void setReal(T r) {
        real = r;
    }

    void setDual(T e) {
        dual = e;
    }

    string getExpression() const {
        return std::to_string(real) + " + " + std::to_string(dual) + "E";
    }

    void operator=(T realVal) {
        real = realVal;
    }

    DualNum operator*(const DualNum& num) const {
        return DualNum(real * num.real, (real * num.dual + dual * num.real));
    }

    void operator*=(const DualNum& num) {
        T tempReal = real * num.real;
        T tempDual = (real * num.dual + dual * num.real);
        real = tempReal;
        dual = tempDual;
    }

    DualNum operator*(const T k) const {
        return DualNum(real * k, dual * k);
    }

    void operator*=(const T k) {
        real *= k;
        dual *= k;
    }

    DualNum operator/(const DualNum& num) const {
        return DualNum(real / num.real, (dual * num.real - real * num.dual) / (num.real * num.real));
    }

    void operator/=(const DualNum& num) {
        T tempReal = real / num.real;
        T tempDual = (dual * num.real - real * num.dual) / (num.real * num.real);
        real = tempReal;
        dual = tempDual;
    }

    DualNum operator/(const T k) const {
        return DualNum(real / k, dual / k);
    }

    void operator/=(const T k) {
        real /= k;
        dual /= k;
    }

    DualNum operator+(const DualNum& num) const {
        return DualNum(real + num.real, dual + num.dual);
    }

    void operator+=(const DualNum& num) {
        real += num.real;
        dual += num.dual;
    }

    DualNum operator+(const T k) const {
        return DualNum(real + k, dual);
    }

    void operator+=(const T k) {
        real += k;
    }

    DualNum operator-(const DualNum& num) const {
        return DualNum(real - num.real, dual - num.dual);
    }

    void operator-=(const DualNum& num) {
        real -= num.real;
        dual -= num.dual;
    }

    DualNum operator-(const T k) const {
        return DualNum(real - k, dual);
    }

    void operator-=(const T k) {
        real -= k;
    }

    DualNum operator^(const DualNum& num) const {
        T realTerm = std::pow(real, num.real);
        return DualNum(realTerm, realTerm * ((dual * num.real / real) + num.dual * std::log(real)));
    }

    void operator^=(const DualNum& num) {
        T realTerm = std::pow(real, num.real);
        real = realTerm;
        dual = realTerm * ((dual * num.real / real) + num.dual * std::log(real));
    }

    DualNum operator^(const T k) const {
        T realTerm = std::pow(real, k);
        return DualNum(realTerm, realTerm * ((dual * k / real)));
    }

    void operator^=(const T k) {
        T realTerm = std::pow(real, k);
        real = realTerm;
        dual = realTerm * ((dual * k / real));
    }
};

// Normal Numbers Operator Overloading

template <typename T>
DualNum<T> operator*(T k, const DualNum<T>& num) {
    return DualNum<T>(k * num.getReal(), k * num.getDual());
}

template <typename T>
DualNum<T> operator/(T k, const DualNum<T>& num) {
    return DualNum<T>(k / num.getReal(), -1 * (k * num.getDual()) / (num.getReal() * num.getReal()));
}

template <typename T>
DualNum<T> operator+(T k, const DualNum<T>& num) {
    return DualNum<T>(k + num.getReal(), num.getDual());
}

template <typename T>
DualNum<T> operator-(T k, const DualNum<T>& num) {
    return DualNum<T>(k - num.getReal(), -1 * num.getDual());
}

template <typename T>
DualNum<T> operator^(T k, const DualNum<T>& num) {
    T realTerm = std::pow(k, num.getReal());
    return DualNum<T>(realTerm, realTerm * num.getDual() * std::log(k));
}

// Some standard functions

namespace Dual {

    DualNum<double> pow(const DualNum<double>& x, const DualNum<double>& y) {
        double realTerm = std::pow(x.getReal(), y.getReal());
        return DualNum<double>(realTerm, realTerm * ((x.getDual() * y.getReal() / x.getReal()) + y.getDual() * std::log(x.getReal())));
    }

    DualNum<double> exp(const DualNum<double>& x) {
        double real = std::exp(x.getReal());
        return DualNum<double>(real, real * x.getDual());
    }

    DualNum<double> log(const DualNum<double>& x) {
        return DualNum<double>(std::log(x.getReal()), x.getDual() / x.getReal());
    }

    template <typename T>
    DualNum<T> relu(DualNum<T> x) {
        if (x.getReal() <= 0) {
            x.setReal(0);
            x.setDual(0);
        }
        return x;
    }

    template <typename T>
    DualNum<T> tanh(DualNum<T> x) {
        DualNum<T> posexp = Dual::exp(x);
        DualNum<T> negexp = Dual::exp(-1.0 * x);
        return ((posexp - negexp) / ((posexp + negexp)));
    }

    template <typename T>
    DualNum<T> sigmoid(DualNum<T> x) {
        return (1.0 / (1.0 + Dual::exp(-1.0 * x)));
    }

    template <typename T>
    DualNum<T> softmax(const std::vector<DualNum<T>>& X, int index, DualNum<T> sum = DualNum<T>(0, 0)) {
        if (sum.getReal() == 0) {
            for (const auto& elem : X) {
                sum += Dual::exp(elem);
            }
        }
        return (Dual::exp(X[index]) / sum);
    }

    template <typename T>
    DualNum<T> mse(const std::vector<DualNum<T>>& y_train, const std::vector<std::vector<DualNum<T>>>& yhat) {
        DualNum<T> sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
        }
        return (sum / static_cast<T>(yhat.size()));
    }

    template <typename T>
    DualNum<T> binary_crossentropy(const std::vector<DualNum<T>>& y_train, const std::vector<std::vector<DualNum<T>>>& yhat) {
        DualNum<T> sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum = sum - (y_train[i] * Dual::log(yhat[i][0]) + (DualNum<T>(1.0) - y_train[i]) * Dual::log(DualNum<T>(1.0) - yhat[i][0]));
        }
        return (sum / static_cast<T>(yhat.size()));
    }

    template <typename T>
    DualNum<T> categorical_crossentropy(const std::vector<DualNum<T>>& y_train, const std::vector<std::vector<DualNum<T>>>& yhat) {
        DualNum<T> sum(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            sum = sum - Dual::log(yhat[i][static_cast<int>(y_train[i].getReal())]);
        }
        return (sum / static_cast<T>(yhat.size()));
    }



    template <typename T>
    DualNum<T> accuracy(const std::vector<DualNum<T>>& y_train, const std::vector<std::vector<DualNum<T>>>& yhat, T threshold = 0.5) {
        DualNum<T> correct(0, 0);
        for (size_t i = 0; i < yhat.size(); i++) {
            if ((y_train[i].getReal() == 1 && yhat[i][0].getReal() >= threshold) ||
                (y_train[i].getReal() == 0 && yhat[i][0].getReal() < threshold)) {
                correct += 1;
            }
        }
        return correct / static_cast<T>(yhat.size());
    }

    template <typename T>
    std::vector<std::vector<DualNum<T>>> originalMatMul(const std::vector<std::vector<DualNum<T>>>& A, const std::vector<std::vector<DualNum<T>>>& B) {
        size_t Arows = A.size();
        size_t Acols = A[0].size();
        size_t Brows = B.size();
        size_t Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match first=(" + std::to_string(Arows) + ", " + std::to_string(Acols) + ") and second = (" + std::to_string(Brows) + ", " + std::to_string(Bcols) + ")");
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

    template <typename T>
    void partialMatMul(const std::vector<std::vector<DualNum<T>>>& A, const std::vector<std::vector<DualNum<T>>>& B, std::vector<std::vector<DualNum<T>>>& result, size_t start, size_t end) {
        size_t Arows = A.size();
        size_t Acols = A[0].size();
        size_t Brows = B.size();
        size_t Bcols = B[0].size();

        if (Acols != Brows) {
            throw std::runtime_error("The shapes of the Matrices don't match first=(" + std::to_string(Arows) + ", " + std::to_string(Acols) + ") and second = (" + std::to_string(Brows) + ", " + std::to_string(Bcols) + ")");
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

    template <typename T>
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
            size_t end = (i == (no_of_threads - 1) ? LargerDimension : start + portion_for_each_thread);
            threads.emplace_back(partialMatMul<T>, std::ref(A), std::ref(B), std::ref(result), start, end);
            start = end;
        }

        for (auto& eachThread : threads) {
            eachThread.join();
        }

        return result;
    }

    template <typename T>
    void displayMatrix(const std::vector<std::vector<DualNum<T>>>& mat) {
        size_t rows = mat.size();
        size_t cols = mat[0].size();

        for (size_t i = 0; i < rows; i++) {
            std::cout << "\n\n";
            for (size_t j = 0; j < cols; j++) {
                std::cout << mat[i][j].getExpression() << "\t";
            }
        }
    }

}

// Some extra functions

template <typename T>
T partialDerivative(DualNum<T>(*func)(const std::vector<DualNum<T>>&), std::vector<DualNum<T>> params, size_t paramIndex = 0, T at = 1) {
    params[paramIndex].setDual(1);
    if (at != 1) {
        params[paramIndex].setReal(at);
    }
    return ((*func)(params)).getDual();
}

template <typename Func, typename T>
T partialDerivative(Func func, std::vector<DualNum<T>> params, size_t paramIndex = 0, T at = 1) {
    params[paramIndex].setDual(1);
    if (at != 1) {
        params[paramIndex].setReal(at);
    }
    return (func(params)).getDual();
}

template <typename T>
T partialDerivative(DualNum<T>(*func)(DualNum<T>), T at = 1) {
    DualNum<T> extX(at, 1);
    DualNum<T> x(at, 0);
    return ((*func)(extX) - (*func)(x)).getDual();
}

template <typename T>
DualNum<T> evaluatePartialDerivative(DualNum<T>(*func)(DualNum<T>), T at = 1) {
    DualNum<T> x(at, 1);
    return ((*func)(x));
}

template <typename T>
T solveUsingNewtonRaphson(DualNum<T>(*func)(DualNum<T>), T initialGuess = 1, size_t max_no_of_iterations = 10000) {
    T threshold = 0.00001;
    T x = initialGuess;
    DualNum<T> result;
    size_t iterations = 0;

    do {
        result = evaluatePartialDerivative(func, x);
        x = x - (result.getReal() / result.getDual());
        iterations++;

        if (iterations > max_no_of_iterations) {
            std::cout << "\nCouldn't converge within given no of iterations." << std::endl;
            break;
        }
    } while (std::abs(result.getReal()) > threshold);

    return x;
}

#endif
