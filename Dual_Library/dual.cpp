#include "dual.h"

// DualNum member functions
DualNum::DualNum(long double r, long double e) {
    real = r;
    dual = e;
}

DualNum::DualNum(const DualNum& num) {
    real = num.real;
    dual = num.dual;
}

long double DualNum::getReal() const {
    return real;
}

long double DualNum::getDual() const {
    return dual;
}

void DualNum::setReal(const long double& r) {
    real = r;
}

void DualNum::setDual(const long double& e) {
    dual = e;
}

string DualNum::getExpression() const {
    return std::to_string(real) + " + " + std::to_string(dual) + "E";
}

void DualNum::operator = (const long double& realVal) {
    real = realVal;
}

DualNum DualNum::operator * (const DualNum& num) const {
    return DualNum(real * num.real, (real * num.dual + dual * num.real));
}

void DualNum::operator *= (const DualNum& num) {

    long double tempReal = real * num.real;

    long double tempDual = (real * num.dual + dual * num.real);

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator * (const long double& k) const {
    return DualNum(real * k, dual * k);
}

void DualNum::operator *= (const long double& k) {

    long double tempReal = real * k;

    long double tempDual = dual * k;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator / (const DualNum& num) const {
    return DualNum(real / num.real, (dual * num.real - real * num.dual) / (num.real * num.real));
}

void DualNum::operator /= (const DualNum& num) {

    long double tempReal = real / num.real;

    long double tempDual = (dual * num.real - real * num.dual) / (num.real * num.real);

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator / (const long double& k) const {
    return DualNum(real / k, dual / k);
}

void DualNum::operator /= (const long double& k) {

    long double tempReal = real / k;

    long double tempDual = dual / k;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator + (const DualNum& num) const {
    return DualNum(real + num.real, dual + num.dual);
}

void DualNum::operator += (const DualNum& num) {

    long double tempReal = real + num.real;

    long double tempDual = dual + num.dual;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator + (const long double& k) const {
    return DualNum(real + k, dual);
}

void DualNum::operator += (const long double& k) {

    long double tempReal = real + k;

    long double tempDual = dual;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator - (const DualNum& num) const {
    return DualNum(real - num.real, dual - num.dual);
}

void DualNum::operator -= (const DualNum& num) {

    long double tempReal = real - num.real;

    long double tempDual = dual - num.dual;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator - (const long double& k) const {
    return DualNum(real - k, dual);
}

void DualNum::operator -= (const long double& k) {

    long double tempReal = real - k;

    long double tempDual = dual;

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator ^ (const DualNum& num) const {
    long double realTerm = std::pow(real, num.real);

    return DualNum(realTerm, realTerm * ((dual * num.real / real) + num.dual * std::log(real)));
}

void DualNum::operator ^= (const DualNum& num) {
    long double realTerm = std::pow(real, num.real);

    long double tempReal = realTerm;

    long double tempDual = realTerm * ((dual * num.real / real) + num.dual * std::log(real));

    real = tempReal;
    dual = tempDual;
}

DualNum DualNum::operator ^ (const long double& k) const {
    long double realTerm = std::pow(real, k);

    return DualNum(realTerm, realTerm * ((dual * k / real)));
}

void DualNum::operator ^= (const long double& k) {
    long double realTerm = std::pow(real, k);

    long double tempReal = realTerm;

    long double tempDual = realTerm * ((dual * k / real));

    real = tempReal;
    dual = tempDual;
}

// Normal Numbers Operator Overloading
DualNum operator * (const long double& k, const DualNum& num) {
    return DualNum(k * num.getReal(), k * num.getDual());
}

DualNum operator / (const long double& k, const DualNum& num) {
    return DualNum(k / num.getReal(), -1 * (k * num.getDual()) / (num.getReal() * num.getReal()));
}

DualNum operator + (const long double& k, const DualNum& num) {
    return DualNum(k + num.getReal(), num.getDual());
}

DualNum operator - (const long double& k, const DualNum& num) {
    return DualNum(k - num.getReal(), -1 * num.getDual());
}

DualNum operator ^ (const long double& k, const DualNum& num) {
    long double realTerm = std::pow(k, num.getReal());

    return DualNum(realTerm, realTerm * num.getDual() * std::log(k));
}

// Dual standard functions
DualNum Dual::pow(const DualNum& x, const DualNum& y) {
    long double realTerm = std::pow(x.getReal(), y.getReal());

    return DualNum(realTerm, realTerm * ((x.getDual() * y.getReal() / x.getReal()) + y.getDual() * std::log(x.getReal())));
}

DualNum Dual::exp(const DualNum& x) {
    long double realTerm = std::exp(x.getReal());

    return DualNum(realTerm, realTerm * x.getDual());
}

DualNum Dual::log(const DualNum& x) {

    return DualNum(std::log(x.getReal()), x.getDual() / x.getReal());
}

DualNum Dual::relu(const DualNum& x) {
    DualNum y = x;

    if (y.getReal() <= 0) {
        y.setReal(0);
        y.setDual(0);
    }

    return y;

    // return DualNum(std::fmax(x.getReal(), 0), x.getReal() > 0 ? x.getDual() : 0);
}

DualNum Dual::tanh(const DualNum& x) {
    DualNum posexp = Dual::exp(x);

    DualNum negexp = Dual::exp(-1 * x);

    return ((posexp - negexp) / ((posexp + negexp)));

    // return DualNum(std::tanh(x.getReal()), x.getDual() / std::pow(std::cosh(x.getReal()), 2));
}

DualNum Dual::sigmoid(const DualNum& x) {
    return (1 / (1 + Dual::exp(-1 * x)));
    // return DualNum(1 / (1 + std::exp(-x.getReal())), std::exp(-x.getReal()) * x.getDual() / std::pow(1 + std::exp(-x.getReal()), 2));
}

DualNum Dual::softmax(const vector<DualNum>& X, int index, DualNum sum) {

    if (sum.getReal() == 0) {

        for (int i = 0; i < X.size(); i++) {
            sum += Dual::exp(X[i]);
        }
    }

    return (Dual::exp(X[index]) / sum);
}

DualNum Dual::mse(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat) {
    DualNum sum(0, 0);

    for (int i = 0; i < yhat.size(); i++) {
        sum += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
    }

    return (sum / (double)(yhat.size()));
}

DualNum Dual::binary_crossentropy(const vector<DualNum>& y_train, const vector<std::vector<DualNum>>& yhat) {
    DualNum sum(0, 0);
    for (int i = 0; i < yhat.size(); i++) {
        sum = sum - (y_train[i] * Dual::log(yhat[i][0]) + (1.0 - y_train[i]) * Dual::log(1.0 - yhat[i][0]));
    }

    return (sum / (yhat.size()));
}

DualNum Dual::categorical_crossentropy(const vector<DualNum>& y_train, const vector<std::vector<DualNum>>& yhat) {
    DualNum sum(0, 0);
    for (int i = 0; i < yhat.size(); i++) {
        sum = sum - Dual::log(yhat[i][static_cast<int>(y_train[i].getReal())]);
    }
    return (sum / (static_cast<int>(yhat.size())));
}

DualNum Dual::accuracy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat, const double threshold) {
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

vector<vector<DualNum>> Dual::originalMatMul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B) {
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

void Dual::partialMatMul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B, vector<vector<DualNum>>& result, int start, int end) {
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

vector<vector<DualNum>> Dual::matmul(const vector<vector<DualNum>>& A, const vector<vector<DualNum>>& B) {
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

void Dual::displayMatrix(const vector<vector<DualNum>>& mat) {
    int rows = mat.size();
    int cols = mat[0].size();

    for (int i = 0; i < rows; i++) {
        std::cout << "\n\n";
        for (int j = 0; j < cols; j++) {
            std::cout << mat[i][j].getExpression() << "\t";
        }
    }
}

// Extra Functions
long double partialDerivative(DualNum(*func)(vector<DualNum>), vector<DualNum> params, int paramIndex, long double at) {
    params[paramIndex].setDual(1);

    if (at != 1) {
        params[paramIndex].setReal(at);
    }

    return ((*func)(params)).getDual();
}

long double partialDerivative(DualNum(*func)(DualNum), long double at) {
    DualNum extX(at, 1);
    DualNum x(at, 0);

    return ((*func)(extX) - (*func)(x)).getDual();
}

DualNum evaluatePartialDerivative(DualNum(*func)(DualNum), long double at) {
    DualNum x(at, 1);

    return ((*func)(x));
}

vector<DualNum> gradient(DualNum(*func)(vector<DualNum>), vector<DualNum> at, coordinate_system system) {

    vector<DualNum> result(at.size());

    if (system == coordinate_system::cartesian) {
        for (int i = 0; i < at.size(); i++) {
            at[i].setDual(1);

            result[i] = partialDerivative(func, at, i, at[i].getReal());

            at[i].setDual(0);
        }
    }

    return result;
}

double solveUsingNewtonRaphson(DualNum(*func)(DualNum), double initialGuess, int max_no_of_iterations) {

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

    } while ((result.getReal() > threshold) || (result.getReal() < -threshold));

    return x;
}
