#include<iostream>
#include<vector>
#include<random>
#include<chrono>

#include "ann_model.h"
#include "dual.h"

#include <iomanip>

using std::vector;

int main() {

    std::random_device rd;

    std::mt19937 gen(rd());

    std::normal_distribution<> dis(0, 1);

    //data generation
    vector<DualNum<double>> Y;
    vector<vector<DualNum<double>>> X;

    int n = 100;

     for (int i = 0; i < n; i++){
        // Y.push_back(DualNum<double>(std::sin((float)i / 10), 0));
        vector<DualNum<double>> vectX = {dis(gen), dis(gen), dis(gen)};

        X.push_back(vectX);

        if (((vectX[0].getReal() >= 0 && vectX[1].getReal() >= 0 && vectX[2].getReal() > 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() < 0 && vectX[2].getReal() >= 0)) || (vectX[0].getReal() > 0 && vectX[1].getReal() < 0 && vectX[2].getReal() < 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() > 0 && vectX[2].getReal() <= 0)){
            Y.push_back(DualNum<double>(1, 0));
        } else {
            Y.push_back(DualNum<double>(0, 0));
        }
    }

    int one = 0;
    for (int i = 0; i < n; i++){
        if (Y[i].getReal() == 1){
            one++;
        }

    }

    vector<vector<DualNum<double>>> X_test = {
        {DualNum<double>(1, 0), DualNum<double>(1, 0), DualNum<double>(1, 0)},
        {DualNum<double>(-1, 0), DualNum<double>(-1, 0), DualNum<double>(1, 0)},
        {DualNum<double>(1, 0), DualNum<double>(-1, 0), DualNum<double>(-1, 0)},
        {DualNum<double>(-1, 0), DualNum<double>(1, 0), DualNum<double>(-1, 0)},
        {DualNum<double>(1, 0), DualNum<double>(-1, 0), DualNum<double>(1, 0)},
        {DualNum<double>(-1, 0), DualNum<double>(1, 0), DualNum<double>(1, 0)},
        {DualNum<double>(1, 0), DualNum<double>(1, 0), DualNum<double>(-1, 0)},
        {DualNum<double>(-1, 0), DualNum<double>(-1, 0), DualNum<double>(-1, 0)},
    };

    std::cout << "\n" << one << " 1s" << " are there in this dataset! \n";

    //model
    ANN_Model<double> model(3, vector<int> {10, 10, 1}, vector<activations> { activations::relu, activations::relu, activations::sigmoid});


    //fitting the model and getting the results
    vector<vector<DualNum<double>>> yhat = model.predictAll(X);

    std::cout << "\nThe current loss is: " << Dual::mse(Y, yhat).getReal() << std::endl;

    vector<vector<DualNum<double>>> y_hat_test = model.predictAll(X_test);

    std::cout << "\nThe point (+, +, +) should be 1, the model gives = " << y_hat_test[0][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, +) should be 1, the model gives = " << y_hat_test[1][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, -) should be 1, the model gives = " << y_hat_test[2][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, -) should be 1, the model gives = " << y_hat_test[3][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, +) should be 0, the model gives = " << y_hat_test[4][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, +) should be 0, the model gives = " << y_hat_test[5][0].getReal() << std::endl;
    std::cout << "\nThe point (+, +, -) should be 0, the model gives = " << y_hat_test[6][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, -) should be 0, the model gives = " << y_hat_test[7][0].getReal() << std::endl;

    auto start = std::chrono::steady_clock::now();

    model.fit(X, Y, lossFunction::binary_crossentropy, 0.1, 1000, 3);

    auto end = std::chrono::steady_clock::now();

    std::cout << "\nTime -> " << std::chrono::duration_cast<std::chrono::milliseconds>((end - start)).count() << " milliseconds\n";

    yhat = model.predictAll(X);

    std::cout << "\nThe loss now is: " << Dual::mse(Y, yhat).getReal() << " and the accuracy is " << Dual::accuracy(Y, yhat).getReal() << std::endl;

    y_hat_test = model.predictAll(X_test);

    std::cout << "\nThe point (+, +, +) should be 1, the model gives = " << y_hat_test[0][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, +) should be 1, the model gives = " << y_hat_test[1][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, -) should be 1, the model gives = " << y_hat_test[2][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, -) should be 1, the model gives = " << y_hat_test[3][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, +) should be 0, the model gives = " << y_hat_test[4][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, +) should be 0, the model gives = " << y_hat_test[5][0].getReal() << std::endl;
    std::cout << "\nThe point (+, +, -) should be 0, the model gives = " << y_hat_test[6][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, -) should be 0, the model gives = " << y_hat_test[7][0].getReal() << std::endl;

    return 0;
}