#include<iostream>
#include<vector>
#include<random>
#include<chrono>

#include "ann_model.h"
#include "dual.h"

#include <iomanip>

using std::vector;

int main() {

    //data generation
    vector<DualNum<double>> Y;
    vector<vector<DualNum<double>>> X;

     for (int i = 0; i < 100; i++){
        Y.push_back(DualNum<double>(std::sin((float)i / 10), 0));
        vector<DualNum<double>> vectX(1, DualNum<double>((float)i / 10, 0));

        X.push_back(vectX);
    }

    //model
    ANN_Model model(1, vector<int> {16, 7, 5, 1}, vector<activations> { activations::relu, activations::relu, activations::sigmoid, activations::linear});


    //fitting the model and getting the results
    vector<vector<DualNum<double>>> yhat = model.predictAll(X);

    std::cout << "\nThe current loss is: " << Dual::mse(Y, yhat).getReal() << std::endl;

    auto start = std::chrono::steady_clock::now();

    model.fit(X, Y, lossFunction::mean_squared_error, 0.1, 100);

    auto end = std::chrono::steady_clock::now();

    std::cout << "\nTime -> " << std::chrono::duration_cast<std::chrono::milliseconds>((end - start)).count() << " milliseconds\n";

    yhat = model.predictAll(X);

    std::cout << "\nThe loss now is: " << Dual::mse(Y, yhat).getReal() << std::endl;


    return 0;
}