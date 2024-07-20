#include<iostream>
#include<vector>
#include<random>
#include<chrono>

#include "ann_model.h"
#include "../Dual_Library/dual.h"
#include "../Utils/grapher.h"
// #include "rough.h"

#include <iomanip>

using std::vector;

int main() {

    //data generation
    vector<DualNum> Y;
    vector<vector<DualNum>> X;

    for (int i = 0; i < 100; i++) {
        Y.push_back(DualNum(std::sin((float)i / 10), 0));
        vector<DualNum> vectX(1, DualNum((float)i / 10, 0));

        X.push_back(vectX);
    }

    vector<DualNum> Y_val;
    vector<vector<DualNum>> X_val;

    for (int i = 0; i < 100; i++) {
        Y_val.push_back(DualNum(std::sin((float)i / 8), 0));
        vector<DualNum> vectX(1, DualNum((float)i / 8, 0));

        X_val.push_back(vectX);
    }

    vector<double> doubleY;
    for (int i = 0; i < Y.size(); i++) {
        doubleY.push_back(Y[i].getReal());
    }

    //model
    ANN_Model model(1, vector<int> {16, 7, 5, 1}, vector<activations> { activations::relu, activations::relu, activations::sigmoid, activations::linear});


    //fitting the model and getting the results
    vector<vector<DualNum>> yhat = model.predictAll(X);

    std::cout << "\nThe current training loss is: " << Dual::mse(Y, yhat).getReal() << std::endl;

    auto start = std::chrono::steady_clock::now();

    model.fit(X, Y, X_val, Y_val, lossFunction::mean_squared_error, 0.1, 1000, 2);

    auto end = std::chrono::steady_clock::now();

    std::cout << "\nTime -> " << std::chrono::duration_cast<std::chrono::milliseconds>((end - start)).count() << " milliseconds\n";

    yhat = model.predictAll(X);

    std::cout << "\nThe training loss now is: " << Dual::mse(Y, yhat).getReal() << std::endl;

    vector<double> doubleYhat;
    for (int i = 0; i < Y.size(); i++) {
        doubleYhat.push_back(yhat[i][0].getReal());
    }

    Graph plotter;

    char x[] = "Epochs";
    char y[] = "Loss";

    plotter.setLegend({ "Training Loss", "Validation Loss" }, { MyColors::point1, MyColors::point2 });

    plotter.pushColoredPoints(model.getLossHistory(), MyColors::point1);
    plotter.pushColoredValPoints(model.getValLossHistory(), MyColors::point2);

    plotter.start(x, y);

    char sineX[] = "X";
    char sineY[] = "Y";

    Graph sinewave;

    sinewave.pushColoredPoints(doubleY, MyColors::point1);

    sinewave.pushColoredValPoints(doubleYhat, MyColors::point2);

    sinewave.setLegend({ "Actual", "Predicted" }, { MyColors::point1, MyColors::point2 });

    sinewave.start(sineX, sineY);


    return 0;
}