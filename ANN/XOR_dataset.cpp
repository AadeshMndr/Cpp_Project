#include<iostream>
#include<vector>
#include<random>
#include<chrono>

#include "ann_model.h"
#include "../Dual_Library/dual.h"
#include "../Utils/grapher.h"

#include <iomanip>

using std::vector;

class XOR_Dataset {
private:
    vector<vector<DualNum>> X_train;
    vector<DualNum> Y_train;
    vector<vector<DualNum>> X_val;
    vector<DualNum> Y_val;
    int size;

public:
    XOR_Dataset(int n) {
        size = n;

        std::random_device rd;

        std::mt19937 gen(rd());

        std::normal_distribution<> dis(0, 1);

        for (int i = 0; i < n; i++) {
            // Y_train.push_back(DualNum(std::sin((float)i / 10), 0));
            vector<DualNum> vectX = { dis(gen), dis(gen), dis(gen) };

            X_train.push_back(vectX);

            if (((vectX[0].getReal() >= 0 && vectX[1].getReal() >= 0 && vectX[2].getReal() > 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() < 0 && vectX[2].getReal() >= 0)) || (vectX[0].getReal() > 0 && vectX[1].getReal() < 0 && vectX[2].getReal() < 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() > 0 && vectX[2].getReal() <= 0)) {
                Y_train.push_back(DualNum(1, 0));
            }
            else {
                Y_train.push_back(DualNum(0, 0));
            }
        }

        for (int i = 0; i < n; i++) {
            // Y_train.push_back(DualNum(std::sin((float)i / 10), 0));
            vector<DualNum> vectX = { dis(gen), dis(gen), dis(gen) };

            X_val.push_back(vectX);

            if (((vectX[0].getReal() >= 0 && vectX[1].getReal() >= 0 && vectX[2].getReal() > 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() < 0 && vectX[2].getReal() >= 0)) || (vectX[0].getReal() > 0 && vectX[1].getReal() < 0 && vectX[2].getReal() < 0) || (vectX[0].getReal() < 0 && vectX[1].getReal() > 0 && vectX[2].getReal() <= 0)) {
                Y_val.push_back(DualNum(1, 0));
            }
            else {
                Y_val.push_back(DualNum(0, 0));
            }
        }
    }

    vector<vector<DualNum>> get_X_train() {
        return X_train;
    }

    vector<vector<DualNum>> get_X_val() {
        return X_val;
    }

    vector<DualNum> get_Y_train() {
        return Y_train;
    }

    vector<DualNum> get_Y_val() {
        return Y_val;
    }

    int getNoOfOnesInTrainingData() {
        int one = 0;
        for (int i = 0; i < size; i++) {
            if (Y_train[i].getReal() == 1) {
                one++;
            }
        }

        return one;
    }

    int getNoOfOnesInValidationData() {
        int one = 0;
        for (int i = 0; i < size; i++) {
            if (Y_train[i].getReal() == 1) {
                one++;
            }
        }

        return one;
    }
};

int main() {

    XOR_Dataset dataset(100);

    vector<vector<DualNum>> X_train = dataset.get_X_train();
    vector<DualNum> Y_train = dataset.get_Y_train();
    vector<vector<DualNum>> X_val = dataset.get_X_val();
    vector<DualNum> Y_val = dataset.get_Y_val();

    vector<vector<DualNum>> X_test = {
        {DualNum(1, 0), DualNum(1, 0), DualNum(1, 0)},
        {DualNum(-1, 0), DualNum(-1, 0), DualNum(1, 0)},
        {DualNum(1, 0), DualNum(-1, 0), DualNum(-1, 0)},
        {DualNum(-1, 0), DualNum(1, 0), DualNum(-1, 0)},
        {DualNum(1, 0), DualNum(-1, 0), DualNum(1, 0)},
        {DualNum(-1, 0), DualNum(1, 0), DualNum(1, 0)},
        {DualNum(1, 0), DualNum(1, 0), DualNum(-1, 0)},
        {DualNum(-1, 0), DualNum(-1, 0), DualNum(-1, 0)},
    };

    std::cout << "\n" << dataset.getNoOfOnesInTrainingData() << " 1s" << " are there in training dataset! \n";
    std::cout << "\n" << dataset.getNoOfOnesInValidationData() << " 1s" << " are there in Validation dataset! \n";

    //model
    ANN_Model model(3, vector<int> {10, 10, 1}, vector<activations> { activations::relu, activations::relu, activations::sigmoid});

    // model.load("XOR_Storage");

    //fitting the model and getting the results
    vector<vector<DualNum>> yhat = model.predictAll(X_train);
    vector<vector<DualNum>> yhat_val = model.predictAll(X_val);

    std::cout << "\nThe Training loss now is: " << Dual::binary_crossentropy(Y_train, yhat).getReal() << " and the accuracy is " << Dual::accuracy(Y_train, yhat).getReal() << std::endl;
    std::cout << "\nThe Validation loss now is: " << Dual::binary_crossentropy(Y_val, yhat_val).getReal() << " and the accuracy is " << Dual::accuracy(Y_val, yhat_val).getReal() << std::endl;

    vector<vector<DualNum>> y_hat_test = model.predictAll(X_test);

    std::cout << "\nThe point (+, +, +) should be 1, the model gives = " << y_hat_test[0][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, +) should be 1, the model gives = " << y_hat_test[1][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, -) should be 1, the model gives = " << y_hat_test[2][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, -) should be 1, the model gives = " << y_hat_test[3][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, +) should be 0, the model gives = " << y_hat_test[4][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, +) should be 0, the model gives = " << y_hat_test[5][0].getReal() << std::endl;
    std::cout << "\nThe point (+, +, -) should be 0, the model gives = " << y_hat_test[6][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, -) should be 0, the model gives = " << y_hat_test[7][0].getReal() << std::endl;

    // char temp;
    // std::cout << "\n\nPress Enter to start training: \n\n";
    // std::cin >> temp;

    auto start = std::chrono::steady_clock::now();

    model.fit(X_train, Y_train, X_val, Y_val, lossFunction::binary_crossentropy, 0.1, 1000, 3);

    // model.save("XOR_Storage");

    auto end = std::chrono::steady_clock::now();

    std::cout << "\nTime -> " << std::chrono::duration_cast<std::chrono::seconds>((end - start)).count() << " seconds\n";

    yhat = model.predictAll(X_train);
    yhat_val = model.predictAll(X_val);

    std::cout << "\nThe Training loss now is: " << Dual::binary_crossentropy(Y_train, yhat).getReal() << " and the accuracy is " << Dual::accuracy(Y_train, yhat).getReal() << std::endl;
    std::cout << "\nThe Validation loss now is: " << Dual::binary_crossentropy(Y_val, yhat_val).getReal() << " and the accuracy is " << Dual::accuracy(Y_val, yhat_val).getReal() << std::endl;

    y_hat_test = model.predictAll(X_test);

    std::cout << "\nThe point (+, +, +) should be 1, the model gives = " << y_hat_test[0][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, +) should be 1, the model gives = " << y_hat_test[1][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, -) should be 1, the model gives = " << y_hat_test[2][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, -) should be 1, the model gives = " << y_hat_test[3][0].getReal() << std::endl;
    std::cout << "\nThe point (+, -, +) should be 0, the model gives = " << y_hat_test[4][0].getReal() << std::endl;
    std::cout << "\nThe point (-, +, +) should be 0, the model gives = " << y_hat_test[5][0].getReal() << std::endl;
    std::cout << "\nThe point (+, +, -) should be 0, the model gives = " << y_hat_test[6][0].getReal() << std::endl;
    std::cout << "\nThe point (-, -, -) should be 0, the model gives = " << y_hat_test[7][0].getReal() << std::endl;

    Graph plotter;

    char x[] = "Epochs";
    char y[] = "Loss";

    plotter.setLegend({ "Training Loss", "Validation Loss" }, { MyColors::point1, MyColors::point2 });

    plotter.pushColoredPoints(model.getLossHistory(), MyColors::point1);
    plotter.pushColoredValPoints(model.getValLossHistory(), MyColors::point2);

    plotter.start(x, y);

    // std::cout << "Number -> " << model.getValAccuracyHistory().size() << "\n\n";

    Graph accuracyPlotter;

    char accX[] = "Epochs";
    char accY[] = "Accuracy";

    accuracyPlotter.setLegend({ "Training Accuracy", "Validation Accuracy" }, { MyColors::point1, MyColors::point2 });

    accuracyPlotter.pushColoredPoints(model.getAccuracyHistory(), MyColors::point1);
    accuracyPlotter.pushColoredValPoints(model.getValAccuracyHistory(), MyColors::point2);

    accuracyPlotter.start(accX, accY);

    return 0;
}
