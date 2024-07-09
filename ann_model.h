#ifndef ANN_MODEL
#define ANN_MODEL

#include<iostream>
#include<vector>
#include<random>
#include<stdexcept>
#include<thread>
#include "dual.h"

using std::vector;

//Some notes:

//Momentum has to be set separately and cannot be set through the .fit function, initial value = 0.1
//The value of the target during categorical_crossentropy must begin with 0, -> because it is used to index the y_train array when calculating the loss

//multi_threading is done for calculation of partial derivatives wrt various parameters. (each thread calculates wrt to some parameters and thus the work is divided !)

//multi_threaded matrix multiplication function (matmul) is also made available but not used, instead single_threaded one is used (originalMatMul)

enum class activations { linear, sigmoid, softmax, relu, tanh };
enum class lossFunction { mean_squared_error, binary_crossentropy, categorical_crossentropy };

template<typename T = double>
class ANN_Model {
private:
    std::vector<DualNum<T>> parameters;
    int inputs;
    std::vector<activations> activationFunctions;
    std::vector<int> nodes;
    std::vector<double> velocity;
    double momentum = 0.1;

public:
    ANN_Model(int inputFeatures, std::vector<int> No_of_Nodes, std::vector<activations> activationFunctionName) {
        inputs = inputFeatures;
        nodes = No_of_Nodes;
        activationFunctions = activationFunctionName;

        if (activationFunctions.size() != nodes.size()) {
            throw std::runtime_error("The number of layers and activation functions do not match.");
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distrib(-0.5, 0.5);

        long int no_of_parameters = inputFeatures * No_of_Nodes[0] + No_of_Nodes[0];

        for (int i = 0; i < No_of_Nodes.size() - 1; i++) {
            no_of_parameters += No_of_Nodes[i] * No_of_Nodes[i + 1] + No_of_Nodes[i + 1];
        }

        for (int i = 0; i < no_of_parameters; i++) {
            double random_number = distrib(gen);
            parameters.push_back(DualNum<T>(random_number));
        }

        velocity = std::vector<double>(parameters.size(), 0);
    }

    std::vector<DualNum<T>>& getParameters() {
        return parameters;
    }

    void setMomentum(double m) {
        momentum = m;
    }

    std::vector<DualNum<T>> predict(std::vector<DualNum<T>> X, std::vector<DualNum<T>> params = {}) {
        if (params.size() == 0) {
            params = parameters;
        }

        if (X.size() != inputs) {
            throw std::runtime_error("Input size does not match the number of input features.");
        }

        std::vector<DualNum<T>> inputVect = X;
        int start = 0;

        for (int i = 0; i < nodes.size(); i++) {
            std::vector<std::vector<DualNum<T>>> currentWeightMatrix;

            for (int j = 0; j < inputVect.size(); j++) {
                int end = start + nodes[i];
                std::vector<DualNum<T>> weightForEachNode(params.begin() + start, params.begin() + end);
                currentWeightMatrix.push_back(weightForEachNode);
                start = end;
            }

            std::vector<DualNum<T>> currentBiasVect(params.begin() + start, params.begin() + start + nodes[i]);
            start += nodes[i];

            std::vector<DualNum<T>> tempVect = Dual::originalMatMul(std::vector<std::vector<DualNum<T>>>(1, inputVect), currentWeightMatrix)[0];
            DualNum<T> expSums(0, 0);

            if (activationFunctions[i] == activations::softmax) {
                for (int j = 0; j < tempVect.size(); j++) {
                    expSums += Dual::exp(tempVect[j] + currentBiasVect[j]);
                }
            }

            for (int j = 0; j < tempVect.size(); j++) {
                tempVect[j] += currentBiasVect[j];

                if (activationFunctions[i] == activations::relu) {
                    tempVect[j] = Dual::relu(tempVect[j]);
                } else if (activationFunctions[i] == activations::sigmoid) {
                    tempVect[j] = Dual::sigmoid(tempVect[j]);
                } else if (activationFunctions[i] == activations::tanh) {
                    tempVect[j] = Dual::tanh(tempVect[j]);
                } else if (activationFunctions[i] == activations::softmax) {
                    tempVect[j] = Dual::softmax(tempVect, j, expSums);
                }
            }

            inputVect = tempVect;
        }

        return inputVect;
    }

    std::vector<std::vector<DualNum<T>>> predictAll(std::vector<std::vector<DualNum<T>>> data, std::vector<DualNum<T>> params = {}) {
        if (data[0].size() != inputs) {
            throw std::runtime_error("Input size does not match the number of input features.");
        }

        std::vector<std::vector<DualNum<T>>> entireOutputs;

        for (int i = 0; i < data.size(); i++) {
            entireOutputs.push_back(predict(data[i], params));
        }

        return entireOutputs;
    }

    template <typename Func>
    void updateParameters(double learning_rate, Func loss, std::vector<DualNum<T>>& newParameters, int start, int end) {
        for (int i = start; i < end; i++) {
            double gradient = partialDerivative(loss, parameters, i, parameters[i].getReal());
            velocity[i] = momentum * velocity[i] - learning_rate * gradient;
            newParameters[i] = parameters[i] + velocity[i];
        }
    }

    void fit(std::vector<std::vector<DualNum<T>>> X_train, std::vector<DualNum<T>> y_train, lossFunction lossFunc = lossFunction::mean_squared_error, double learning_rate = 0.01, int epochs = 100, int verbose = 2) {
        if (y_train.size() != X_train.size()) {
            throw std::runtime_error("The number of training samples does not match the number of labels.");
        }

        auto loss = [X_train, y_train, this, lossFunc](std::vector<DualNum<T>> params) mutable {
            std::vector<std::vector<DualNum<T>>> yhat = this->predictAll(X_train, params);
            DualNum<T> result(0, 0);

            for (int i = 0; i < yhat.size(); i++) {
                if (lossFunc == lossFunction::mean_squared_error) {
                    result += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
                } else if (lossFunc == lossFunction::binary_crossentropy) {
                    result -= ((y_train[i] * Dual::log(yhat[i][0]) + (1 - y_train[i]) * Dual::log(1 - yhat[i][0])));
                } else if (lossFunc == lossFunction::categorical_crossentropy) {
                    result -= Dual::log(yhat[i][y_train[i].getReal()]);
                }
            }

            return result / static_cast<double>(yhat.size());
        };

        int no_of_threads = std::thread::hardware_concurrency();

        for (int epoch = 1; epoch <= epochs; epoch++) {
            std::vector<DualNum<T>> newParameters(parameters.size(), DualNum<T>(0, 0));
            std::vector<std::thread> threads;

            if (no_of_threads >= parameters.size()) {
                for (int i = 0; i < parameters.size(); i++) {
                    threads.emplace_back(&ANN_Model::updateParameters<decltype(loss)>, this, learning_rate, loss, std::ref(newParameters), i, i + 1);
                }
            } else {
                int start = 0;
                int portion_for_each = parameters.size() / no_of_threads;

                for (int i = 0; i < no_of_threads; i++) {
                    int end = (i == (no_of_threads - 1) ? parameters.size() : start + portion_for_each);
                    threads.emplace_back(&ANN_Model::updateParameters<decltype(loss)>, this, learning_rate, loss, std::ref(newParameters), start, end);
                    start = end;
                }
            }

            for (auto& thread : threads) {
                thread.join();
            }

            parameters = newParameters;

            auto loss_function = Dual::mse;

            if (verbose == 2) {
                std::cout << "\nEpoch no: " << epoch << " completed!, Loss = " << loss_function(y_train, this->predictAll(X_train)).getReal() << std::endl;
            } else if (verbose == 1) {
                std::cout << "\nEpoch no: " << epoch << " completed!" << std::endl;
            } else if (verbose == 3) {
                std::vector<std::vector<DualNum<T>>> yhat = this->predictAll(X_train);
                std::cout << "\nEpoch no: " << epoch << " completed!, Loss = " << loss_function(y_train, yhat).getReal() << " accuracy = " << Dual::accuracy(y_train, yhat).getReal() << std::endl;
            }
        }
    }
};


#endif