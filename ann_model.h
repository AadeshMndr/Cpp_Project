#ifndef ANN_MODEL
#define ANN_MODEL

#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>
#include <thread>
#include "dual.h"

using std::vector;

enum class activations { linear, sigmoid, softmax, relu, tanh };
enum class lossFunction { mean_squared_error, binary_crossentropy, categorical_crossentropy };

template <typename T>
class ANN_Model {
private:
    vector<DualNum<T>> parameters;
    int inputs;
    vector<activations> activationFunctions;
    vector<int> nodes;
    vector<T> velocity;
    T momentum = 0.1;

public:
    ANN_Model(int inputFeatures, vector<int> No_of_Nodes, vector<activations> activationFunctionName) {
        inputs = inputFeatures;
        nodes = No_of_Nodes;
        activationFunctions = activationFunctionName;

        if (activationFunctions.size() != nodes.size()) {
            throw std::runtime_error("The No. of Layers and the No. of activation functions don't match!");
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distrib(-0.5, 0.5);

        long int no_of_parameters = inputFeatures * No_of_Nodes[0] + No_of_Nodes[0];

        for (size_t i = 0; i < No_of_Nodes.size() - 1; i++) {
            no_of_parameters += No_of_Nodes[i] * No_of_Nodes[i + 1] + No_of_Nodes[i + 1];
        }

        for (long int i = 0; i < no_of_parameters; i++) {
            T random_number = distrib(gen);
            parameters.push_back(DualNum<T>(random_number));
        }

        velocity = vector<T>(parameters.size(), 0);
    }

    vector<DualNum<T>>& getParameters() {
        return parameters;
    }

    void setMomentum(T m) {
        momentum = m;
    }

    vector<DualNum<T>> predict(vector<DualNum<T>> X, vector<DualNum<T>> params = {}) {
        if (params.empty()) {
            params = parameters;
        }

        if (X.size() != inputs) {
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        vector<DualNum<T>> inputVect = X;
        int start = 0;

        for (size_t i = 0; i < nodes.size(); i++) {
            vector<vector<DualNum<T>>> currentWeightMatrix;

            for (size_t j = 0; j < inputVect.size(); j++) {
                int end = start + nodes[i];
                vector<DualNum<T>> weightForEachNode(params.begin() + start, params.begin() + end);
                currentWeightMatrix.push_back(weightForEachNode);
                start = end;
            }

            vector<DualNum<T>> currentBiasVect(params.begin() + start, params.begin() + start + nodes[i]);
            start += nodes[i];

            vector<DualNum<T>> tempVect = Dual::originalMatMul(vector<vector<DualNum<T>>>(1, inputVect), currentWeightMatrix)[0];
            DualNum<T> expSums(0, 0);

            if (activationFunctions[i] == activations::softmax) {
                for (size_t j = 0; j < tempVect.size(); j++) {
                    expSums = expSums + Dual::exp(tempVect[j] + currentBiasVect[j]);
                }
            }

            for (size_t j = 0; j < tempVect.size(); j++) {
                tempVect[j] = tempVect[j] + currentBiasVect[j];

                switch (activationFunctions[i]) {
                    case activations::relu:
                        tempVect[j] = Dual::relu(tempVect[j]);
                        break;
                    case activations::sigmoid:
                        tempVect[j] = Dual::sigmoid(tempVect[j]);
                        break;
                    case activations::tanh:
                        tempVect[j] = Dual::tanh(tempVect[j]);
                        break;
                    case activations::softmax:
                        tempVect[j] = Dual::softmax(tempVect, j, expSums);
                        break;
                    default:
                        break;
                }
            }

            inputVect = tempVect;
        }

        return inputVect;
    }

    vector<vector<DualNum<T>>> predictAll(vector<vector<DualNum<T>>> data, vector<DualNum<T>> params = {}) {
        if (data[0].size() != inputs) {
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        vector<vector<DualNum<T>>> entireOutputs;
        for (const auto& datum : data) {
            entireOutputs.push_back(predict(datum, params));
        }

        return entireOutputs;
    }

    template <typename Func>
    void updateParameters(T learning_rate, Func loss, vector<DualNum<T>>& newParameters, int start, int end) {
        for (int i = start; i < end; i++) {
            T gradient = partialDerivative(loss, parameters, i, parameters[i].getReal());
            velocity[i] = momentum * velocity[i] - learning_rate * gradient;
            newParameters[i] = parameters[i] + velocity[i];
        }
    }

    void fit(vector<vector<DualNum<T>>> X_train, vector<DualNum<T>> y_train, lossFunction lossFunc = lossFunction::mean_squared_error, T learning_rate = 0.01, int epochs = 100, int verbose = 2) {
        if (y_train.size() != X_train.size()) {
            throw std::runtime_error("The shape of the X_train and y_train don't match!");
        }

        auto loss = [X_train, y_train, this, lossFunc](vector<DualNum<T>> params) mutable {
            vector<vector<DualNum<T>>> yhat = this->predictAll(X_train, params);
            DualNum<T> result(0, 0);

            for (size_t i = 0; i < yhat.size(); i++) {
                switch (lossFunc) {
                    case lossFunction::mean_squared_error:
                        result = result + (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
                        break;
                    case lossFunction::binary_crossentropy:
                        result = result - (y_train[i] * Dual::log(yhat[i][0]) + (DualNum<T>(1.0) - y_train[i]) * Dual::log(DualNum<T>(1.0) - yhat[i][0]));
                        break;
                    case lossFunction::categorical_crossentropy:
                        result = result - Dual::log(yhat[i][static_cast<int>(y_train[i].getReal())]);
                        break;
                }
            }

            return result / static_cast<T>(yhat.size());
        };

        int no_of_threads = std::thread::hardware_concurrency();

        for (int epoch = 1; epoch <= epochs; epoch++) {
            vector<DualNum<T>> newParameters(parameters.size(), DualNum<T>(0, 0));
            vector<std::thread> threads;

            if (no_of_threads >= parameters.size()) {
                for (size_t i = 0; i < parameters.size(); i++) {
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

            T loss_value = 0;
            switch (lossFunc) {
                case lossFunction::mean_squared_error:
                    loss_value = Dual::mse(y_train, this->predictAll(X_train)).getReal();
                    break;
                case lossFunction::binary_crossentropy:
                    loss_value = Dual::binary_crossentropy(y_train, this->predictAll(X_train)).getReal();
                    break;
                case lossFunction::categorical_crossentropy:
                    loss_value = Dual::categorical_crossentropy(y_train, this->predictAll(X_train)).getReal();
                    break;
            }

            if (verbose == 2) {
                std::cout << "\n Epoch no: " << epoch << " completed!, Loss = " << loss_value << std::endl;
            } else if (verbose == 1) {
                std::cout << "\n Epoch no: " << epoch << " completed!" << std::endl;
            } else if (verbose == 3) {
                vector<vector<DualNum<T>>> yhat = this->predictAll(X_train);
                T accuracy = Dual::accuracy(y_train, yhat).getReal();
                std::cout << "\n Epoch no: " << epoch << " completed!, Loss = " << loss_value << " accuracy = " << accuracy << std::endl;
            }
        }
    }
};

#endif
