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

class ANN_Model {
private:
    vector<DualNum> parameters;
    int inputs;
    vector<activations> activationFunctions;
    vector<int> nodes;
    vector<double> velocity;
    double momentum = 0.1;

public:
    ANN_Model(int inputFeatures, vector<int> No_of_Nodes, vector<activations> activationFunctionName) {
        inputs = inputFeatures;

        nodes = No_of_Nodes;

        activationFunctions = activationFunctionName;

        if (activationFunctions.size() != nodes.size()) {
            throw std::runtime_error("The No. of Layers and the No. of activation functions don't match! ");
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

            parameters.push_back(DualNum(random_number));
        }

        velocity = vector<double>(parameters.size(), 0);
    }

    vector<DualNum>& getParameters(){   //returns by reference so cannot be made a const member function
        return parameters;
    }

    void setMomentum(double m) {
        momentum = m;
    }

    vector<DualNum> predict(const vector<DualNum>& X, vector<DualNum> params = {}) const {
        if (params.size() == 0) {
            params = parameters;
        }

        if (X.size() != inputs) {
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        vector<DualNum> inputVect = X;

        int start = 0;

        for (int i = 0; i < nodes.size(); i++) {
            vector<vector<DualNum>> currentWeightMatrix;


            for (int j = 0; j < inputVect.size(); j++) {
                int end = start + nodes[i];

                vector<DualNum> weightForEachNode(params.begin() + start, params.begin() + end);

                currentWeightMatrix.push_back(weightForEachNode);

                start = end;
            }

            vector<DualNum> currentBiasVect(params.begin() + start, params.begin() + start + nodes[i]);

            start += nodes[i];

            vector<DualNum> tempVect = Dual::originalMatMul(vector<vector<DualNum>>(1, inputVect), currentWeightMatrix)[0];

            DualNum expSums(0, 0);
            if (activationFunctions[i] == activations::softmax) {
                for (int j = 0; j < tempVect.size(); j++) {
                    expSums += Dual::exp(tempVect[j] + currentBiasVect[j]);
                }
            }

            for (int j = 0; j < tempVect.size(); j++) {
                tempVect[j] += currentBiasVect[j];

                if (activationFunctions[i] == activations::relu) {
                    tempVect[j] = Dual::relu(tempVect[j]);
                }
                else if (activationFunctions[i] == activations::sigmoid) {
                    tempVect[j] = Dual::sigmoid(tempVect[j]);
                }
                else if (activationFunctions[i] == activations::tanh) {
                    tempVect[j] = Dual::tanh(tempVect[j]);
                }
                else if (activationFunctions[i] == activations::softmax) {
                    tempVect[j] = Dual::softmax(tempVect, j, expSums);
                }
            }

            inputVect = tempVect;
        }

        return inputVect;
    }

    vector<vector<DualNum>> predictAll(const vector<vector<DualNum>>& data, vector<DualNum> params = {}) {
        if (data[0].size() != inputs) {
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        vector<vector<DualNum>> entireOutputs;

        for (int i = 0; i < data.size(); i++) {
            entireOutputs.push_back(predict(data[i], params));
        }

        return entireOutputs;
    }

    template <typename Func>
    void updateParameters(double learning_rate, Func loss, vector<DualNum>& newParameters, int start, int end) {
        for (int i = start; i < end; i++) {
            double gradient = partialDerivative(loss, parameters, i, parameters[i].getReal());
            velocity[i] = momentum * velocity[i] - learning_rate * gradient;
            newParameters[i] = parameters[i] + velocity[i];
        }
    }

    void fit(const vector<vector<DualNum>>& X_train, const vector<DualNum>& y_train, lossFunction lossFunc = lossFunction::mean_squared_error, double learning_rate = 0.01, int epochs = 100, int verbose = 2) {
        if (y_train.size() != X_train.size()) {
            throw std::runtime_error("The shape of the X_train and y_train don't match !");
        }

        auto loss = [X_train, y_train, this, lossFunc](vector<DualNum> params) mutable {
            vector<vector<DualNum>> yhat = this->predictAll(X_train, params);

            DualNum result(0, 0);

            for (int i = 0; i < yhat.size(); i++) {
                if (lossFunc == lossFunction::mean_squared_error) {
                    result += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
                }
                else if (lossFunc == lossFunction::binary_crossentropy) {
                    result -= ((y_train[i] * Dual::log(yhat[i][0]) + (1 - y_train[i]) * Dual::log(1 - yhat[i][0])));
                }
                else if (lossFunc == lossFunction::categorical_crossentropy) {
                    result -= Dual::log(yhat[i][y_train[i].getReal()]);
                }
            }

            return result / (double)(yhat.size());
            };

        int no_of_threads = std::thread::hardware_concurrency();

        for (int epoch = 1; epoch <= epochs; epoch++) {
            vector<DualNum> newParameters(parameters.size(), DualNum(0, 0));

            vector<std::thread> threads;

            if (no_of_threads >= parameters.size()) {
                for (int i = 0; i < parameters.size(); i++) {
                    threads.emplace_back(&ANN_Model::updateParameters<decltype(loss)>, this, learning_rate, loss, std::ref(newParameters), i, i + 1);
                }
            }
            else {
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

            double loss_value = 0;
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
            }
            else if (verbose == 1) {
                std::cout << "\n Epoch no: " << epoch << " completed!" << std::endl;
            }
            else if (verbose == 3) {
                vector<vector<DualNum>> yhat = this->predictAll(X_train);
                std::cout << "\n Epoch no: " << epoch << " completed!, Loss = " << loss_value << " accuracy = " << Dual::accuracy(y_train, yhat).getReal() << std::endl;
            }

        }
    }
};

#endif