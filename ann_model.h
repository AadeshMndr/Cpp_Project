#ifndef ANN_MODEL
#define ANN_MODEL

#include<iostream>
#include<vector>
#include<random>
#include<stdexcept>
#include<thread>
#include "dual.h"
// #include "rough.h"

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
    vector<DualNum> parameters;                     // to store all parameters of the network
    int inputs;                                     // to store no of inputs
    vector<activations> activationFunctions;        // to store the "enums of the activation functions" applied to each layer 
    vector<int> nodes;                              // to store the "no of nodes" in each layer
    vector<double> velocity;                        // to store the velocity value for each parameter
    double momentum = 0.1;
    vector<double> lossHistory;
    vector<double> accuracyHistory;
    vector<double> lossHistory_val;
    vector<double> accuracyHistory_val;

public:
    ANN_Model(int inputFeatures, vector<int> No_of_Nodes, vector<activations> activationFunctionName) {
        inputs = inputFeatures;

        nodes = No_of_Nodes;

        activationFunctions = activationFunctionName;

        if (activationFunctions.size() != nodes.size()) {
            throw std::runtime_error("The No. of Layers and the No. of activation functions don't match! ");
        }

        //Some prerequiste initialization that we have to do to generate random numbers.
        // The generator / chooser
        std::random_device rd;
        std::mt19937 gen(rd());
        // The range / area for choosing
        std::uniform_real_distribution<> distrib(-0.5, 0.5);   //within this range  (a normal distribution)


        //Calculation of no of total parameters required according to the network prescribed.
        long int no_of_parameters = inputFeatures * No_of_Nodes[0] + No_of_Nodes[0];
        for (int i = 0; i < No_of_Nodes.size() - 1; i++) {
            no_of_parameters += No_of_Nodes[i] * No_of_Nodes[i + 1] + No_of_Nodes[i + 1];
        }


        //Initializing the parameters with random numbers from normal distribution. (In a linear form)
        for (int i = 0; i < no_of_parameters; i++) {
            double random_number = distrib(gen);

            parameters.push_back(DualNum(random_number));
        }

        //Initializing the velocity parameter corresponding to each parameter with 0.
        velocity = vector<double>(parameters.size(), 0);
    }

    vector<double> getLossHistory() {
        return lossHistory;
    }

    vector<double> getAccuracyHistory() {
        return accuracyHistory;
    }

    vector<double> getValLossHistory() {
        return lossHistory_val;
    }

    vector<double> getValAccuracyHistory() {
        return accuracyHistory_val;
    }

    vector<DualNum>& getParameters() {   //returns by reference so cannot be made a const member function
        return parameters;
    }

    void setMomentum(double m) {
        momentum = m;
    }

    // X is supposed to be in the form:
    /*        eg:  for input = 3
        (feature_1, feature_2 feature_3)   -> A single input vector.
    */
    vector<DualNum> predict(const vector<DualNum>& X, vector<DualNum> params = {}) const {

        // If it is specified some parameters then it works only with those parameters,
        // else it works with all parameters available
        if (params.size() == 0) {
            params = parameters;
        }

        if (X.size() != inputs) {      //inputs = No of input features
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        vector<DualNum> inputVect = X;


        //Now we try to send this inputVect through the entire network.

        int start = 0; // pointer / reading head for the parameter list

        // For each layer we do some calculation (This runs "no of layers" of times)
        for (int i = 0; i < nodes.size(); i++) {
            vector<vector<DualNum>> currentWeightMatrix;

            //we generate the weight matrix from the parameters stored in linear form.
            /*
                    eg: for input = 3, no of nodes in this layer = 5
                [
                    (param, param, param, param, param),
                    (param, param, param, param, param),
                    (param, param, param, param, param),
                ]
            */
            for (int j = 0; j < inputVect.size(); j++) {
                int end = start + nodes[i];

                vector<DualNum> weightForEachNode(params.begin() + start, params.begin() + end);

                currentWeightMatrix.push_back(weightForEachNode);

                start = end;
            }

            //we generate a bias vector
            vector<DualNum> currentBiasVect(params.begin() + start, params.begin() + start + nodes[i]);

            //to let the reading head know that the bias vector was initialized.
            start += nodes[i];


            // This requires (1 x inputs) and (inputs x no of nodes) shapes to get ---> (1 x no of nodes)
            vector<DualNum> tempVect = Dual::originalMatMul(vector<vector<DualNum>>(1, inputVect), currentWeightMatrix)[0];


            //just runs if activation function is softmax (for optimization purposes)
            DualNum expSums(0, 0);
            if (activationFunctions[i] == activations::softmax) {
                for (int j = 0; j < tempVect.size(); j++) {
                    expSums += Dual::exp(tempVect[j] + currentBiasVect[j]);
                }
            }


            // Adding the bias term and applying the activation function.
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


            // Now after the forward pass through this "single" layer is done, we provide the result of this
            // layer, (ie tempVect) as the inputVect for the next layer.
            inputVect = tempVect;
        }

        //finally we return the result after the forward pass of this "single" Input vector through the entire network.
        return inputVect;
    }


    /*
        The data matrix is like:

            eg: for inputs = 3, and no of data = 6.

        [
            (feature_1, feature_2, feature_3), --> a single input vector
            (feature_1, feature_2, feature_3), --> a single input vector
            (feature_1, feature_2, feature_3), --> a single input vector
            (feature_1, feature_2, feature_3), --> a single input vector
            (feature_1, feature_2, feature_3), --> a single input vector
            (feature_1, feature_2, feature_3), --> a single input vector
        ]
    */
    // This function runs predict for all data. (all input vectors in data)
    vector<vector<DualNum>> predictAll(const vector<vector<DualNum>>& data, vector<DualNum> params = {}) {
        if (data[0].size() != inputs) {
            throw std::runtime_error("The shape of the inputs isn't right!");
        }

        //is to store the final output of all data.
        vector<vector<DualNum>> entireOutputs;

        //runs for each input vector, (ie, runs "no of input vectors present in data" number of times)
        for (int i = 0; i < data.size(); i++) {
            entireOutputs.push_back(predict(data[i], params));
        }

        //Mind that each output is also in a vector form, even if that output is a vector of size 1.
        // so it will look like:
        /*
                eg: for data = 6

            [
                output vector,
                output vector,
                output vector,
                output vector,
                output vector,
                output vector,
            ]
        */

        //the size of the output vector = no of nodes in the last layer

        return entireOutputs;
    }


    // The temp Func was added so that this function can later be passed to .emplace_back in order to work with threads
    // This function is able to run gradient descent for specified range of parameters (eg: parameter 4 to parameter 12) and place the updated parameter in the storage container provided (ie: newParameters)
    template <typename Func>
    void updateParameters(double learning_rate, Func loss, vector<DualNum>& newParameters, int start, int end) {
        for (int i = start; i < end; i++) {

            //Gradient descent algorithm with momentum and velocity.
            double gradient = partialDerivative(loss, parameters, i, parameters[i].getReal());   //<---- Only here is the loss function actually executed so that we get f(x) along with f'(x)
            velocity[i] = momentum * velocity[i] - learning_rate * gradient;
            newParameters[i] = parameters[i] + velocity[i];
        }
    }

    //Here we run the epochs. (ie: predict, calculate loss and get parameter updates, actually make those updates, calculate loss / accuracy for display purposes, and do it all over again)
    void fit(const vector<vector<DualNum>>& X_train, const vector<DualNum>& y_train, const vector<vector<DualNum>>& X_val = {}, const vector<DualNum>& y_val = {}, lossFunction lossFunc = lossFunction::mean_squared_error, double learning_rate = 0.01, int epochs = 100, int verbose = 2) {
        if (y_train.size() != X_train.size()) {
            throw std::runtime_error("The shape of the X_train and y_train don't match !");
        }

        // The loss function as a lambda function
        auto loss = [X_train, y_train, this, lossFunc](vector<DualNum> params) {

            //Predict for all data.
            vector<vector<DualNum>> yhat = this->predictAll(X_train, params);

            DualNum result(0, 0); // ---> holds loss value

            // Calculate loss according to the predicted results
            for (int i = 0; i < yhat.size(); i++) {
                if (lossFunc == lossFunction::mean_squared_error) {
                    result += (y_train[i] - yhat[i][0]) * (y_train[i] - yhat[i][0]);
                }
                else if (lossFunc == lossFunction::binary_crossentropy) {
                    result -= ((y_train[i] * Dual::log(yhat[i][0]) + (1.0 - y_train[i]) * Dual::log(1.0 - yhat[i][0])));
                }
                else if (lossFunc == lossFunction::categorical_crossentropy) {
                    result -= Dual::log(yhat[i][y_train[i].getReal()]);
                }
            }

            return result / (double)(yhat.size());
            };


        // check no of threads
        int no_of_threads = std::thread::hardware_concurrency();

        accuracyHistory.clear();
        lossHistory.clear();
        lossHistory_val.clear();
        accuracyHistory_val.clear();

        for (int epoch = 1; epoch <= epochs; epoch++) {
            //Initialize an empty parameter list container, to later store the updated values
            vector<DualNum> newParameters(parameters.size(), DualNum(0, 0));

            //Initialize a container for the threads
            vector<std::thread> threads;

            if (no_of_threads >= parameters.size()) {
                // if we have more threads than parameters,
                // then make each thread deal with "one single" parameter
                // ie: we make each thread run the "updateParameters" function with the loss function we just defined as an argument, for a "single parameter".

                for (int i = 0; i < parameters.size(); i++) {
                    threads.emplace_back(&ANN_Model::updateParameters<decltype(loss)>, this, learning_rate, loss, std::ref(newParameters), i, i + 1);
                }
            }
            else {
                // If we have more no of parameters then threads then,
                // We make each thread deal with a portion of the parameters, (equally divided, with the last thread getting any extra remaining ones)

                int start = 0;
                int portion_for_each = parameters.size() / no_of_threads;

                for (int i = 0; i < no_of_threads; i++) {
                    int end = (i == (no_of_threads - 1) ? parameters.size() : start + portion_for_each);
                    threads.emplace_back(&ANN_Model::updateParameters<decltype(loss)>, this, learning_rate, loss, std::ref(newParameters), start, end);
                    start = end;
                }
            }


            // We wait for each thread to complete it's job.
            for (auto& thread : threads) {
                thread.join();
            }

            //actually update the parameters with it's new values
            parameters = newParameters;


            // calculate the loss again for displaying purposes
            double loss_value = 0;

            double loss_value_val = 0;

            vector<vector<DualNum>> yhat = this->predictAll(X_train);
            vector<vector<DualNum>> yhat_val;

            if (X_val.size() != 0) {
                yhat_val = this->predictAll(X_val);
            }

            switch (lossFunc) {
            case lossFunction::mean_squared_error:
                loss_value = Dual::mse(y_train, yhat).getReal();
                if (X_val.size() != 0) {
                    loss_value_val = Dual::mse(y_val, yhat_val).getReal();
                }
                break;
            case lossFunction::binary_crossentropy:
                loss_value = Dual::binary_crossentropy(y_train, yhat).getReal();
                if (X_val.size() != 0) {
                    loss_value_val = Dual::binary_crossentropy(y_val, yhat_val).getReal();
                }
                break;
            case lossFunction::categorical_crossentropy:
                loss_value = Dual::categorical_crossentropy(y_train, yhat).getReal();
                if (X_val.size() != 0) {
                    loss_value_val = Dual::categorical_crossentropy(y_val, yhat_val).getReal();
                }
                break;
            }

            lossHistory.push_back(loss_value);
            if (X_val.size() != 0) {
                lossHistory_val.push_back(loss_value_val);
            }


            double accuracyValue;
            double accuracyValue_val;
            if (lossFunc != lossFunction::mean_squared_error) {
                accuracyValue = Dual::accuracy(y_train, yhat).getReal();
                accuracyHistory.push_back(accuracyValue);
                if (X_val.size() != 0) {
                    accuracyValue_val = Dual::accuracy(y_val, yhat_val).getReal();

                    accuracyHistory_val.push_back(accuracyValue_val);
                }
            }

            if (verbose == 2) {
                if (X_val.size() != 0 ){
                    std::cout << "\n Epoch no: " << epoch << " completed! \n Training Loss = " << loss_value << " Validation Loss = " << loss_value_val << "\n" << std::endl;
                } else {
                    std::cout << "\n Epoch no: " << epoch << " completed!, Loss = " << loss_value << std::endl;
                }
            }
            else if (verbose == 1) {
                std::cout << "\n Epoch no: " << epoch << " completed!" << std::endl;
            }
            else if (verbose == 3) {

                if (X_val.size() != 0) {
                    std::cout << "\n Epoch no: " << epoch << " completed! \n Training Loss = " << loss_value << " Validation Loss = " << loss_value_val << "\n Training accuracy = " << accuracyValue << " Validation accuracy = " << accuracyValue_val << "\n" << std::endl;
                }
                else {
                    std::cout << "\n Epoch no: " << epoch << " completed! \n Loss = " << loss_value << " accuracy = " << accuracyValue << std::endl;
                }

            }

        }
    }
};

#endif