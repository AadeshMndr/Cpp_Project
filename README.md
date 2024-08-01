# Dual Numbers Library

This project provides a C++ library for dual numbers, including both standard and extended representations. Dual numbers are useful in various fields such as automatic differentiation, optimization, and computational geometry. The library supports basic arithmetic operations, standard mathematical functions, and matrix operations involving dual numbers.

## **Table of Contents**
- Installation
- Usage
- DualNum Class
- Standard Functions
- ExtendedDualNum Class
- Compiling instructions

## Installation
To use this library, include the header and source files present inside the Dual_Library folder inside your C++ project. Make sure your build system compiles the dual.cpp file if it uses DualNum class and the extendedDual.cpp file if it uses the ExtendedDualNum class.

This library internally uses threads, which is compatible with gcc version 11 and after, so make sure you have a compiler that is equivalent to gcc version 11 or newer.

## Usage
Here is a basic example of how to use the DualNum class:

```cpp
#include "dual.h"

int main() {
    DualNum a(3.0, 1.0); // Create a dual number with real part 3.0 and dual part 1.0
    DualNum b(2.0, 0.0); // Create a dual number with real part 2.0 and dual part 0.0

    DualNum c = a * b;   // Perform multiplication

    std::cout << "Real part: " << c.getReal() << std::endl;
    std::cout << "Dual part: " << c.getDual() << std::endl;

    return 0;
}
```
### DualNum Class
The DualNum class represents dual numbers with a real part and a dual part.

#### Constructors

- `DualNum(long double r = 0, long double e = 0);`

  Creates a dual number with a real part r and a dual part e.
  
- `DualNum(const DualNum& num);`
  
  Copy constructor.
  
#### Member Functions
- `long double getReal() const;`

  Returns the real part of the dual number.
- `long double getDual() const;`
  
  Returns the dual part of the dual number.
  
- `void setReal(const long double& r);`
  
  Sets the real part of the dual number.
  
- `void setDual(const long double& e);`
  
  Sets the dual part of the dual number.
- `string getExpression() const;`
  
  Returns a string representation of the dual number.

  
#### Operator Overloads
The DualNum class supports standard arithmetic operations (+, -, *, /, ^) and comparison operators (==, >, <, >=, <=) for dual numbers and long double scalars.

#### Standard Functions
The library includes several standard functions that operate on dual numbers, such as:

- `Dual::pow(const DualNum& x, const DualNum& y);`
- `Dual::exp(const DualNum& x);`
- `Dual::log(const DualNum& x);`
- `Dual::relu(const DualNum& x);`
- `Dual::tanh(const DualNum& x);`
- `Dual::sigmoid(const DualNum& x);`
- `Dual::softmax(const vector<DualNum>& X, int index, DualNum sum = DualNum(0, 0));`
- `Dual::mse(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);`
- `Dual::binary_crossentropy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);`
- `Dual::categorical_crossentropy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat);`
- `Dual::accuracy(const vector<DualNum>& y_train, const vector<vector<DualNum>>& yhat, const double threshold = 0.5);`

*The above functions work exactly as the normally understood versions of it, but just compatible with dual numbers.*

#### Extra Functions
Additional functionalities include calculating partial derivatives, gradients, and solving equations using the Newton-Raphson method:

- `partialDerivative(DualNum(*func)(vector<DualNum>), vector<DualNum> params, int paramIndex = 0, long double at = 1);`
    
    This function returns the first derivative of the function passes evaluated at the value specified in 'at' parameter.

```cpp
#include <iostream>
#include <vector>
#include "dual.h"

// Define the function f(x, y) = x^2 + 3xy + y^2 using DualNum
DualNum function(std::vector<DualNum> params) {
    DualNum x = params[0];
    DualNum y = params[1];
    return (x * x) + (3 * x * y) + (y * y);
}

int main() {
    // Initialize the point (x, y) = (1, 2)
    std::vector<DualNum> params = { DualNum(1.0, 0.0), DualNum(2.0, 0.0) };

    // Calculate the partial derivative with respect to x (index 0) at the point (1, 2)
    long double partialDeriv_x = partialDerivative(function, params, 0);

    std::cout << "Partial derivative with respect to x at (1, 2): " << partialDeriv_x << std::endl;

    // Calculate the partial derivative with respect to y (index 1) at the point (1, 2)
    long double partialDeriv_y = partialDerivative(function, params, 1);

    std::cout << "Partial derivative with respect to y at (1, 2): " << partialDeriv_y << std::endl;

    return 0;
}
```

- `evaluatePartialDerivative(DualNum(*func)(DualNum), long double at = 1);`

    This function simply evaluates the given function using dual numbers set up for 1st derivative and returns the entire result as it is.

- `vector<DualNum> gradient(DualNum(*func)(vector<DualNum>), vector<DualNum> at, coordinate_system system = coordinate_system::cartesian);`
- `double solveUsingNewtonRaphson(DualNum(*func)(DualNum), double initialGuess = 1, int max_no_of_iterations = 10000);`

### ExtendedDualNum Class

The ExtendedDualNum class is designed for advanced numerical computation, particularly focusing on dual numbers and extended dual numbers. It uses dual numbers in it's matrix form and it can be used to calculate the derivative upto the nth term.

#### Constructors
- `ExtendedDualNum(double num, int orderOfMat = 2);`

    Initializes the dual number with a given real part and matrix order, the order of the matrix is the number of terms you want to have in your dual number.

- `ExtendedDualNum(vector<vector<double>> dataMatrix);`

    Initializes using a data matrix.

- `ExtendedDualNum(const ExtendedDualNum& another);`

    Copy constructor.


An example code of how we can use the ExtendedDualNum class.
```cpp
#include<iostream>
#include<vector>
#include "extendedDual.h"

using std::cout;
using std::endl;
using namespace ExtendedDual;

// Function to demonstrate usage
ExtendedDualNum sampleFunction(vector<ExtendedDualNum> params) {
    return params[0] * params[0] + params[1] * params[1];
}

int main() {
    // Initialize two ExtendedDualNum objects with a value and matrix size (order)
    ExtendedDualNum a(3.0, 2);
    ExtendedDualNum b(4.0, 2);

    // Perform arithmetic operations
    ExtendedDualNum sum = a + b;
    ExtendedDualNum difference = a - b;
    ExtendedDualNum product = a * b;
    ExtendedDualNum quotient = a / b;

    // Display results of arithmetic operations
    cout << "Sum: " << sum << endl;
    cout << "Difference: " << difference << endl;
    cout << "Product: " << product << endl;
    cout << "Quotient: " << quotient << endl;

    // Compute exponential, logarithm, and trigonometric functions
    ExtendedDualNum exponential = exp(a);
    ExtendedDualNum logarithm = log(a);
    ExtendedDualNum sine = sin(a);
    ExtendedDualNum cosine = cos(a);
    ExtendedDualNum tangent = tan(a);

    // Display results of functions
    cout << "Exponential: " << exponential << endl;
    cout << "Logarithm: " << logarithm << endl;
    cout << "Sine: " << sine << endl;
    cout << "Cosine: " << cosine << endl;
    cout << "Tangent: " << tangent << endl;

    return 0;
} 
```

#### Public Methods
- `vector<vector<double>> getData() const; `

    Returns the data matrix of the dual number.

- `double getNum() const;`

    Returns the real part of the dual number.

- `int getSize() const;`

    Returns the size of the data matrix.

- `vector<vector<double>> getDualMatrix() const;`

    Returns the matrix part of the dual number.

- `void setData(vector<vector<double>> newData);`

    Sets the data matrix.

- `void setDual(double dual = 1);`

    Sets the dual part of the number.

- `void displayData() const;`

    Displays the dual number data.

#### Operator Overloads
The ExtendedDualNum class supports standard arithmetic operations (+, -, *, /, ^) and comparison operators (==, >, <, >=, <=) for dual numbers and long double scalars.

#### Standard Functions
The library includes several standard functions that operate on dual numbers, exactly as the ones described for the DualNum class, the only difference being that now we have to use the ExtendedDual namespace instead of the Dual namespace for the standard functions.


#### Extra Functions
- `vector<double> partialDerivative(ExtendedDualNum (*func)(vector<ExtendedDualNum>);`
- `vector<ExtendedDualNum> params, int paramIndex = 0);`

    Computes the partial derivative of a vector function.

- `vector<double> partialDerivative(ExtendedDualNum (*func)(ExtendedDualNum), double at = 1, int order = 3);`

    Computes the partial derivative of a function.

- `ExtendedDualNum evaluatePartialDerivative(ExtendedDualNum(*func)(ExtendedDualNum), long double at = 1);`

    Evaluates the partial derivative at a given point.

An example on how we can use the partialDerivative functions:
```cpp
#include<iostream>
#include<iomanip>
#include"extendedDual.h"
#include<limits>

ExtendedDualNum someFunc(ExtendedDualNum x){
    return (ExtendedDual::log(x, x/9)) * (x^0) / ((x ^ (-1)) / 99);
}

int main() {

    std::cout << std::setprecision(10);

    int orderOfExtendedDualNumber = 9; //upto which derivative

    //You get an array (vector) that has the terms of the dual number
    vector<double> answers = ExtendedDual::partialDerivative(someFunc, 0.5, orderOfExtendedDualNumber);

    for (int i = 0; i < orderOfExtendedDualNumber; i++){
        std::cout << "\nThe " << i << " order derivative is: " << answers[i] << std::endl;
    }

    return 0;
}
```

#### Advanced Calculus
- `vector<double> gradient(VectorFunctionPointer, vector<ExtendedDualNum> at, coordinate_system system = coordinate_system::cartesian);`

    Computes the gradient of a function.

- `double laplacian(VectorFunctionPointer, vector<ExtendedDualNum> at, coordinate_system system = coordinate_system::cartesian);`

    Computes the Laplacian of a function.

- `vector<vector<double>> jacobian(vectorExtendedDual::VectorFunctionPointer, vector<ExtendedDualNum> at);`

    Computes the Jacobian matrix.

- `double solveUsingNewtonRaphson(ExtendedDualNum(*func)(ExtendedDualNum), double initialGuess = 1, int max_no_of_iterations = 10000);`

    Solves equations using the Newton-Raphson method.

An example of using these functions:
```cpp
#include<iostream>
#include "extendedDual.h"

ExtendedDualNum potential(vector<ExtendedDualNum> params){
    return (((params[0])^2) + params[1] + params[2]);
}

ExtendedDualNum fx(vector<ExtendedDualNum> params){
    return (2 * params[0] + params[1]);
}

ExtendedDualNum fy(vector<ExtendedDualNum> params){
    return (3 * params[0] + params[1]);
}

ExtendedDualNum someFunction(ExtendedDualNum x){
    return ((x ^ 2) + 5 * x - 50);
}

int main(){
    vector<double> gradientValue = ExtendedDual::gradient(potential, {1, 1, 1});

    std::cout << "\n\nf(x, y, z) = x^2 + y + z   at  (1, 1, 1)" << std::endl;
    std::cout << "\nThe gradient is: " << "(" << gradientValue[0] << ", " << gradientValue[1] << ", " << gradientValue[2] << ")" << std::endl;

    std::cout << "\nThe laplacian is: " << ExtendedDual::laplacian(potential, {ExtendedDualNum(1, 3), ExtendedDualNum(1, 3), ExtendedDualNum(1, 3)}) << std::endl;

    vector<vector<double>> jacobianMatrix = ExtendedDual::jacobian({fx, fy}, {1, 1});

    std::cout << "\n\nf(x, y) = (2x + y, 3x + y)   at  (1, 1)" << std::endl;
    std::cout << "\nThe Jacobian Matrix is: \n\n";
    for (int i = 0; i < jacobianMatrix.size(); i++){
        for (int j = 0; j < jacobianMatrix[0].size(); j++){
            std::cout << jacobianMatrix[i][j] << "\t";
        }
        std::cout << "\n\n";
    }

    double answer = ExtendedDual::solveUsingNewtonRaphson(someFunction);


    std::cout << "\nequation: x^2 + 5x - 50 = 0" << std::endl;
    std::cout << "\nThus the solution of the equation is " << answer << "\n\n" << std::endl;
    
    return 0;
}   
```

Using this library we can make differentiation logic for entire machine learning models as well, for demonstration purposes we have included an ANN model made using this library in this same repository, to test it, make sure you use a compatible version of the compiler as specified in the Usage section. 
Also make use you have SFML installed in your PC and it's path setup as specified in the make file, to view graphs, or else just comment out the parts that is related to generation of a graph. 

## Compiling Instructions

To use the exemplary applications that we have included in our repository / project:

- you need to make sure that you install SFML-2.6.1 or higher ->  https://www.sfml-dev.org/download/sfml/2.6.1
- then make sure to set the path for SFML in the make file as shown below (replace with <> parts with your paths):

```MakeFile
# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -g -I< Complete path to -> SFML-2.6.1/include >

# Linker flags
LDFLAGS = -L<Complete path to -> SFML-2.6.1/lib > -lsfml-graphics -lsfml-window -lsfml-system


.... the rest is the same as the make file as provided ....
```
- Then to run the program to compare automatic differentiation with numerical differentiation:
run `make run_DER` in terminal in the location of the project folder.

- Then to run the program that showcases evaluation of different mathematical operations:
run `make run_GRAD` in terminal in the location of the project folder.

- Then to run the program that showcases that ANN model can be trained for the XOR Dataset using this library:
run `make run_XOR` in terminal in the location of the project folder.

- Then to run the program that showcases that ANN model can be trained for regression using this library:
run `make run_REG` in terminal in the location of the project folder.

Happy Hacking !


