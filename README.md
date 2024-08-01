# Dual Numbers Library

This project provides a C++ library for dual numbers, including both standard and extended representations. Dual numbers are useful in various fields such as automatic differentiation, optimization, and computational geometry. The library supports basic arithmetic operations, standard mathematical functions, and matrix operations involving dual numbers.

## **Table of Contents**
- Installation
- Usage
- DualNum Class
- Standard Functions

## Installation
To use this library, include the header and source files in your C++ project. Make sure your build system compiles the dual.cpp file if it uses DualNum class and the extendedDual.cpp file if it uses the ExtendedDualNum class.

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

- `evaluatePartialDerivative(DualNum(*func)(DualNum), long double at = 1);`

    This function simply evaluates the given function using dual numbers set up for 1st derivative and returns the entire result as it is.

- `vector<DualNum> gradient(DualNum(*func)(vector<DualNum>), vector<DualNum> at, coordinate_system system = coordinate_system::cartesian);`
- `double solveUsingNewtonRaphson(DualNum(*func)(DualNum), double initialGuess = 1, int max_no_of_iterations = 10000);`


