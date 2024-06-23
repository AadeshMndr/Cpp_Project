#include<iostream>

class ExtendedDualNum {
    private:
        double real;
        double firstDual;
        double secondDual;

    public:
        ExtendedDualNum(double r = 0, double f = 0, double s = 0): real(r), firstDual(f), secondDual(s){
        }

        double getReal(){
            return real;
        }
        double getFirstDual(){
            return firstDual;
        }
        double getSecondDual(){
            return secondDual;
        }

        ExtendedDualNum operator + (ExtendedDualNum num) {
            return (ExtendedDualNum(real + num.real, firstDual + num.firstDual, secondDual + num.secondDual));
        }

        //Similarly overload all other operators
        // -, *, / , ^ (power), += , -=, *= , /=
        //do the above for all ExtendedDualNum * ExtendedDualNum, double * ExtendedDualNum (outside the class), ExtendedDualNum * double 
};

int main(){
    ExtendedDualNum n(4);

    std::cout << n.getReal();

    return 0;
}