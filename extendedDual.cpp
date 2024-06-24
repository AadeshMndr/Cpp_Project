#include<stdio.h>
#include<vector>
#include<cmath>
#include<string>

using std::vector;

class extendedDual {
    private:
        long double real;
        long double dual1;
        long double dual2;

    public:

        extendedDual(long double r, long double d1, long double d2) {
            real = r;
            dual1 = d1;
            dual2 = d2;
        }

        extendedDual(const extendedDual& ed) {
            real = ed.real;
            dual1 = ed.dual1;
            dual2 = ed.dual2;
        }

        long double getReal() {
            return real;
        }

        long double getDual1() {
            return dual1;
        }

        long double getDual2() {
            return dual2;
        }

        void setReal(long double r) {
            real = r;
        }

        void setDual1(long double d1) {
            dual1 = d1;
        }

        void setDual2(long double d2) {
            dual2 = d2;
        }

        extendedDual operator=(const extendedDual& ed) {
            real = ed.real;
            dual1 = ed.dual1;
            dual2 = ed.dual2;
            return *this;
        }

        extendedDual operator+(const extendedDual& ed) {
            extendedDual result(real + ed.real, dual1 + ed.dual1, dual2 + ed.dual2);
            return result;
        }

        extendedDual operator-(const extendedDual& ed) {
            extendedDual result(real - ed.real, dual1 - ed.dual1, dual2 - ed.dual2);
            return result;
        }
        
        extendedDual operator*(const extendedDual& ed) {
            extendedDual result(real * ed.real, real * ed.dual1 + dual1 * ed.real, real*ed.dual2 + 2 * dual1 * ed.dual1 + dual2 * ed.real);
            return result;
        }




};