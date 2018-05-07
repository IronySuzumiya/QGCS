#include <math.h>
#include <assert.h>

#include "Complex.h"
#include "Const.h"

Complex complex_rect(double real, double imag) {
    return (Complex) { real, imag };
}

Complex complex_polar(double mod, double arg) {
    return (Complex) { mod * cos(arg), mod * sin(arg) };
}

Complex complex_neg(Complex a) {
    return (Complex) { -REAL(a), -IMAG(a) };
}

Complex complex_conj(Complex a) {
    return (Complex) { REAL(a), -IMAG(a) };
}

Complex complex_add(Complex a, Complex b) {
    return (Complex) { REAL(a) + REAL(b), IMAG(a) + IMAG(b) };
}

Complex complex_sub(Complex a, Complex b) {
    return (Complex) { REAL(a) - REAL(b), IMAG(a) - IMAG(b) };
}

Complex complex_mul(Complex a, Complex b) {
    return (Complex) { REAL(a) * REAL(b) - IMAG(a) * IMAG(b) , REAL(a) * IMAG(b) + REAL(b) * IMAG(a) };
}

Complex complex_div(Complex a, Complex b) {
    assert(!complex_is_zero(b));
    Complex c = complex_div_real(complex_mul(a, complex_conj(b)), complex_norm(b));
    //assert(complex_equal(complex_mul(b, c), a));
    return c;
}

Complex complex_pow(Complex a, int b) {
    if (b == 0) {
        return (Complex) { 1, 0 };
    }
    else if (b == 1) {
        return a;
    }
    else if (b > 1) {
        Complex temp = a;
        for (int i = 0; i < b - 1; ++i) {
            temp = complex_mul(temp, a);
        }
        return temp;
    }
    else {
        Complex temp = complex_rect(1, 0);
        for (int i = 0; i < -b; ++i) {
            temp = complex_div(temp, a);
        }
        return temp;
    }
}

Complex complex_sqrt(Complex a) {
    double alpha = sqrt(complex_abs(a) + REAL(a));
    return (Complex) { alpha / sqrt(2), IMAG(a) / (sqrt(2) * alpha) };
}

Complex complex_mul_real(Complex a, double b) {
    return (Complex) { REAL(a) * b, IMAG(a) * b };
}

Complex complex_div_real(Complex a, double b) {
    //assert(!double_is_zero(b));
    return (Complex) { REAL(a) / b, IMAG(a) / b };
}

double complex_norm(Complex a) {
    return REAL(a) * REAL(a) + IMAG(a) * IMAG(a);
}

double complex_abs(Complex a) {
    return sqrt(complex_norm(a));
}

int complex_is_real(Complex a) {
    return IMAG(a) < EPSILON && IMAG(a) > -EPSILON;
}

int complex_is_imag(Complex a) {
    return REAL(a) < EPSILON && REAL(a) > -EPSILON;
}

int complex_is_zero(Complex a) {
    return complex_is_real(a) && complex_is_imag(a);
}

int complex_equal(Complex a, Complex b) {
    return complex_is_zero(complex_sub(a, b));
}

int double_is_zero(double a) {
    return a < EPSILON && a > -EPSILON;
}

int double_equal(double a, double b) {
    return double_is_zero(a - b);
}
