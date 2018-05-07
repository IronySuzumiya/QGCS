#ifndef _COMPLEX_H
#define _COMPLEX_H

typedef struct __complex {
    double real;
    double imag;
} Complex;

#define REAL(x) ((x).real)
#define IMAG(x) ((x).imag)

Complex complex_rect(double real, double imag);
Complex complex_polar(double mod, double arg);
Complex complex_neg(Complex a);
Complex complex_conj(Complex a);
Complex complex_add(Complex a, Complex b);
Complex complex_sub(Complex a, Complex b);
Complex complex_mul(Complex a, Complex b);
Complex complex_div(Complex a, Complex b);
Complex complex_pow(Complex a, int b);
Complex complex_sqrt(Complex a);
Complex complex_mul_real(Complex a, double b);
Complex complex_div_real(Complex a, double b);
double complex_norm(Complex a);
double complex_abs(Complex a);
int complex_is_real(Complex a);
int complex_is_imag(Complex a);
int complex_is_zero(Complex a);
int complex_equal(Complex a, Complex b);
int double_is_zero(double a);
int double_equal(double a, double b);

#endif
