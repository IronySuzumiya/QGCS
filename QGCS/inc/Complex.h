#ifndef _COMPLEX_H
#define _COMPLEX_H

typedef struct __complex {
    float real;
    float imag;
} Complex;

#define REAL(x) ((x).real)
#define IMAG(x) ((x).imag)

Complex complex_rect(float real, float imag);
Complex complex_polar(float mod, float arg);
Complex complex_neg(Complex a);
Complex complex_conj(Complex a);
Complex complex_add(Complex a, Complex b);
Complex complex_sub(Complex a, Complex b);
Complex complex_mul(Complex a, Complex b);
Complex complex_div(Complex a, Complex b);
Complex complex_pow(Complex a, int b);
Complex complex_sqrt(Complex a);
Complex complex_mul_real(Complex a, float b);
Complex complex_div_real(Complex a, float b);
float complex_norm(Complex a);
float complex_abs(Complex a);
int complex_is_real(Complex a);
int complex_is_imag(Complex a);
int complex_is_zero(Complex a);
int complex_equal(Complex a, Complex b);
int float_is_zero(float a);
int float_equal(float a, float b);

#endif
