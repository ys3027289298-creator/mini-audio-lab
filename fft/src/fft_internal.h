/* fft_internal.h — shared radix-2 complex FFT used by fft.c and spectral.c. */
#ifndef FFT_INTERNAL_H
#define FFT_INTERNAL_H

typedef struct { float re, im; } C;

void fft_forward(C *a, int n);
void fft_inverse(C *a, int n);

#endif
