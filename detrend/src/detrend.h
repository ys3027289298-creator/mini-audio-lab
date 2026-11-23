/* detrend: 多项式去趋势(离散 Legendre 正交基投影) + 窗函数法 FIR 高通 + FFT */
#ifndef DETREND_H
#define DETREND_H

/* 从 x 中减去最小二乘拟合的 p 阶多项式,结果写入 out(可与 x 不同址) */
void detrend_poly(const double* x, int n, int p, double* out);

/* 窗函数法(Hamming)线性相位高通:截止 fc Hz,过渡带 trans Hz,采样率 fs。
   返回实际抽头数(奇数,不超过 maxtaps),系数写入 h。 */
int fir_highpass_design(double fc, double fs, double trans, double* h, int maxtaps);

/* 零相位(群延迟已补偿)FIR 卷积,边界外补零 */
void fir_apply_centered(const double* x, int n, const double* h, int m, double* y);

/* 基 2 迭代 FFT,n 须为 2 的幂;inverse=0 正变换 */
void fft_double(double* re, double* im, int n, int inverse);

#endif
