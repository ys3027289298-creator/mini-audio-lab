/* recur: measure the float twiddle-recurrence error W(len) in isolation.
 * For each stage length len, run the repo's float recurrence
 *   w <- w * wl,  wl = (float)cos(2pi/len) + i*(float)sin(2pi/len)
 * for k = 0..len/2-1, against the exact (double) twiddle wl_exact^k.
 * Also reports |wl_float|-1 (systematic magnitude drift per step) and
 * whether cos(2pi/len) rounds to exactly 1.0f.
 *
 * Output CSV: len,Wmax,Wrms,magdrift_per_step,cos_is_one,drift_model
 * where drift_model = (len/2)*| |wl_fl|-1 | (linear-in-k systematic bound).
 */
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("len,Wmax,Wrms,magdrift_per_step,cos_is_one,drift_model\n");
    for (int p = 1; p <= 22; p++) {
        long len = 1L << p;
        double ang = 2.0 * 3.14159265358979 / len;
        float cfl = (float)cos(ang), sfl = (float)sin(ang);
        double cex = cos(ang), sex = sin(ang);
        double magfl = sqrt((double)cfl * cfl + (double)sfl * sfl);

        float wfr = 1.0f, wfi = 0.0f;
        double wer = 1.0, wei = 0;
        double wmax = 0.0, wsq = 0.0;
        long half = len / 2;
        for (long k = 0; k < half; k++) {
            double dr = (double)wfr - wer, di = (double)wfi - wei;
            double d = sqrt(dr * dr + di * di);
            if (d > wmax) wmax = d;
            wsq += d * d;
            float nr = wfr * cfl - wfi * sfl;
            float ni = wfr * sfl + wfi * cfl;
            wfr = nr; wfi = ni;
            double er = wer * cex - wei * sex;
            double ei = wer * sex + wei * cex;
            wer = er; wei = ei;
        }
        printf("%ld,%.9g,%.9g,%.9g,%d,%.9g\n",
               len, wmax, sqrt(wsq / half), magfl - 1.0,
               cfl == 1.0f ? 1 : 0, (half) * fabs(magfl - 1.0));
    }
    return 0;
}
