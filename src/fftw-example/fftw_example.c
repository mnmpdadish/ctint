/* Start reading here */

#include <fftw3.h>

#define NUM_POINTS 65536


/* Never mind this bit */

#include <stdio.h>
#include <math.h>

#define REAL 0
#define IMAG 1

void acquire_from_somewhere(fftw_complex* signal) {
    /* Generate two sine waves of different frequencies and
     * amplitudes.
     */
    
    int i;
    for (i = 0; i < NUM_POINTS; ++i) {
        double theta = (double)i /1000.0  ;

        signal[i][REAL] = 1.0 * cos(10.0 * theta) +
                          0.5 * cos(25.0 * theta);

        signal[i][IMAG] = 0.0 * sin(10.0 * theta) +
                          0.0 * sin(25.0 * theta);
    }
}

void do_something_with(fftw_complex* result) {
    int i;
    for (i = 0; i < NUM_POINTS; ++i) {
        double mag = sqrt(result[i][REAL] * result[i][REAL] +
                          result[i][IMAG] * result[i][IMAG]);

        printf("%g\n", mag);
    }
}


/* Resume reading here */

int main() {
    fftw_complex signal[NUM_POINTS];
    fftw_complex result[NUM_POINTS];

    fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
                                      signal,
                                      result,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE);

    acquire_from_somewhere(signal);
    fftw_execute(plan);
    do_something_with(result);

    fftw_destroy_plan(plan);
    printf("NUM_POINTS*sizeof(fftw_complex)= %d\n",NUM_POINTS*sizeof(fftw_complex));
    return 0;
}
