#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

// Bit reversal function
int bit_reversal(int p, int j) {
    int j1 = j;
    int k0 = 0;
    for (int q = 1; q <= p; q++) {
        int j2 = j1 / 2;
        k0 = k0 * 2 + (j1 - 2 * j2);
        j1 = j2;
    }
    return k0;
}

// FFT function
void fft(int p, double *real, double *img) {
    int n = 1 << p; // n = 2^p
    int n2 = n / 2;
    int f1 = p - 1;
    int k = 0;

    // Main FFT loop
    for (int L = 1; L <= p; L++) {
        for (int i = 1; i <= n2; i++) {
            int j = bit_reversal(p, k) / (1 << f1);
            double a = 2 * PI * j / n;
            double c = cos(a);
            double s = sin(a);
            double u = real[k + n2] * c + img[k + n2] * s;
            double v = img[k + n2] * c - real[k + n2] * s;
            real[k + n2] = real[k] - u;
            img[k + n2] = img[k] + v;
            k++;
        }
        k += n2;
        if (k < n) continue;
        k = 0;
        f1--;
        n2 /= 2;
    }

    // Bit reversal for output
    for (k = 0; k < n; k++) {
        int j = bit_reversal(p, k);
        if (j > k) {
            double u = real[k] / n;
            double v = img[k] / n;
            real[k] = real[j] / n;
            img[k] = img[j] / n;
            real[j] = u;
            img[j] = v;
        }
    }
}

int main(int argc, char **argv) {
    int p = 3; // Example power of 2
    int n = 1 << p; // n = 2^p

    double real[n]; // Array to store real parts
    double img[n]; // Array to store imaginary parts

    // Initialize real and img
    for (int i = 0; i < n; i++) {
        real[i] = i + 1; // Real part initialized to 1, 2, ..., n
        img[i] = 0;     // Imaginary part initialized to 0
    }

    // Perform FFT
    fft(p, real, img);

    // Print results
    printf("Results:\n");
    for (int i = 0; i < n; i++) {
        printf("real[%d] = %f, img[%d] = %f\n", i, real[i], i, img[i]);
    }

    return 0;
}
