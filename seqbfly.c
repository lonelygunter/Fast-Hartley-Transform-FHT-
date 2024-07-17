#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

/* Function prototypes */
double cosN(int n, int l);
double sinN(int n, int l);
void transform_array(double *data, int n);


int main(int argc, char **argv){
    double data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}; // Data array
    int n = sizeof(data)/sizeof(data[0]); // lenght of data array
    int logN = (int)log2(n);

    // create level 0: f_0 (tau)
    transform_array(data, n);

    // print level 0
    printf("\nf0:\t");
    for (int i = 0; i < n; i++){
        printf("%.2f\t", data[i]);
    }


    // create other levels
    double *temp = (double*)malloc(n * sizeof(double));

    for (int s = 0; s < logN; s++){
        int m = pow(2, s);
        int u = m;

        for (int k = 0; k < n; k += 2*m){
            for (int j = 0; j < m; j++){
                double t  = data[k + j];        // first element
                double u1 = data[k + j + m];    // second element
                double u2 = data[u + k];        // third element

                // butterfly operations
                temp[k + j] = t + (u1 * cosN(j, s+1)) + (u2 * sinN(j, s+1));
                temp[k + j + m] = t + (u1 * cosN(j+m, s+1)) + (u2 * sinN(j+m, s+1));

                // update u
                if (u == m) {
                    u = pow(2, s+1)-1;
                } else {
                    u--;
                }
            }
        }

        // copy temp into data
        for (int i = 0; i < n; i++) {
            data[i] = temp[i];
        }

        // print level data
        printf("\nf%d:\t", s+1);
        for (int i = 0; i < n; i++){
            printf("%.2f\t", data[i]);
        }
    }

    free(temp);


    // H(v)
    for (int k = 0; k < n; k++){
        data[k] /= n;
    }

    // print H(v)
    printf("\n\nH(v):\t");
    for (int i = 0; i < n; i++){
        printf("%.2f\t", data[i]);
    }

    return 0;
}



/* Function definitions */


/* Function to calculate cos((2 * PI * n) / pow(2, l)) 
    n: index
    l: level
*/
double cosN(int n, int l){
    return cos((2 * PI * n) / pow(2, l));
}


/* Function to calculate sin((2 * PI * n) / pow(2, l)) 
    n: index
    l: level
*/
double sinN(int n, int l){
    return sin((2 * PI * n) / pow(2, l));
}


/* Function to transform array
    data: data array
    n: number of data points
*/
void transform_array(double *data, int n) {
    // base case
    if (n == 2) {
        return;
    }

    double *temp = (double*)malloc(n * sizeof(double));
    int mid = n / 2;

    // separate even and odd index elements
    for (int i = 0; i < mid; i++) {
        temp[i] = data[2 * i];
        temp[mid + i] = data[2 * i + 1];
    }

    for (int i = 0; i < mid; i++) {
        data[i] = temp[i];
        data[mid + i] = temp[mid + i];
    }

    // recursive call
    transform_array(data, mid);
    transform_array(data + mid, mid);

    free(temp);
}