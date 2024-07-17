#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

/* Function prototypes */
int take_input(FILE *inputfp, char *filename, int *n, double **input);
void transform_array(double *input, int n);
void evaluateLevels(double *input, int n);
double cosN(int n, int l);
double sinN(int n, int l);
void evaluateHv(double *input, int n);
int print_output(FILE *outputfp, char *filename, int n, double *output);


int main(int argc, char **argv){
    // double input[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}; // Input array
    FILE *inputfp = NULL;       // Input file pointer
    FILE *outputfp = NULL;      // Output file pointer
    double *input = NULL;       // Input array

    int n = 0;                  // Size of the input array

    // Check if the input file is provided
    if (argc < 3 || argc > 3){
        printf("\nUsage: %s <input filename> <output filename>", argv[0]);
        return 1;
    }

    if (take_input(inputfp, argv[1], &n, &input)){
        return 1;
    }

    // start the clock
    clock_t t = clock();

    // create level 0: f_0 (tau)
    transform_array(input, n);

    // create other levels
    evaluateLevels(input, n);

    // Evaluate H(v)
    evaluateHv(input, n);

    // stop the clock
    t = clock() - t;
    printf("\nSequential time: %f", (double) t / CLOCKS_PER_SEC);

    // Write the output to the file
    if (print_output(outputfp, argv[2], n, input)){
        return 1;
    }

    return 0;
}



/* Function definitions */


/* Read the input from the file
    * @param inputfp: file pointer
    * @param filename: name of the file
    * @param n: size of the input array
    * @param input: input array
    * @return 0 if the input is read successfully, 1 otherwise
*/
int take_input(FILE *inputfp, char *filename, int *n, double **input){
    // Open file for reading
    inputfp = fopen(filename, "r");
    if (inputfp == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the size of the array
    fscanf(inputfp, "%d", n);

    // Allocate memory for the input array
    *input = (double *)malloc(*n * sizeof(double));
    if (input == NULL) {
        perror("Memory allocation failed");
        return 1;
    }

    // Read the input array
    for (int i = 0; i < *n; i++) {
        fscanf(inputfp, "%lf", &(*input)[i]);
    }

    // Close the file
    fclose(inputfp);

    return 0;
}


/* Function to transform the input array into the tau level
    * @param input: input array
    * @param n: number of input points
    * @return void
*/
void transform_array(double *input, int n) {
    // base case
    if (n == 2) {
        return;
    }

    double *temp = (double*)malloc(n * sizeof(double));
    int mid = n / 2;

    // separate even and odd index elements
    for (int i = 0; i < mid; i++) {
        temp[i] = input[2 * i];
        temp[mid + i] = input[2 * i + 1];
    }

    for (int i = 0; i < mid; i++) {
        input[i] = temp[i];
        input[mid + i] = temp[mid + i];
    }

    // recursive call
    transform_array(input, mid);
    transform_array(input + mid, mid);

    free(temp);
}


/* Function to evaluate all the f(tau) levels
    * @param input: input array
    * @param n: size of the input array
    * @return void
*/
void evaluateLevels(double *input, int n){
    double *temp = (double*)malloc(n * sizeof(double));
    int logN = (int)log2(n);

    for (int s = 0; s < logN; s++){
        int m = pow(2, s);
        int u = m;

        for (int k = 0; k < n; k += 2*m){
            for (int j = 0; j < m; j++){
                double t  = input[k + j];        // first element
                double u1 = input[k + j + m];    // second element
                double u2 = input[u + k];        // third element

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

        // copy temp into input
        for (int i = 0; i < n; i++) {
            input[i] = temp[i];
        }
    }

    free(temp);
}


/* Function to calculate cos((2 * PI * n) / pow(2, l)) 
    * @param n: index
    * @param l: level
    * @return cos value
*/
double cosN(int n, int l){
    return cos((2 * PI * n) / pow(2, l));
}


/* Function to calculate sin((2 * PI * n) / pow(2, l)) 
    * @param n: index
    * @param l: level
    * @return sin value
*/
double sinN(int n, int l){
    return sin((2 * PI * n) / pow(2, l));
}


/* Function to evaluate H(v)
    * @param input: input array
    * @param n: size of the input array
    * @return void
*/
void evaluateHv(double *input, int n){
    for (int k = 0; k < n; k++){
        input[k] /= n;
    }
}


/* Write the output to the file
    * @param inputfp: file pointer
    * @param filename: name of the file
    * @param n: size of the output array
    * @param output: output array
    * @return 0 if the output is written successfully, 1 otherwise
*/
int print_output(FILE *outputfp, char *filename, int n, double *output){
    // Open file for writing
    outputfp = fopen(filename, "w");
    if (outputfp == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the input array
    for (int i = 0; i < n; i++) {
        fprintf(outputfp, "%f ", output[i]);
    }

    // Close the file
    fclose(outputfp);

    return 0;
}