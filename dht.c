#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

/* Function prototypes */ 
int take_input(FILE *fp, char *filename, int *n, float **input);
void dht(int n, float *input, int norm, float *output);
void init(float *arr, int n);

/* Main function */
int main(int argc, char **argv){
    FILE *fp = NULL;        // Input file pointer
    int n = 0;              // Size of the input array
    float* input = NULL;    // Input array
    int norm = 1;           // Normalization factor

    // Check if the input file is provided
    if (argc < 2 || argc > 3){
        printf("\nUsage: %s <file> [norm]", argv[0]);
        return 1;
    }

    if (take_input(fp, argv[1], &n, &input)){
        return 1;
    }

    // print the input
    printf("\nInput: ");
    for (int i = 0; i < n; i++){
        printf("%.1f ", input[i]);
    }

    // if insert the denom of the normalization factor
    if (argc == 3 && atoi(argv[2]) != 0) {
        norm = n / atoi(argv[2]);
        printf("\nNormalization factor %d", norm);
    } else {
        norm = 1;
        printf("\nNormalization factor %d", norm);
    }

    // inizialized the output array
    float* output = (float *)malloc(n * sizeof(float));
    init(output, n);

    // call the dht function
    dht(n, input, norm, output);

    // print the output
    printf("\nOutput: ");
    for (int i = 0; i < n; i++){
        printf("%.1f ", output[i]);
    }

    return 0;
}


/* Function definitions */


/* Read the input from the file
    * @param fp: file pointer
    * @param filename: name of the file
    * @param n: size of the input array
    * @param input: input array
    * @return 0 if the input is read successfully, 1 otherwise
*/
int take_input(FILE *fp, char *filename, int *n, float **input){
    // Open file for reading
    fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the size of the array
    fscanf(fp, "%d", n);

    // Allocate memory for the input array
    *input = (float *)malloc(*n * sizeof(float));
    if (input == NULL) {
        perror("Memory allocation failed");
        return 1;
    }

    // Read the input array
    for (int i = 0; i < *n; i++) {
        fscanf(fp, "%f", &(*input)[i]);
    }

    // Close the file
    fclose(fp);

    return 0;
}


/* Calculate the Discrete Hartley Transform
    * @param n: size of the input array
    * @param input: input array
    * @param norm: normalization factor
    * @param output: output array
*/
void dht(int n, float *input, int norm, float *output){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            float a = 2 * PI * i * j / n;

            output[i] += (input[j] * (cos(a) + sin(a)) / n) * norm;
        }
    }
}


/* Initialize an array with zeros
    * @param arr: array to be initialized
    * @param n: size of the array
*/
void init(float *arr, int n){
    for (int i = 0; i < n; i++){
        arr[i] = 0;
    }
}