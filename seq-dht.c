#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

/* Function prototypes */ 
int take_input(FILE *inputfp, char *filename, int *n, float **input);
int print_output(FILE *outputfp, char *filename, int n, float *output);
void dht(int n, float *input, int norm, float *output);
void init(float *arr, int n);


/* Main function */
int main(int argc, char **argv){
    FILE *inputfp = NULL;   // Input file pointer
    FILE *outputfp = NULL;  // Output file pointer
    int n = 0;              // Size of the input array
    float* input = NULL;    // Input array
    int norm = 1;           // Normalization factor

    // Check if the input file is provided
    if (argc < 2 || argc > 4){
        printf("\nUsage: %s ./input/<input filename> ./output/<output filename> [norm]", argv[0]);
        return 1;
    }

    if (take_input(inputfp, argv[1], &n, &input)){
        return 1;
    }

    // print the input
    // printf("\nInput: ");
    // for (int i = 0; i < n; i++){
    //     printf("%.1f ", input[i]);
    // }

    // if insert the denom of the normalization factor
    if (argc == 4 && atoi(argv[3]) != 0) {
        norm = n / atoi(argv[3]);
    } else {
        norm = 1;
    }

    // inizialized the output array
    float* output = (float *)malloc(n * sizeof(float));
    init(output, n);

    // start the clock
    clock_t t = clock();

    // call the dht function
    dht(n, input, norm, output);

    t = clock() - t;

    printf("\nSequential time: %f", (double) t / CLOCKS_PER_SEC);

    if (print_output(outputfp, argv[2], n, output)){
        return 1;
    }

    // print the output
    // printf("\nOutput: ");
    // for (int i = 0; i < n; i++){
    //     printf("%.1f ", output[i]);
    // }

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
int take_input(FILE *inputfp, char *filename, int *n, float **input){
    // Open file for reading
    inputfp = fopen(filename, "r");
    if (inputfp == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the size of the array
    fscanf(inputfp, "%d", n);

    // Allocate memory for the input array
    *input = (float *)malloc(*n * sizeof(float));
    if (input == NULL) {
        perror("Memory allocation failed");
        return 1;
    }

    // Read the input array
    for (int i = 0; i < *n; i++) {
        fscanf(inputfp, "%f", &(*input)[i]);
    }

    // Close the file
    fclose(inputfp);

    return 0;
}


/* Write the output to the file
    * @param inputfp: file pointer
    * @param filename: name of the file
    * @param n: size of the output array
    * @param output: output array
    * @return 0 if the output is written successfully, 1 otherwise
*/
int print_output(FILE *outputfp, char *filename, int n, float *output){
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