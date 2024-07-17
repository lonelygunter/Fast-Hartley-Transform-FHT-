#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

/* Function prototypes */
int take_input(FILE *inputfp, char *filename, int *n, float **input);

int main(int argc, char **argv){
    FILE *inputfp = NULL;   // Input file pointer
    FILE *outputfp = NULL;  // Output file pointer
    int n = 0;              // Size of the input array
    float* input = NULL;    // Input array
    int norm = 1;           // Normalization factor

    /* MPI variables */
    int id;                 // Rank of the process
    int p;                  // Number of processes
    double elapsed_time;    // Execution time
    float h = 0;            // H[i]

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

    // output array


    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Start timer
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    // Get the rank and number of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // Calculate H[i]
    for (int j = 0; j < n; j++){
        float a = 2 * PI * id * j / n;

        h += (input[j] * (cos(a) + sin(a)) / n) * norm;
    }

    // Stop timer
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    // print output
    printf ("Process %d H[i] = %f\n", id, h);
    fflush (stdout);

    // Print the execution time of the process 0
    if (!id) {
        printf ("Execution time %8.6f\n", elapsed_time);
        fflush (stdout);
    }

    MPI_Finalize();

    // print output
    

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


