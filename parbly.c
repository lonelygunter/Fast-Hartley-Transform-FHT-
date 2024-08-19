#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846

/* Function prototypes */
int take_input(FILE *inputfp, char *filename, int *n, double **input);
void transform_array(double *input, int n);
int eval_ebp(int n, int p, int level);
int aggregate_data(double *input, int p, int rank, int ebp);
void eval_levels(double *input, int n, int rank, int size);
double cosN(int n, int l);
double sinN(int n, int l);
void eval_Hv(double *input, int n);
int print_output(FILE *outputfp, char *filename, int n, double *output);


int main(int argc, char **argv) {
    FILE *inputfp = NULL;       // Input file pointer
    FILE *outputfp = NULL;      // Output file pointer
    double *input = NULL;       // Input array

    int n = 0;                  // Size of the input array

    // Check if the input file is provided
    if (argc < 3 || argc > 3){
        printf("\nUsage: %s <input filename> <output filename>", argv[0]);
        return 1;
    }

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank and size
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("\nrank: %d", rank);

    if (!rank) {
        if (take_input(inputfp, argv[1], &n, &input)){
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // print the input array
        printf("\n%d; ar: ", rank);
        for (int i = 0; i < n; i++) {
            printf("%.2f\t", input[i]);
        }
        printf("\n");

        fflush(stdout);
    }

    // Broadcast the size of the input array to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for the input array in non-root processes
    if (rank) {
        input = (double *)malloc(n * sizeof(double));

        if (input == NULL) {
            perror("Memory allocation failed");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
    }

    // start the clock
    MPI_Barrier(MPI_COMM_WORLD);
    double time = -MPI_Wtime();

    // Aggregate the data from the root process
    int ebp = 0;
    ebp = eval_ebp(n, size, 1);
    printf("\nebp: %d", ebp);

    if (!rank) {
        // create level 0: f_0 (tau)
        transform_array(input, n);

        if (ebp == 1){
            // TODO
        }

        // Aggregate the data
        if (aggregate_data(input, size, rank, ebp)){
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        fflush(stdout);
    } else {
        // Receive the data from the root process
        MPI_Recv(input, ebp, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    // print level 0
    printf("\n%d; f0:\t", rank);
    for (int i = 0; i < n; i++) {
        printf("%.2f\t", input[i]);
    }

    // create other levels
    eval_levels(input, n, rank, size);
    
    // aggregate the data for H(v)


    // Evaluate H(v)
    if(!rank){
        eval_Hv(input, n);

        // print H(v)
        printf("\n%d; Hv:\t", rank);
        for (int i = 0; i < n; i++) {
            printf("%.2f\t", input[i]);
        }
        printf("\n");
    }

    // stop the clock
    MPI_Barrier(MPI_COMM_WORLD);
    time += MPI_Wtime();

    if (rank == 0){
        printf("\nSequential time: %f", time);

        printf("\n----3-------");
    }

    MPI_Finalize();

    // Write the output to the file
    if (rank == 0) {
        if (print_output(outputfp, argv[2], n, input)){
            return 1;
        }

        printf("\n----4-------");
    }

    free(input);

    return 0;
}



/* Function definitions */


/** Read the input from the file
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
    if (*input == NULL) {
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


/** Function to transform the input array into the tau level
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


/** Function to aggregate the data
    * @param n: size of the input array
    * @param p: number of processes
    * @param level: current level
    * @return element by process
*/
int eval_ebp(int n, int p, int level){
    int p_ideal = 0;

    p_ideal = n / pow(2, level);

    printf("\np: %d", p);
    printf("\np_ideal: %d = %d / 2^%d", p_ideal, n, level);

    if (!(p <= p_ideal)){
        return 1;
    }

    printf("\nebp = (%d/%d) * 2^%d", p_ideal, p, level);

    return (p_ideal/p) * pow(2, level);
}


/** Function to aggregate the data
    * @param input: input array
    * @param p: number of processes
    * @param ebp: element by process
    * @return 0 if the data is aggregated successfully, 1 otherwise
*/
int aggregate_data(double *input, int p, int rank, int ebp){
    for (int i = 0; i < p; i += ebp){

        for (int j = 1; j < ebp; j++){
            if (MPI_Send(&input[i+j-1+ebp], ebp, MPI_DOUBLE, j, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
                return 1;
            }
        }
    }

    return 0;
}


/** Function to evaluate all the f(tau) levels
    * @param input: input array
    * @param n: size of the input array
    * @param rank: rank of the current process
    * @param size: total number of processes
    * @return void
*/
void eval_levels(double *input, int n, int rank, int size){
    double *temp = (double*)malloc(n * sizeof(double));
    int logN = (int)log2(n);

    for (int s = 0; s < logN; s++){
        int m = pow(2, s);
        int u = m;
        int j = rank % m;

        printf("\n%d; ev: ", rank);
        for (int i = 0; i < n; i++) {
            printf("%.2f\t", input[i]);
        }

        // Parallelize this loop

        if (j + m < n) {
            double t  = input[j];        // first element
            double u1 = input[j + m];    // second element
            double u2 = input[u];        // third element

            // butterfly operations
            temp[j] = t + (u1 * cosN(j, s+1)) + (u2 * sinN(j, s+1));
            temp[j + m] = t + (u1 * cosN(j+m, s+1)) + (u2 * sinN(j+m, s+1));

            // update u
            if (u == m) {
                u = pow(2, s+1)-1;
            } else {
                u--;
            }
        }

        // copy temp into input
        for (int i = 0; i < n; i++) {
            input[i] = temp[i];
        }

        // print level data
        printf("\n%d; f%d:\t", rank, s+1);
        for (int i = 0; i < n; i++){
            printf("%.2f\t", input[i]);
        }

        // Gather the data from all processes
        


        MPI_Barrier(MPI_COMM_WORLD);
    }

    // for (int s = 0; s < logN; s++){
    //     int m = pow(2, s);
    //     int u = m;

    //     // Parallelize this loop
    //     for (int k = rank * (n / size); k < (rank + 1) * (n / size); k += 2 * m) {
    //         for (int j = 0; j < m; j++){
    //             if (k + j + m < n) {
    //                 double t  = input[k + j];        // first element
    //                 double u1 = input[k + j + m];    // second element
    //                 double u2 = input[u + k];        // third element

    //                 // butterfly operations
    //                 temp[k + j] = t + (u1 * cosN(j, s+1)) + (u2 * sinN(j, s+1));
    //                 temp[k + j + m] = t + (u1 * cosN(j+m, s+1)) + (u2 * sinN(j+m, s+1));

    //                 // update u
    //                 if (u == m) {
    //                     u = pow(2, s+1)-1;
    //                 } else {
    //                     u--;
    //                 }
    //             }
    //         }
    //     }

    //     // copy temp into input
    //     for (int i = 0; i < n; i++) {
    //         input[i] = temp[i];
    //     }

    //     // print level data
    //     printf("\n%d; f%d:\t", rank, s+1);
    //     for (int i = 0; i < n; i++){
    //         printf("%.2f\t", input[i]);
    //     }

    //     // Gather the data from all processes
    //     // MPI_Allgather(MPI_IN_PLACE, n / size, MPI_DOUBLE, input, n / size, MPI_DOUBLE, MPI_COMM_WORLD);
    // }

    free(temp);
}


/** Function to calculate cos((2 * PI * n) / pow(2, l)) 
    * @param n: index
    * @param l: level
    * @return cos value
*/
double cosN(int n, int l){
    return cos((2 * PI * n) / pow(2, l));
}


/** Function to calculate sin((2 * PI * n) / pow(2, l)) 
    * @param n: index
    * @param l: level
    * @return sin value
*/
double sinN(int n, int l){
    return sin((2 * PI * n) / pow(2, l));
}


/** Function to evaluate H(v)
    * @param input: input array
    * @param n: size of the input array
    * @return void
*/
void eval_Hv(double *input, int n){
    for (int k = 0; k < n; k++){
        input[k] /= n;
    }
}


/** Write the output to the file
    * @param outputfp: file pointer
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

    // Write the output array
    for (int i = 0; i < n; i++) {
        fprintf(outputfp, "%f ", output[i]);
    }

    // Close the file
    fclose(outputfp);

    return 0;
}