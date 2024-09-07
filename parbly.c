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
void eval_levels(double *input, int n, int rank, int p);
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
    int p = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (!rank) {
        if (take_input(inputfp, argv[1], &n, &input)){
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // print the input array
        // printf("\n%d; ar: ", rank);
        // for (int i = 0; i < n; i++) {
        //     printf("%.2f\t", input[i]);
        // }
        // printf("\n");

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
    ebp = eval_ebp(n, p, 1);

    if (!rank) {
        // create level 0: f_0 (tau)
        transform_array(input, n);

        if (ebp == 1){
            // TODO
        }

        // Aggregate the data
        if (aggregate_data(input, p, rank, ebp)){
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // fflush(stdout);
    } else {
        // Receive the data from the root process
        MPI_Recv(input, ebp, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // print level 0
    printf("\n%d; f0:\t", rank);
    for (int i = 0; i < 2; i++) {
        printf("%.2f\t", input[i]);
    }


    // create other levels
    eval_levels(input, n, rank, p);

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
    }

    MPI_Finalize();

    // Write the output to the file
    if (rank == 0) {
        if (print_output(outputfp, argv[2], n, input)){
            return 1;
        }
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

    if (!(p <= p_ideal)){
        return 1;
    }

    return (p_ideal/p) * pow(2, level);
}


/** Function to aggregate the data
    * @param input: input array
    * @param p: number of processes
    * @param ebp: element by process
    * @return 0 if the data is aggregated successfully, 1 otherwise
*/
int aggregate_data(double *input, int p, int rank, int ebp){
    double *couple_input = (double*)malloc(ebp * sizeof(double));

    for (int i = 1; i < p; i++){
        for (int j = 0; j < ebp; j++){
            couple_input[j] = input[(ebp*i)+j];
        }
        // printf("\n%d; MPI_Send: %f,%f to %d", rank, couple_input[0], couple_input[1], i);

        MPI_Send(couple_input, ebp, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }

    return 0;
}


/** Function to evaluate all the f(tau) levels
    * @param input: input array
    * @param n: size of the input array
    * @param rank: rank of the current process
    * @param p: total number of processes
    * @return void
*/
void eval_levels(double *input, int n, int rank, int p){
    int its = pow(2, 1);  // items to send
    double *temp = (double*)malloc(its*2 * sizeof(double));
    int logN = (int)log2(n);

    // initialize temp
    for (int i = 0; i < its*2; i++) {
        temp[i] = NAN;
    }

    // copy the data
    for (int i = 0; i < its; i++){
        temp[i] = input[i];
    }

    // iterate through all levels
    for (int l = 1; l <= logN; l++){
        int m = pow(2, l-1);
        int mod_l = rank % m;
        int u = 0;
        its = pow(2, l);

        // temp = (double*)realloc(temp, (its) * sizeof(double));
        printf("\n%d; f%d -> its: %d", rank, l, its);

        // print temp
        printf("\n%d; temp:\t", rank);
        for (int i = 0; i < its; i++){
            printf("%.2f\t", temp[i]);
        }

        // update u
        if (!mod_l) {
            u = m;
        } else {
            u = (m * 2) - mod_l;
        }

        // evaluate the butterfly operations
        double t  = temp[mod_l];        // first element
        printf("\n%d; f%d -> t=temp[%d]=%f", rank, l, mod_l, t);
        double u1 = temp[mod_l + m];    // second element
        printf("\n%d; f%d -> u1=temp[%d]=%f", rank, l, mod_l+m, u1);
        double u2 = temp[u];        // third element
        printf("\n%d; f%d -> u2=temp[%d]=%f", rank, l, u, u2);

        // initialize memory for temp
        for (int i = 0; i < its; i++) {
            temp[i] = NAN;
        }

        // butterfly operations
        temp[mod_l] = t + (u1 * cosN(mod_l, l)) + (u2 * sinN(mod_l, l));
        printf("\n%d; f%d -> temp[%d]=%f", rank, l, mod_l, temp[mod_l]);
        temp[mod_l + m] = t + (u1 * cosN(mod_l+m, l)) + (u2 * sinN(mod_l+m, l));
        printf("\n%d; f%d -> temp[%d]=%f", rank, l, mod_l+m, temp[mod_l+m]);

        // print level data
        printf("\n%d; f%d:\t", rank, l);
        for (int i = 0; i < its; i++){
            if (temp[i] != NAN){
                printf("%.2f\t", temp[i]);
            } else {
                printf("NAN\t");
            }
        }

        // Gather the data from all processes
        int mod_ts = rank % its;
        double *temp_recv = (double*)malloc((its*2) * sizeof(double));

        if (!mod_ts){
            int p_gap = (pow(2, l) < p) ? pow(2, l) : p;

            for (int i = rank+1; i < rank+p_gap && i < p; i++) {
                printf("\n%d; +++++its %d, mod_ts %d++++++recv1 from: %d", rank, its, mod_ts, i);
                // print temp
                printf("\n%d; pre recv temp:\t", rank);
                for (int j = 0; j < its*2; j++){
                    printf("%.2f\t", temp[j]);
                }
                printf("\n");

                MPI_Recv(temp_recv, its*2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // print temp_recv
                printf("\n%d; post temp_recv:\t", rank);
                for (int j = 0; j < its*2; j++){
                    if (temp_recv[j] != NAN){
                        printf("%.2f\t", temp_recv[j]);
                    } else {
                        printf("NAN\t");
                    }
                }

                // unisci tutto
                if (l == 1){
                    for (int i = its; i < its*2; i++){
                        temp[i] = temp_recv[i-its];
                    }
                } else {
                    for (int i = 0; i < its; i++){
                        if (temp_recv[i] != NAN){
                            temp[i] = temp_recv[i];
                        }
                    }
                }
            }

            printf("\n%d; temp1:\t", rank);
            for (int i = 0; i < its*2; i++){
                printf("%.2f\t", temp[i]);
            }

            if (l != logN) {
                int p_gap = (pow(2, l) < p) ? pow(2, l) : p;

                for (int i = rank+1; i < rank+p_gap && i < p; i++) {
                    printf("\n%d; send2 to: %d -> ", rank, i);
                    for (int i = 0; i < its*2; i++) {
                        printf("%.2f\t", temp[i]);
                    }

                    MPI_Send(temp, its*2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                }
            }
        } else {
            printf("\n%d; send1 to: %d -> ", rank, rank-mod_ts);
            for (int i = 0; i < its*2; i++) {
                printf("%.2f\t", temp[i]);
            }
            MPI_Send(temp, its*2, MPI_DOUBLE, rank-mod_ts, 0, MPI_COMM_WORLD);

            if (l != logN) {
                printf("\n%d; temp3:\t", rank);
                for (int i = 0; i < its*2; i++){
                    printf("%.2f\t", temp[i]);
                }

                printf("\n%d; +++++its %d, mod_ts %d++++++recv2 from: %d", rank, its, mod_ts, rank-mod_ts);
                MPI_Recv(temp, its*2, MPI_DOUBLE, rank-mod_ts, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                printf("\n%d; temp2:\t", rank);
                for (int i = 0; i < its*2; i++){
                    printf("%.2f\t", temp[i]);
                }
            }
        }

        free(temp_recv);
    }

    // copy temp into input
    for (int i = 0; i < n; i++) {
        input[i] = temp[i];
    }

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