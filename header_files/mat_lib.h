// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar

#ifndef MATH_LIBRARY_H
#define MATH_LIBRARY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  

////////////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////////////

void create_matrix(double ***A, int n_rows, int n_cols);
void create_matrix_int(int ***A, int n_rows, int n_cols);
double** create_matrix1(int n_rows, int n_cols);
double* create_vector(int n_rows);
void free_matrix(double **A, int n_rows);
void free_matrix_int(int **A, int n_rows);
void free_vector(double *A);
void print_matrix(double **A, int n_rows, int n_cols);
void multiply_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A, int n_cols_B);
void multiply_matrix_vector(double **A, double *B, double *C, int n_rows_A, int n_cols_A);
void multiply_vector_matrix(double *B, double **A, double **C, int n_rows_A, int n_cols_A);
void multiply_scalar_matrix(double scalar, double **A, double **B, int n_rows_A, int n_cols_A);
void multiply_scalar_vector(double scalar, double *A, double *B, int n_rows_A);
void multiply_vector_matrix_columnwise(double *B, double **A, double *C, int n_rows_A, int n_cols_A);
void add_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void add_matrices_to_first(double **A, double **B, int n_rows_A, int n_cols_A);
void add_vectors(double *A, double *B, double *C, int n_rows_A);
void subtract_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void subtract_vectors(double *A, double *B, double *C, int n_rows_A);
void transpose_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_vector(double *A, double *B, int n_rows_A);
void swap_vectors(double *A, double *B, int n_rows_A);
void vector_outer_product(double *A, double *B, double **C, int n_rows_A, int n_rows_B);
double vector_inner_product(double *A, double *B, int n_rows_A);
double vector_norm(double *A, int n_rows_A);
double** identity_matrix(int n_rows);
void matrixInverse_Gauss_Jordan(double** matrix1, double** matrix, int order);
void luDecomposition(double **A, int n, int *P);
void forwardSubstitution(double **L, double *y, double *b, int n, int *P);
void backwardSubstitution(double **U, double *x, double *y, int n);
void matrixInverse_LU(double **A, double **A_inv, int n);
void printMatrix(double **A, int n);
void matrixInverse_Gauss_Jordan2(double** matrix1, double** inverse, int order);
void write_matrix_to_file(double **A, int n_rows, int n_cols, char *filename);
double l2_norm(double *A, double *B, int n_rows_A);
void linear_system_solver(double** A, double* x, double* b, int n);
void linear_system_solver_jacobi(double** A, double* x, double* b, int n);

void multiply_sparse_matrix_vector(double** D_coeff, double* f, double* dfdx, int** cloud, bool* tag, int n_rows_D, int n_cols_D);
void multiply_sparse_vector_matrix(double* f, double** D_coeff, double** ftimesD, int n_rows_D, int n_cols_D);

////////////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////////////
void multiply_sparse_matrix_vector(double** D_coeff, double* f, double* dfdx, int** cloud, bool* tag, int n_rows_D, int n_cols_D){
    // dfdx = create_matrix1(n_rows_D, n_cols_D);
    double results = 0.0;
    for (int i = 0; i < n_rows_D; i++){
        if (tag[i]==false){
            for (int j = 0; j < n_cols_D; j++){
                results += D_coeff[i][j] * f[cloud[i][j]];
            }
        }
        dfdx[i] = results;
        results = 0;
    }
}

void multiply_sparse_vector_matrix(double* f, double** D_coeff, double** ftimesD, int n_rows_D, int n_cols_D){
    for (int i = 0; i < n_rows_D; i++){
            for (int j = 0; j < n_cols_D; j++){
                ftimesD[i][j] = f[i] * D_coeff[i][j];
            }
        }
    }



double** create_matrix1(int n_rows, int n_cols)
{
    int i;
    double **A;
    A = (double **)malloc(n_rows * sizeof(double *));
    for (i = 0; i < n_rows; i++)
    {
        A[i] = (double *)malloc(n_cols * sizeof(double));
    }
    for (i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            A[i][j] = 0;
        }
    }
    return A;
}

void create_matrix(double ***A, int n_rows, int n_cols)
{
    int i;
    *A = (double **)malloc(n_rows * sizeof(double *));
    for (i = 0; i < n_rows; i++)
    {
        (*A)[i] = (double *)malloc(n_cols * sizeof(double));
    }
    for (i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            (*A)[i][j] = 0;
        }
    }
}

void create_matrix_int(int ***A, int n_rows, int n_cols)
{
    int i;
    *A = (int **)malloc(n_rows * sizeof(int *));
    for (i = 0; i < n_rows; i++)
    {
        (*A)[i] = (int *)malloc(n_cols * sizeof(int));
    }
}

double* create_vector(int n_rows)
{
    double *A;
    A = (double *)malloc(n_rows * sizeof(double));
    for (int i = 0; i < n_rows; i++)
    {
        A[i] = 0;
    }
    return A;
}

void free_matrix(double **A, int n_rows)
{
    int i;
    for (i = 0; i < n_rows; i++)
    {
        free(A[i]);
    }
    free(A);
}

void free_matrix_int(int **A, int n_rows)
{
    int i;
    for (i = 0; i < n_rows; i++)
    {
        free(A[i]);
    }
    free(A);
}

void free_vector(double *A)
{
    free(A);
}

void print_matrix(double **A, int n_rows, int n_cols)
{
    int i, j;
    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_cols; j++)
        {
            printf("%g ", A[i][j]);
        }
        printf("\n");
    }
}

void write_matrix_to_file(double **A, int n_rows, int n_cols, char *filename)
{
    FILE *f = fopen(filename,   "w");   
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    int i, j;
    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_cols; j++)
        {
            fprintf(f, "%g ", A[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void multiply_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A, int n_cols_B)
{
    int i, j, k;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_B; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < n_cols_A; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void multiply_matrix_vector(double **A, double *B, double *C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = 0;
        for (j = 0; j < n_cols_A; j++)
        {
            C[i] += A[i][j] * B[j];
        }
    }
}

void multiply_vector_matrix(double *B, double **A, double **C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            C[i][j] = B[i] * A[i][j];
        }
    }
}

void multiply_scalar_matrix(double scalar, double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[i][j] = scalar * A[i][j];
        }
    }
}

void multiply_scalar_vector(double scalar, double *A, double *B, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        B[i] = scalar * A[i];
    }
}

void add_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void add_matrices_to_first(double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            A[i][j] += B[i][j];
        }
    }
}

void add_vectors(double *A, double *B, double *C, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = A[i] + B[i];
    }
}

void subtract_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void subtract_vectors(double *A, double *B, double *C, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = A[i] - B[i];
    }
}

void transpose_matrix(double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[j][i] = A[i][j];
        }
    }
}

void copy_matrix(double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[i][j] = A[i][j];
        }
    }
}

void copy_vector(double *A, double *B, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        B[i] = A[i];
    }
}

void vector_outer_product(double *A, double *B, double **C, int n_rows_A, int n_rows_B)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_rows_B; j++)
        {
            C[i][j] = A[i] * B[j];
        }
    }
}

double vector_inner_product(double *A, double *B, int n_rows_A)
{
    int i;
    double result = 0;
    for (i = 0; i < n_rows_A; i++)
    {
        result += A[i] * B[i];
    }
    return result;
}

double vector_norm(double *A, int n_rows_A)
{
    int i;
    double result = 0;
    for (i = 0; i < n_rows_A; i++)
    {
        result += A[i] * A[i];
    }
    return sqrt(result);
}

double** identity_matrix(int n_rows)
{
    int i, j;
    double **A;
    A = (double **)malloc(n_rows * sizeof(double *));
    for (i = 0; i < n_rows; i++)
    {
        A[i] = (double *)malloc(n_rows * sizeof(double));
    }
    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_rows; j++)
        {
            if (i == j)
            {
                A[i][j] = 1;
            }
            else
            {
                A[i][j] = 0;
            }
        }
    }
    return A;
}

void matrixInverse_Gauss_Jordan(double** matrix1, double** inverse, int order)
{
    double temp;
    double** matrix;
    create_matrix(&matrix, order, 2 * order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            matrix[i][j] = matrix1[i][j];
        }
    }
    // Create the augmented matrix
    for (int i = 0; i < order; i++) {
        for (int j = order; j < 2 * order; j++) {
            // Add '1' at the diagonal places of
            // the matrix to create a identity matrix
            if (j == (i + order))
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0;
        }
    }

    // Interchange the row of matrix,
    for (int i = order - 1; i > 0; i--) {
        // Directly swapping the rows using pointers saves
        // time
 
        if (matrix[i - 1][0] < matrix[i][0]) {
            double* temp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = temp;
        }
    }
 
    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 0; i < order; i++) {
 
        for (int j = 0; j < order; j++) {
 
            if (j != i) {
 
                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
 
                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {
 
        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {
 
            matrix[i][j] = matrix[i][j] / temp;
        }
    }

    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            inverse[i][j] = matrix[i][j + order];
        }
    }
    free_matrix(matrix, order);
}

// Function to perform LU decomposition
void luDecomposition(double **A, int n, int *P) {
    for (int i = 0; i < n; i++) {
        P[i] = i; // Initialize permutation vector
    }

    for (int k = 0; k < n - 1; k++) {
        int pivot_row = k;
        double pivot = A[k][k];

        // Find the row with the largest pivot
        for (int i = k + 1; i < n; i++) {
            if (abs(A[i][k]) > abs(pivot)) {
                pivot = A[i][k];
                pivot_row = i;
            }
        }

        // Swap rows in permutation vector
        int temp = P[k];
        P[k] = P[pivot_row];
        P[pivot_row] = temp;

        // Swap rows in A
        double *temp_row = A[k];
        A[k] = A[pivot_row];
        A[pivot_row] = temp_row;

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            A[i][k] = factor;
        }
    }
}

// Function to solve Ly = b
void forwardSubstitution(double **L, double *y, double *b, int n, int *P) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[P[i]] - sum) / L[i][i];
    }
}

// Function to solve Ux = y
void backwardSubstitution(double **U, double *x, double *y, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
}

// Function to perform matrix inversion using LU decomposition
void matrixInverse_LU(double **A, double **A_inv, int n) {
    int *P = malloc(n * sizeof(int));

    luDecomposition(A, n, P);

    double **L = malloc(n * sizeof(double *));
    double **U = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        L[i] = malloc(n * sizeof(double));
        U[i] = malloc(n * sizeof(double));
    }

    // Extract L and U matrices
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i > j)
                L[i][j] = A[i][j];
            else if (i == j)
                L[i][j] = 1.0;
            else
                L[i][j] = 0.0;

            if (i <= j)
                U[i][j] = A[i][j];
            else
                U[i][j] = 0.0;
        }
    }

    // Solve Ly = b for each column of inverse
    double *y = malloc(n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                b[j] = 1.0;
            else
                b[j] = 0.0;
        }
        forwardSubstitution(L, y, b, n, P);
        backwardSubstitution(U, A_inv[i], y, n);
    }

    free(P);
    free(y);
    free(b);
    for (int i = 0; i < n; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
}

// Function to print a matrix
void printMatrix(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.2f ", A[i][j]);
        }
        printf("\n");
    }
}


void matrixInverse_Gauss_Jordan2(double** matrix1, double** inverse, int order)
{
    double temp;
    double** matrix;
    create_matrix(&matrix, order, 2 * order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            matrix[i][j] = matrix1[i][j];
        }
    }

    // Create the augmented matrix
    for (int i = 0; i < order; i++) {
        for (int j = order; j < 2 * order; j++) {
            // Add '1' at the diagonal places of
            // the matrix to create a identity matrix
            if (j == (i + order))
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0;
        }
    }

    // Interchange the row of matrix,
    for (int i = order - 1; i > 0; i--) {
        // Directly swapping the rows using pointers saves
        // time
 
        if (matrix[i - 1][0] < matrix[i][0]) {
            double* temp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = temp;
        }
    }
 
    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 0; i < order; i++) {
 
        for (int j = 0; j < order; j++) {
 
            if (j != i) {
 
                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
 
                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {
 
        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {
 
            matrix[i][j] = matrix[i][j] / temp;
        }
    }

    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            inverse[i][j] = matrix[i][j + order];
        }
    }
    free_matrix(matrix, order);
}

double l2_norm(double *A, double *B, int n_rows_A)
{
    int i;
    double result = 0;
    for (i = 0; i < n_rows_A; i++)
    {
        result += (A[i] - B[i]) * (A[i] - B[i]);
    }
    return sqrt(result/n_rows_A);
}

void linear_system_solver(double** A, double* x, double* b, int n)
{
    double** A_inv = create_matrix1(n, n);
    matrixInverse_Gauss_Jordan(A, A_inv, n);
    multiply_matrix_vector(A_inv, b, x, n, n);
    free_matrix(A_inv, n);
}

void swap_vectors(double *A, double *B, int n_rows_A)
{
    double temp;
    for (int i = 0; i < n_rows_A; i++)
    {
        temp = A[i];
        A[i] = B[i];
        B[i] = temp;
    }
}

void multiply_vector_matrix_columnwise(double *B, double **A, double *C, int n_rows_A, int n_cols_A){
    int i, j;
    for (i = 0; i < n_cols_A; i++)
    {
        C[i] = 0;
        for (j = 0; j < n_rows_A; j++)
        {
            C[i] += B[j] * A[j][i];
        }
    }
}


#endif