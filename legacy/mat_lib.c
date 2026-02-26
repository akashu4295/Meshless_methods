// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include "functions.h"
#include "structures.h"

////////////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////////////
double min3(double a, double b, double c){
    double min = a;
    if (b < min)
        min = b;
    if (c < min)
        min = c;
    return min;
}

void multiply_sparse_matrix_vector_vectorised(const double *D_coeff, const double *f, double *dfdx, const int *cloud, int n_rows_D, int n_cols_D)
{
    for (int i = 0; i < n_rows_D; i++) {
        double result = 0.0;
        int row_start = i * n_cols_D;
        for (int j = 0; j < n_cols_D; j++) {
            int idx = cloud[row_start + j];
            result += D_coeff[row_start + j] * f[idx];
        }
        dfdx[i] = result;
    }
}

void multiply_sparse_matrix_vector_vectorised_gpu(double *D_coeff, double *f, double *dfdx,
                                                  int *cloud, int n_rows_D, int n_cols_D)
{
    // Expectation: the host caller must ensure D_coeff, f, dfdx, cloud are present on the device
    // (or the function is called inside an acc data region where they are present).

    // Run the outer loop on the device in parallel
    #pragma acc parallel loop gang present(D_coeff, f, dfdx, cloud)
    for (int i = 0; i < n_rows_D; ++i) {
        double result = 0.0;
        int row_start = i * n_cols_D;

        // inner loop is sequential for each i (you can change to "vector" if beneficial)
        #pragma acc loop seq
        for (int j = 0; j < n_cols_D; ++j) {
            int idx = cloud[row_start + j];
            // sanity: idx should be within bounds of f[] on device
            result += D_coeff[row_start + j] * f[idx];
        }
        dfdx[i] = result;
    }
}

void multiply_sparse_matrix_vector_vectorised_gpu_async(double *D_coeff, double *f, double *dfdx,
                                                  int *cloud, int n_rows_D, int n_cols_D, int async_queue)
{
    // Run the outer loop on the device in parallel
    #pragma acc parallel loop gang async(async_queue) present(D_coeff, f, dfdx, cloud)
    for (int i = 0; i < n_rows_D; ++i) {
        double result = 0.0;
        int row_start = i * n_cols_D;

        // inner loop is sequential for each i (you can change to "vector" if beneficial)
        #pragma acc loop seq
        for (int j = 0; j < n_cols_D; ++j) {
            int idx = cloud[row_start + j];
            // sanity: idx should be within bounds of f[] on device
            result += D_coeff[row_start + j] * f[idx];
        }
        dfdx[i] = result;
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

double* create_matrix_vectorised(int n_rows, int n_cols)
{
    int i;
    double *A;
    A = (double *)malloc(n_rows * n_cols * sizeof(double));
    for (i = 0; i < n_rows*n_cols; i++)
            A[i] = 0;
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

double* create_vector(int n_rows)
{
    double *A;
    A = (double *)malloc(n_rows * sizeof(double));
    for (int i = 0; i < n_rows; i++)
    {
        A[i] = 0.0;
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
void multiply_matrices_vectorised(const double *A, const double *B, double *C, int n_rows_A, int n_cols_A, int n_cols_B)
{
    int t1, t2; double sum;
    for (int i = 0; i < n_rows_A; i++) {
        t1 = i * n_cols_B;
        t2 = i * n_cols_A;
        for (int j = 0; j < n_cols_B; j++) {
            sum = 0.0;
            for (int k = 0; k < n_cols_A; k++) {
                sum += A[t2 + k] * B[k * n_cols_B + j];
            }
            C[t1 + j] = sum;
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

void matrixInverse_Gauss_Jordan(double** matrix1, double** inverse, int order)
{
    double temp1;
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
 
                temp1 = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
 
                    matrix[j][k] -= matrix[i][k] * temp1;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {
 
        temp1 = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {
 
            matrix[i][j] = matrix[i][j] / temp1;
        }
    }

    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            inverse[i][j] = matrix[i][j + order];
        }
    }
    free_matrix(matrix, order);
}

void matrixInverse_Gauss_Jordan_vectorised(double* A, double* Ainv, int n){
    // Create augmented matrix [A | I]
    double *aug = (double *)malloc(sizeof(double) * n * 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * 2 * n + j] = A[i * n + j];
            aug[i * 2 * n + (j + n)] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss–Jordan elimination
    for (int i = 0; i < n; i++) {
        // Pivot selection (partial)
        double maxEl = fabs(aug[i * 2 * n + i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            double val = fabs(aug[k * 2 * n + i]);
            if (val > maxEl) {
                maxEl = val;
                maxRow = k;
            }
        }

        // Swap rows if needed
        if (maxRow != i) {
            for (int k = 0; k < 2 * n; k++) {
                double tmp = aug[i * 2 * n + k];
                aug[i * 2 * n + k] = aug[maxRow * 2 * n + k];
                aug[maxRow * 2 * n + k] = tmp;
            }
        }

        // Normalize pivot row
        double pivot = aug[i * 2 * n + i];
        if (fabs(pivot) < 1e-14) {
            fprintf(stderr, "Matrix is singular or nearly singular.\n");
            free(aug);
            return;
        }
        for (int k = 0; k < 2 * n; k++)
            aug[i * 2 * n + k] /= pivot;

        // Eliminate all other rows
        for (int r = 0; r < n; r++) {
            if (r != i) {
                double factor = aug[r * 2 * n + i];
                for (int c = 0; c < 2 * n; c++) {
                    aug[r * 2 * n + c] -= factor * aug[i * 2 * n + c];
                }
            }
        }
    }

    // Extract inverse
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Ainv[i * n + j] = aug[i * 2 * n + (j + n)];

    free(aug);
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

void multiply_vector_matrix_columnwise_vectorised(double *B, double *A, double *C, int n_rows_A, int n_cols_A) {
    for (int i = 0; i < n_cols_A; i++) {
        C[i] = 0.0;
        for (int j = 0; j < n_rows_A; j++) {
            C[i] += B[j] * A[j * n_cols_A + i];
        }
    }
}