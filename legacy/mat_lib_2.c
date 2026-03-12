// Author : Akash Unnikrishnan
// Cleaned version (safety + numerical robustness improvements)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "functions.h"
#include "structures.h"

#define PIVOT_TOL 1e-14   /* TODO: make tolerance a function argument */

////////////////////////////////////////////////////////////////////////
// Utility
////////////////////////////////////////////////////////////////////////

static void* safe_malloc(size_t size)
{
    void *ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

////////////////////////////////////////////////////////////////////////
// Basic helpers
////////////////////////////////////////////////////////////////////////

double min3(double a, double b, double c)
{
    double min = a;
    if (b < min) min = b;
    if (c < min) min = c;
    return min;
}

////////////////////////////////////////////////////////////////////////
// Sparse Matrix × Vector
////////////////////////////////////////////////////////////////////////

void multiply_sparse_matrix_vector_vectorised(
    const double *D_coeff,
    const double *f,
    double *dfdx,
    const int *cloud,
    int n_rows_D,
    int n_cols_D)
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

void multiply_sparse_matrix_vector_vectorised_gpu(
    double *D_coeff,
    double *f,
    double *dfdx,
    int *cloud,
    int n_rows_D,
    int n_cols_D)
{
#pragma acc parallel loop gang present(D_coeff,f,dfdx,cloud)
    for (int i = 0; i < n_rows_D; ++i) {

        double result = 0.0;
        int row_start = i * n_cols_D;

#pragma acc loop vector
        for (int j = 0; j < n_cols_D; ++j) {
            int idx = cloud[row_start + j];
            result += D_coeff[row_start + j] * f[idx];
        }

        dfdx[i] = result;
    }
}

void multiply_sparse_matrix_vector_vectorised_gpu_async(
    double *D_coeff,
    double *f,
    double *dfdx,
    int *cloud,
    int n_rows_D,
    int n_cols_D,
    int async_queue)
{
#pragma acc parallel loop gang async(async_queue) present(D_coeff,f,dfdx,cloud)
    for (int i = 0; i < n_rows_D; ++i) {

        double result = 0.0;
        int row_start = i * n_cols_D;

#pragma acc loop vector
        for (int j = 0; j < n_cols_D; ++j) {
            int idx = cloud[row_start + j];
            result += D_coeff[row_start + j] * f[idx];
        }

        dfdx[i] = result;
    }
}

////////////////////////////////////////////////////////////////////////
// Matrix creation
////////////////////////////////////////////////////////////////////////

double** create_matrix1(int n_rows, int n_cols)
{
    double **A = safe_malloc(n_rows * sizeof(*A));

    for (int i = 0; i < n_rows; i++) {
        A[i] = safe_malloc(n_cols * sizeof(**A));
        for (int j = 0; j < n_cols; j++)
            A[i][j] = 0.0;
    }
    return A;
}

double* create_matrix_vectorised(int n_rows, int n_cols)
{
    size_t N = (size_t)n_rows * n_cols;
    double *A = safe_malloc(N * sizeof(*A));

    for (size_t i = 0; i < N; i++)
        A[i] = 0.0;

    return A;
}

void create_matrix(double ***A, int n_rows, int n_cols)
{
    *A = create_matrix1(n_rows, n_cols);
}

double* create_vector(int n_rows)
{
    double *A = safe_malloc(n_rows * sizeof(*A));

    for (int i = 0; i < n_rows; i++)
        A[i] = 0.0;

    return A;
}

void free_matrix(double **A, int n_rows)
{
    for (int i = 0; i < n_rows; i++)
        free(A[i]);

    free(A);
}

////////////////////////////////////////////////////////////////////////
// Matrix Multiplication
////////////////////////////////////////////////////////////////////////

void multiply_matrices(
    double **A,
    double **B,
    double **C,
    int n_rows_A,
    int n_cols_A,
    int n_cols_B)
{
    for (int i = 0; i < n_rows_A; i++)
        for (int j = 0; j < n_cols_B; j++) {

            double sum = 0.0;

            for (int k = 0; k < n_cols_A; k++)
                sum += A[i][k] * B[k][j];

            C[i][j] = sum;
        }
}

void multiply_matrices_vectorised(
    const double *A,
    const double *B,
    double *C,
    int n_rows_A,
    int n_cols_A,
    int n_cols_B)
{
    for (int i = 0; i < n_rows_A; i++) {

        int rowA = i * n_cols_A;
        int rowC = i * n_cols_B;

        for (int j = 0; j < n_cols_B; j++) {

            double sum = 0.0;

            for (int k = 0; k < n_cols_A; k++)
                sum += A[rowA + k] * B[k * n_cols_B + j];

            C[rowC + j] = sum;
        }
    }
}

////////////////////////////////////////////////////////////////////////
// Matrix–Vector
////////////////////////////////////////////////////////////////////////

void multiply_matrix_vector(
    double **A,
    double *B,
    double *C,
    int n_rows_A,
    int n_cols_A)
{
    for (int i = 0; i < n_rows_A; i++) {

        double sum = 0.0;

        for (int j = 0; j < n_cols_A; j++)
            sum += A[i][j] * B[j];

        C[i] = sum;
    }
}

void multiply_vector_matrix(
    double *B,
    double **A,
    double **C,
    int n_rows_A,
    int n_cols_A)
{
    /* TODO: rename -> scale_rows() */
    for (int i = 0; i < n_rows_A; i++)
        for (int j = 0; j < n_cols_A; j++)
            C[i][j] = B[i] * A[i][j];
}

////////////////////////////////////////////////////////////////////////
// Norms
////////////////////////////////////////////////////////////////////////

double vector_norm(double *A, int n)
{
    double result = 0.0;

    for (int i = 0; i < n; i++)
        result += A[i] * A[i];

    return sqrt(result);
}

double l2_norm(double *A, double *B, int n)
{
    /* TODO: this is RMS error, not true L2 norm */
    double result = 0.0;

    for (int i = 0; i < n; i++) {
        double d = A[i] - B[i];
        result += d * d;
    }

    return sqrt(result / n);
}

////////////////////////////////////////////////////////////////////////
// Gauss–Jordan Inverse (ROBUST VERSION)
////////////////////////////////////////////////////////////////////////

void matrixInverse_Gauss_Jordan_vectorised(
    double *A,
    double *Ainv,
    int n)
{
    size_t width = 2 * (size_t)n;
    double *aug = safe_malloc(sizeof(double) * n * width);

    /* Build augmented matrix */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {

            aug[i*width + j]     = A[i*n + j];
            aug[i*width + j + n] = (i == j);
        }

    for (int i = 0; i < n; i++) {

        /* Partial pivoting */
        int maxRow = i;
        double maxVal = fabs(aug[i*width + i]);

        for (int r = i+1; r < n; r++) {
            double val = fabs(aug[r*width + i]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = r;
            }
        }

        if (maxVal < PIVOT_TOL) {
            fprintf(stderr,"Matrix singular.\n");
            free(aug);
            return;
        }

        if (maxRow != i)
            for (size_t c = 0; c < width; c++) {
                double tmp = aug[i*width+c];
                aug[i*width+c] = aug[maxRow*width+c];
                aug[maxRow*width+c] = tmp;
            }

        double pivot = aug[i*width+i];

        for (size_t c = 0; c < width; c++)
            aug[i*width+c] /= pivot;

        for (int r = 0; r < n; r++) {

            if (r == i) continue;

            double factor = aug[r*width+i];

            for (size_t c = 0; c < width; c++)
                aug[r*width+c] -= factor * aug[i*width+c];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Ainv[i*n+j] = aug[i*width + j + n];

    free(aug);
}

void matrixInverse_Gauss_Jordan(
    double **matrix1,
    double **inverse,
    int order)
{
    /* TODO: consider removing pointer-based version */
    double *A    = create_matrix_vectorised(order, order);
    double *Ainv = create_matrix_vectorised(order, order);

    for (int i=0;i<order;i++)
        for (int j=0;j<order;j++)
            A[i*order+j] = matrix1[i][j];

    matrixInverse_Gauss_Jordan_vectorised(A,Ainv,order);

    for (int i=0;i<order;i++)
        for (int j=0;j<order;j++)
            inverse[i][j] = Ainv[i*order+j];

    free(A);
    free(Ainv);
}

////////////////////////////////////////////////////////////////////////

void multiply_vector_matrix_columnwise_vectorised(
    double *B,
    double *A,
    double *C,
    int n_rows_A,
    int n_cols_A)
{
    for (int i = 0; i < n_cols_A; i++) {

        double sum = 0.0;

        for (int j = 0; j < n_rows_A; j++)
            sum += B[j] * A[j*n_cols_A + i];

        C[i] = sum;
    }
}