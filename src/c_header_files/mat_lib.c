// Author : Akash Unnikrishnan
// Clean SIMD + GPU safe version

#include "functions.h"

#define ALIGNMENT 64
#define PIVOT_TOL 1e-14   /* TODO: expose as user parameter */

////////////////////////////////////////////////////////////////////////
// Safe aligned allocation (CPU SIMD friendly, GPU safe)
////////////////////////////////////////////////////////////////////////

static void* safe_malloc(size_t size)
{
    void *ptr = NULL;

    if (posix_memalign(&ptr, ALIGNMENT, size) != 0) {
        fprintf(stderr,"Allocation failed\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

////////////////////////////////////////////////////////////////////////
// Basic utilities
////////////////////////////////////////////////////////////////////////

double min3(double a,double b,double c)
{
    double m=a;
    if(b<m) m=b;
    if(c<m) m=c;
    return m;
}

////////////////////////////////////////////////////////////////////////
// Sparse matrix-vector
////////////////////////////////////////////////////////////////////////

void multiply_sparse_matrix_vector_vectorised(
    const double * restrict D_coeff,
    const double * restrict f,
    double * restrict dfdx,
    const int * restrict cloud,
    int n_rows_D,
    int n_cols_D)
{
    for(int i=0;i<n_rows_D;i++){

        double sum=0.0;
        int base=i*n_cols_D;

#pragma omp simd reduction(+:sum)
        for(int j=0;j<n_cols_D;j++)
            sum+=D_coeff[base+j]*f[cloud[base+j]];

        dfdx[i]=sum;
    }
}

////////////////////////////////////////////////////////////////////////
// GPU versions (unchanged behaviour)
////////////////////////////////////////////////////////////////////////

void multiply_sparse_matrix_vector_vectorised_gpu(
    double *D_coeff,double *f,double *dfdx,
    int *cloud,int n_rows_D,int n_cols_D)
{
#pragma acc parallel loop gang present(D_coeff,f,dfdx,cloud)
    for(int i=0;i<n_rows_D;i++){

        double sum=0.0;
        int base=i*n_cols_D;

#pragma acc loop vector
        for(int j=0;j<n_cols_D;j++)
            sum+=D_coeff[base+j]*f[cloud[base+j]];

        dfdx[i]=sum;
    }
}

void multiply_sparse_matrix_vector_vectorised_gpu_async(
    double *D_coeff,double *f,double *dfdx,
    int *cloud,int n_rows_D,int n_cols_D,
    int async_queue)
{
#pragma acc parallel loop gang async(async_queue) present(D_coeff,f,dfdx,cloud)
    for(int i=0;i<n_rows_D;i++){

        double sum=0.0;
        int base=i*n_cols_D;

#pragma acc loop vector
        for(int j=0;j<n_cols_D;j++)
            sum+=D_coeff[base+j]*f[cloud[base+j]];

        dfdx[i]=sum;
    }
}

////////////////////////////////////////////////////////////////////////
// Matrix creation
////////////////////////////////////////////////////////////////////////

double** create_matrix1(int n_rows,int n_cols)
{
    double **A=safe_malloc(n_rows*sizeof(*A));

    for(int i=0;i<n_rows;i++){
        A[i]=safe_malloc(n_cols*sizeof(**A));
        for(int j=0;j<n_cols;j++)
            A[i][j]=0.0;
    }
    return A;
}

double* create_matrix_vectorised(int n_rows,int n_cols)
{
    size_t N=(size_t)n_rows*n_cols;
    double *A=safe_malloc(N*sizeof(*A));

    for(size_t i=0;i<N;i++) A[i]=0.0;
    return A;
}

void create_matrix(double ***A,int n_rows,int n_cols)
{
    *A=create_matrix1(n_rows,n_cols);
}

double* create_vector(int n)
{
    double *v=safe_malloc(n*sizeof(*v));
    for(int i=0;i<n;i++) v[i]=0.0;
    return v;
}

void free_matrix(double **A,int n_rows)
{
    for(int i=0;i<n_rows;i++) free(A[i]);
    free(A);
}

////////////////////////////////////////////////////////////////////////
// Dense matrix multiply
////////////////////////////////////////////////////////////////////////

void multiply_matrices(
    double **A,double **B,double **C,
    int n_rows_A,int n_cols_A,int n_cols_B)
{
    for(int i=0;i<n_rows_A;i++)
        for(int j=0;j<n_cols_B;j++){

            double sum=0.0;

#pragma omp simd reduction(+:sum)
            for(int k=0;k<n_cols_A;k++)
                sum+=A[i][k]*B[k][j];

            C[i][j]=sum;
        }
}

void multiply_matrices_vectorised(
    const double * restrict A,
    const double * restrict B,
    double * restrict C,
    int n_rows_A,int n_cols_A,int n_cols_B)
{
    for(int i=0;i<n_rows_A;i++){

        int rowA=i*n_cols_A;
        int rowC=i*n_cols_B;

        for(int j=0;j<n_cols_B;j++){

            double sum=0.0;

#pragma omp simd reduction(+:sum)
            for(int k=0;k<n_cols_A;k++)
                sum+=A[rowA+k]*B[k*n_cols_B+j];

            C[rowC+j]=sum;
        }
    }
}

////////////////////////////////////////////////////////////////////////
// Matrix-vector
////////////////////////////////////////////////////////////////////////

void multiply_matrix_vector(
    double **A,double *B,double *C,
    int n_rows_A,int n_cols_A)
{
    for(int i=0;i<n_rows_A;i++){

        double sum=0.0;

#pragma omp simd reduction(+:sum)
        for(int j=0;j<n_cols_A;j++)
            sum+=A[i][j]*B[j];

        C[i]=sum;
    }
}

/* TODO: rename -> scale_rows */
void multiply_vector_matrix(
    double *B,double **A,double **C,
    int n_rows_A,int n_cols_A)
{
    for(int i=0;i<n_rows_A;i++)
        for(int j=0;j<n_cols_A;j++)
            C[i][j]=B[i]*A[i][j];
}

////////////////////////////////////////////////////////////////////////
// Norms
////////////////////////////////////////////////////////////////////////

double vector_norm(double *A,int n)
{
    double sum=0.0;

#pragma omp simd reduction(+:sum)
    for(int i=0;i<n;i++)
        sum+=A[i]*A[i];

    return sqrt(sum);
}

/* TODO: actually RMS norm */
double l2_norm(double *A,double *B,int n)
{
    double sum=0.0;

#pragma omp simd reduction(+:sum)
    for(int i=0;i<n;i++){
        double d=A[i]-B[i];
        sum+=d*d;
    }

    return sqrt(sum/n);
}

////////////////////////////////////////////////////////////////////////
// Robust Gauss–Jordan inverse
////////////////////////////////////////////////////////////////////////

void matrixInverse_Gauss_Jordan_vectorised(
    double *A,double *Ainv,int n)
{
    size_t W=2*(size_t)n;
    double *aug=safe_malloc(sizeof(double)*n*W);

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){
            aug[i*W+j]=A[i*n+j];
            aug[i*W+j+n]=(i==j);
        }

    for(int i=0;i<n;i++){

        int pivot=i;
        double max=fabs(aug[i*W+i]);

        for(int r=i+1;r<n;r++){
            double v=fabs(aug[r*W+i]);
            if(v>max){max=v;pivot=r;}
        }

        if(max<PIVOT_TOL){
            fprintf(stderr,"Singular matrix\n");
            free(aug);
            return;
        }

        if(pivot!=i)
            for(size_t c=0;c<W;c++){
                double t=aug[i*W+c];
                aug[i*W+c]=aug[pivot*W+c];
                aug[pivot*W+c]=t;
            }

        double piv=aug[i*W+i];

#pragma omp simd
        for(size_t c=0;c<W;c++)
            aug[i*W+c]/=piv;

        for(int r=0;r<n;r++){
            if(r==i) continue;

            double f=aug[r*W+i];

#pragma omp simd
            for(size_t c=0;c<W;c++)
                aug[r*W+c]-=f*aug[i*W+c];
        }
    }

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            Ainv[i*n+j]=aug[i*W+j+n];

    free(aug);
}

void matrixInverse_Gauss_Jordan(
    double **matrix1,double **inverse,int n)
{
    /* TODO: eventually remove pointer-based interface */

    double *A=create_matrix_vectorised(n,n);
    double *Ai=create_matrix_vectorised(n,n);

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            A[i*n+j]=matrix1[i][j];

    matrixInverse_Gauss_Jordan_vectorised(A,Ai,n);

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            inverse[i][j]=Ai[i*n+j];

    free(A);
    free(Ai);
}

////////////////////////////////////////////////////////////////////////

void multiply_vector_matrix_columnwise_vectorised(
    double *B,double *A,double *C,
    int n_rows_A,int n_cols_A)
{
    for(int i=0;i<n_cols_A;i++){

        double sum=0.0;

#pragma omp simd reduction(+:sum)
        for(int j=0;j<n_rows_A;j++)
            sum+=B[j]*A[j*n_cols_A+i];

        C[i]=sum;
    }
}