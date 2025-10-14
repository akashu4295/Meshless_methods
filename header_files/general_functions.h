// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// Functions used to write the output to files

#ifndef GENERAL_FUNCTIONS_H
#define GENERAL_FUNCTIONS_H

#include "structures.h"
#include "mat_lib.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

/////////////////////////////////////////////////////////////////////////////////
// Function Declarations
/////////////////////////////////////////////////////////////////////////////////

void calculate_errors_2d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* lapf, double* F, double* Fx, double* Fy, double* lapF);
void calculate_errors_3d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* fz, double* lapf, double* F, double* Fx, double* Fy, double* Fz, double* lapF);
void calculate_test_vectors_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf);
void calculate_test_vectors_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf);
void set_manufactured_solution_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf, int k);
void set_manufactured_solution_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf, int k);
void test_derivatives(PointStructure* myPointStruct, short num_levels, short domain_dimension);
void free_test_vectors(double** f, double** fx, double** fy, double** fz, double** lapf, double*** Dx, double*** Dy, double*** Dz, double*** lap, int num_nodes);
double l2_norm_gen(PointStructure* myPointStruct, double *A, double *B, int n_rows_A);
void make_directory(const char* name);
double calculate_dt(PointStructure* myPointStruct);
/////////////////////////////////////////////////////////////////////////////////
// Function Definitions
/////////////////////////////////////////////////////////////////////////////////
double calculate_dt(PointStructure* myPointStruct){
    double dt = 1e10;
    double d = myPointStruct->d_avg;
    double termd =    2*parameters.nu/(d*d);
    double termc =    1/d;
    dt = 0.5/(termd + termc);
    return (dt*parameters.courant_number);
}


void test_derivatives(PointStructure* myPointStruct, short num_levels, short domain_dimension){
    for (int ii = 0; ii < num_levels; ii++) {
        // Manufactured RHS function
        printf("Testing level %d\n", ii);
        int k = 1; // wave number
        if (domain_dimension > 2) {
            double *f, *fx, *fy, *fz, *lapf;
            f = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            fx = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            fy = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            fz = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            lapf = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            
            // Ensure the data is copied to the device (GPU)
            
            set_manufactured_solution_3d(&myPointStruct[ii], &f, &fx, &fy, &fz, &lapf, k);
            calculate_test_vectors_3d(&myPointStruct[ii], &f, &fx, &fy, &fz, &lapf);
            
            double *F, *Fx, *Fy, *Fz, *lapF;
            F = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            Fx = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            Fy = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            Fz = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            lapF = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            
            // Allocate memory for the new data on the device
            
            set_manufactured_solution_3d(&myPointStruct[ii], &F, &Fx, &Fy, &Fz, &lapF, k);
            calculate_errors_3d(&myPointStruct[ii], f, fx, fy, fz, lapf, F, Fx, Fy, Fz, lapF);
            
            // Free data on the device
            
            free(f); free(fx); free(fy); free(fz); free(lapf); free(F); free(Fx); free(Fy); free(Fz); free(lapF);
        } else {
            double *f, *fx, *fy, *lapf;
            f = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            fx = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            fy = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            lapf = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            
            
            set_manufactured_solution_2d(&myPointStruct[ii], &f, &fx, &fy, &lapf, k);
            calculate_test_vectors_2d(&myPointStruct[ii], &f, &fx, &fy, &lapf);
            
            double *F, *Fx, *Fy, *lapF;
            F = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            Fx = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            Fy = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            lapF = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
            
            
            set_manufactured_solution_2d(&myPointStruct[ii], &F, &Fx, &Fy, &lapF, k);
            calculate_errors_2d(&myPointStruct[ii], f, fx, fy, lapf, F, Fx, Fy, lapF);
            
            free(f); free(fx); free(fy); free(lapf); free(F); free(Fx); free(Fy); free(lapF);
        }
    }
}

void calculate_test_vectors_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf){
    multiply_sparse_matrix_vector(myPointStruct->Dx, *f, *fx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->Dy, *f, *fy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->Dz, *f, *fz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->lap, *f, *lapf, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
}

void calculate_test_vectors_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf){    
    multiply_sparse_matrix_vector(myPointStruct->Dx, *f, *fx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->Dy, *f, *fy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->lap, *f, *lapf, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
}

void set_manufactured_solution_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf, int k) {
    int num_nodes = myPointStruct->num_nodes;

    // Offload the loop to the GPU
    for (int i = 0; i < num_nodes; i++) {
        (*f)[i] = sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                  sin(2 * 3.14 * k * myPointStruct->y[i]) * 
                  sin(2 * 3.14 * k * myPointStruct->z[i]);

        (*fx)[i] = 2 * 3.14 * k * cos(2 * 3.14 * k * myPointStruct->x[i]) * 
                   sin(2 * 3.14 * k * myPointStruct->y[i]) * 
                   sin(2 * 3.14 * k * myPointStruct->z[i]);

        (*fy)[i] = 2 * 3.14 * k * sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                   cos(2 * 3.14 * k * myPointStruct->y[i]) * 
                   sin(2 * 3.14 * k * myPointStruct->z[i]);

        (*fz)[i] = 2 * 3.14 * k * sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                   sin(2 * 3.14 * k * myPointStruct->y[i]) * 
                   cos(2 * 3.14 * k * myPointStruct->z[i]);

        (*lapf)[i] = -12 * 3.14 * 3.14 * k * k * 
                     sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                     sin(2 * 3.14 * k * myPointStruct->y[i]) * 
                     sin(2 * 3.14 * k * myPointStruct->z[i]);
    }
}

void set_manufactured_solution_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf, int k) {
    int num_nodes = myPointStruct->num_nodes;
    
    // Offload the loop to the GPU
    for (int i = 0; i < num_nodes; i++) {
        (*f)[i] = sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                  sin(2 * 3.14 * k * myPointStruct->y[i]);

        (*fx)[i] = 2 * 3.14 * k * cos(2 * 3.14 * k * myPointStruct->x[i]) * 
                   sin(2 * 3.14 * k * myPointStruct->y[i]);

        (*fy)[i] = 2 * 3.14 * k * sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                   cos(2 * 3.14 * k * myPointStruct->y[i]);

        (*lapf)[i] = -8 * 3.14 * 3.14 * k * k * 
                     sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                     sin(2 * 3.14 * k * myPointStruct->y[i]);
    }
}

void calculate_errors_2d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* lapf, double* F, double* Fx, double* Fy, double* lapF){
    double l2_error_f = l2_norm_gen(myPointStruct, F, f, myPointStruct->num_nodes);
    double l2_error_fx = l2_norm_gen(myPointStruct, Fx, fx, myPointStruct->num_nodes);
    double l2_error_fy = l2_norm_gen(myPointStruct, Fy, fy, myPointStruct->num_nodes);
    double l2_error_lapf = l2_norm_gen(myPointStruct, lapF, lapf, myPointStruct->num_nodes);
    printf("L2 error in f, fx, fy, lapf:\n%lf %lf %lf %lf\n", l2_error_f, l2_error_fx, l2_error_fy, l2_error_lapf);
}

void calculate_errors_3d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* fz, double* lapf, double* F, double* Fx, double* Fy, double* Fz, double* lapF){
    double l2_error_f = l2_norm_gen(myPointStruct, F, f, myPointStruct->num_nodes);
    double l2_error_fx = l2_norm_gen(myPointStruct, Fx, fx, myPointStruct->num_nodes);
    double l2_error_fy = l2_norm_gen(myPointStruct, Fy, fy, myPointStruct->num_nodes);
    double l2_error_fz = l2_norm_gen(myPointStruct, Fz, fz, myPointStruct->num_nodes);
    double l2_error_lapf = l2_norm_gen(myPointStruct, lapF, lapf, myPointStruct->num_nodes);
    printf("L2 error in f, fx, fy, fz, lapf:\n%lf %lf %lf %lf %lf\n", l2_error_f, l2_error_fx, l2_error_fy, l2_error_fz, l2_error_lapf);
}

double l2_norm_gen(PointStructure* myPointStruct, double *A, double *B, int n_rows_A) {
    double result = 0.0;
    // Offload the loop to the GPU
    for (int i = 0; i < n_rows_A; i++) {
        if (myPointStruct->boundary_tag[i] == 0) {
            double diff = A[i] - B[i];
            result += diff * diff;
        }
    }
    // Normalize the result after summing up all the differences
    return sqrt(result / n_rows_A);
}

void free_test_vectors(double** f, double** fx, double** fy, double** fz, double** lapf, double*** Dx, double*** Dy, double*** Dz, double*** lap, int num_nodes)
{
    free(*f);
    free(*fx);
    free(*fy);
    free(*fz);
    free(*lapf);
    free_matrix(*Dx, num_nodes);
    free_matrix(*Dy, num_nodes);
    free_matrix(*Dz, num_nodes);
    free_matrix(*lap, num_nodes);
}

void set_source_term(PointStructure* myPointStruct, double* source, int k) {
    if (parameters.dimension == 2) {
        // Parallelize the 2D case
        for (int i = 0; i < myPointStruct->num_nodes; i++) {
            source[i] = -8 * 3.14 * 3.14 * k * k * 
                        sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                        sin(2 * 3.14 * k * myPointStruct->y[i]);
        }
    } else {
        // Parallelize the 3D case
        for (int i = 0; i < myPointStruct->num_nodes; i++) {
            source[i] = -12 * 3.14 * 3.14 * k * k * 
                        sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                        sin(2 * 3.14 * k * myPointStruct->y[i]) * 
                        sin(2 * 3.14 * k * myPointStruct->z[i]);
        }
    }
}

void make_directory(const char* name) {
   #ifdef __linux__
       mkdir(name, 777); 
   #else
       mkdir(name);
   #endif
}

void read_solution_csv(FieldVariables* field, PointStructure* myPointStruct, const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        return;
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &myPointStruct->x[i], &myPointStruct->y[i], &myPointStruct->z[i], &field->u[i], &field->v[i], &field->w[i], &field->p[i]);
    }
    fclose(file);
}

#endif
