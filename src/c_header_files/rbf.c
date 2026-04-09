// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
#include "structures.h"
#include "functions.h"
#include "kdtree.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// Function definitions
double calculate_phs_rbf(double *pt1, double *pt2, int phs, int dimension) {
    double sum = 0;
    for (int i = 0; i < dimension; i++) {
        sum += pow(pt1[i] - pt2[i], 2);
    }
    if (sum == 0)
        return 0;
    return pow(sum, phs*0.5);
}

void create_A_matrix_from_cloud_indices_vectorised(PointStructure* myPointStruct, double* A, int cloud_index) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    int mpn = m+n;
    
    int* cloud = (int*) malloc (myPointStruct->num_cloud_points * sizeof(int));
    for (int i = 0; i<m; i++)
        cloud[i] = myPointStruct->cloud_index[cloud_index*m +i];
        
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};
    // Parallelize the outer loop and the inner loop for the phi matrix with PHS RBF function
    for (int i = 0; i < m; i++) {
        A[i*mpn+i] = 0; // Diagonal element
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = i + 1; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            A[i*mpn+j] = calculate_phs_rbf(pt1, pt2, parameters.phs_degree, parameters.dimension);
            A[j*mpn+i] = A[i*mpn+j];
        }
    }

    // Parallelize the loop for adding polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A[i*mpn +j + m] = 1;
            if (myPointStruct->pow_x[j] != 0)
                A[i*mpn +j + m] *= pow(myPointStruct->x[cloud[i]] - seed_pt[0], myPointStruct->pow_x[j]);
            if (myPointStruct->pow_y[j] != 0)
                A[i*mpn +j + m] *= pow(myPointStruct->y[cloud[i]] - seed_pt[1], myPointStruct->pow_y[j]);
            if (myPointStruct->pow_z[j] != 0)
                A[i*mpn +j + m] *= pow(myPointStruct->z[cloud[i]] - seed_pt[2], myPointStruct->pow_z[j]);

            A[(j + m)*mpn + i] = A[i*mpn +j + m];
        }
    }
    free(cloud);
}

void gradx_matrix_vectorised(PointStructure* myPointStruct, double* grad, int cloud_index) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    int mpn = m+n;
    
    int* cloud = (int*) malloc (myPointStruct->num_cloud_points * sizeof(int));
    for (int i = 0; i<m; i++)
        cloud[i] = myPointStruct->cloud_index[cloud_index*m +i];
        
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i*mpn +j] = parameters.phs_degree * (pt1[0] - pt2[0]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree - 2, parameters.dimension);
        }
    }

    // Parallelizing the polynomial terms loop
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        // grad[i*mpn +m] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_x[j] == 0) {
                grad[i*mpn +j + m] = 0;
            } else {
                grad[i*mpn +j + m] = myPointStruct->pow_x[j] * pow(pt1[0] - seed_pt[0], myPointStruct->pow_x[j] - 1);
                grad[i*mpn +j + m] *= pow(pt1[1] - seed_pt[1], myPointStruct->pow_y[j]);
                grad[i*mpn +j + m] *= pow(pt1[2] - seed_pt[2], myPointStruct->pow_z[j]);
            }
            grad[(j + m)*mpn +i] = grad[i*mpn +j + m];
        }
    }
    free(cloud);
}

void grady_matrix_vectorised(PointStructure* myPointStruct, double* grad, int cloud_index) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    int mpn = m+n;
    
    int* cloud = (int*) malloc (myPointStruct->num_cloud_points * sizeof(int));
    for (int i = 0; i<m; i++)
        cloud[i] = myPointStruct->cloud_index[cloud_index*m +i];

    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]],myPointStruct->y[cloud[0]],myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) { 
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i*mpn +j] = parameters.phs_degree*(pt1[1]-pt2[1]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree-2, parameters.dimension);
        }
    }

    // add polynomial terms
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]]; 
        // grad[i*mpn +0] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_y[j] == 0)
                grad[i*mpn +j+m] = 0;
            else{
                grad[i*mpn +j+m] = myPointStruct->pow_y[j]*pow(pt1[1]-seed_pt[1], myPointStruct->pow_y[j]-1);
                grad[i*mpn +j+m] *= pow(pt1[0]-seed_pt[0], myPointStruct->pow_x[j]);
                grad[i*mpn +j+m] *= pow(pt1[2]-seed_pt[2], myPointStruct->pow_z[j]);
            }
            grad[(j+m)*mpn +i] = grad[i*mpn +j+m];
            // printf("grad[%d][%d] = %lf\n", i, j+m, grad[i][j+m]);
        }
    }
    free(cloud);
}

void gradz_matrix_vectorised(PointStructure* myPointStruct, double* grad, int cloud_index) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    int mpn = m+n;
    
    int* cloud = (int*) malloc (myPointStruct->num_cloud_points * sizeof(int));
    for (int i = 0; i<m; i++)
        cloud[i] = myPointStruct->cloud_index[cloud_index*m +i];
    
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]],myPointStruct->y[cloud[0]],myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++){ 
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i*mpn +j] = parameters.phs_degree*(pt1[2]-pt2[2]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree-2, parameters.dimension);
        }
    }

    // add polynomial terms
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        // grad[i*mpn +0] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_z[j] == 0)
                grad[i*mpn +j+m] = 0;
            else{
                grad[i*mpn +j+m] = myPointStruct->pow_z[j]*pow(pt1[2]-seed_pt[2], myPointStruct->pow_z[j]-1);
                grad[i*mpn +j+m] *= pow(pt1[1]-seed_pt[1], myPointStruct->pow_y[j]);
                grad[i*mpn +j+m] *= pow(pt1[0]-seed_pt[0], myPointStruct->pow_x[j]);
            }
            grad[(j+m)*mpn +i] = grad[i*mpn +j+m];
        }
    }
    free(cloud);
}

void laplacian_matrix_vectorised(PointStructure* myPointStruct, double* lap, int cloud_index) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    int mpn = m+n;
    
    int* cloud = (int*) malloc (myPointStruct->num_cloud_points * sizeof(int));
    for (int i = 0; i<m; i++)
        cloud[i] = myPointStruct->cloud_index[cloud_index*m +i];

    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};

    // Parallelizing the outer loop for calculating the Laplacian matrix
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            lap[i*mpn +j] = parameters.phs_degree * (parameters.phs_degree + parameters.dimension -2) *
                         calculate_phs_rbf(pt1, pt2, parameters.phs_degree - 2, parameters.dimension);
        }
    }

    // Adding polynomial terms
    double dtempx, dtempy, dtempz;
    
    if (parameters.dimension == 3) {
        for (int i = 0; i < m; i++) {
            pt1[0] = myPointStruct->x[cloud[i]] - seed_pt[0]; 
            pt1[1] = myPointStruct->y[cloud[i]] - seed_pt[1]; 
            pt1[2] = myPointStruct->z[cloud[i]] - seed_pt[2];
            for (int j = 0; j < n; j++) {      
                dtempx = 0; dtempy = 0; dtempz = 0;
                if (myPointStruct->pow_x[j] > 1) {  
                    dtempx =  myPointStruct->pow_x[j] * (myPointStruct->pow_x[j] - 1) *
                              pow(pt1[0], myPointStruct->pow_x[j] - 2);
                    dtempx *= pow(pt1[1], myPointStruct->pow_y[j]);
                    dtempx *= pow(pt1[2], myPointStruct->pow_z[j]);
                }
                if (myPointStruct->pow_y[j] > 1) {  
                    dtempy =  myPointStruct->pow_y[j] * (myPointStruct->pow_y[j] - 1) *
                              pow(pt1[1], myPointStruct->pow_y[j] - 2);
                    dtempy *= pow(pt1[0], myPointStruct->pow_x[j]);
                    dtempy *= pow(pt1[2], myPointStruct->pow_z[j]);
                }
                if (myPointStruct->pow_z[j] > 1) {  
                    dtempz =  myPointStruct->pow_z[j] * (myPointStruct->pow_z[j] - 1) *
                              pow(pt1[2], myPointStruct->pow_z[j] - 2);
                    dtempz *= pow(pt1[0], myPointStruct->pow_x[j]);
                    dtempz *= pow(pt1[1], myPointStruct->pow_y[j]);
                }
                lap[i*mpn +j + m] = dtempx + dtempy + dtempz;
                lap[(j + m)*mpn +i] = lap[i*mpn +j + m];
            }
        }
    }

    if (parameters.dimension == 2) {
        for (int i = 0; i < m; i++) {
            pt1[0] = myPointStruct->x[cloud[i]] - seed_pt[0]; 
            pt1[1] = myPointStruct->y[cloud[i]] - seed_pt[1]; 
            pt1[2] = myPointStruct->z[cloud[i]] - seed_pt[2];
            for (int j = 0; j < n; j++) {        
                dtempx = 0; dtempy = 0;
                // if (myPointStruct->pow_x[j] > 1) {  
                //     dtempx =  myPointStruct->pow_x[j] * (myPointStruct->pow_x[j] - 1) *
                //               pow(pt1[0], myPointStruct->pow_x[j] - 2);
                //     dtempx *= pow(pt1[1], myPointStruct->pow_y[j]);
                // }
                // else if (myPointStruct->pow_y[j] > 1) {  
                //     dtempy =  myPointStruct->pow_y[j] * (myPointStruct->pow_y[j] - 1) *
                //               pow(pt1[1], myPointStruct->pow_y[j] - 2);
                //     dtempy *= pow(pt1[0], myPointStruct->pow_x[j]);
                // }
                if (myPointStruct->pow_x[j] > 1) {  
                    dtempx = myPointStruct->pow_x[j] * (myPointStruct->pow_x[j] - 1) *
                            pow(pt1[0], myPointStruct->pow_x[j] - 2) *
                            pow(pt1[1], myPointStruct->pow_y[j]);
                }
                if (myPointStruct->pow_y[j] > 1) {  // NOT else if
                    dtempy = myPointStruct->pow_y[j] * (myPointStruct->pow_y[j] - 1) *
                            pow(pt1[1], myPointStruct->pow_y[j] - 2) *
                            pow(pt1[0], myPointStruct->pow_x[j]);
                }
                lap[i*mpn + j + m] = dtempx + dtempy;
                // lap[i*mpn +j + m] = dtempx + dtempy;
                lap[(j + m)*mpn +i] = lap[i*mpn +j + m];
            }
        }
    }
    free(cloud);
}

// Following function creates the full derivative matrices for the mesh

void create_full_gradx_matrix_vectorised(PointStructure* myPointStruct) {
    double *grad, *A, *A_inv, *B1;
    int n = myPointStruct->num_cloud_points;
    int m = n + myPointStruct->num_poly_terms;
    
    // Allocate memory for matrices
    grad = create_matrix_vectorised(m, m);
    A = create_matrix_vectorised(m, m);
    A_inv = create_matrix_vectorised(m, m);
    B1 = create_matrix_vectorised(m, m);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        // Sequential parts remain within the loop
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct, A, myPointStruct->cloud_index[i*n]);
        gradx_matrix_vectorised(myPointStruct, grad, myPointStruct->cloud_index[i*n]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m);
        multiply_matrices_vectorised(grad, A_inv, B1, m, m, m);
        
        // Parallelize the inner loop (assignment) if it's large enough
        for (int j = 0; j < myPointStruct->num_cloud_points; j++) 
            myPointStruct->Dx[i*n+j] = B1[j];
    }

    // Free matrices
    free(grad);
    free(A);
    free(A_inv);
    free(B1);
}

void create_full_grady_matrix_vectorised(PointStructure* myPointStruct) {
    double *grad, *A, *A_inv, *B1;
    int n = myPointStruct->num_cloud_points;
    int m = n + myPointStruct->num_poly_terms;
    
    // Allocate memory for matrices
    grad = create_matrix_vectorised(m, m);
    A = create_matrix_vectorised(m, m);
    A_inv = create_matrix_vectorised(m, m);
    B1 = create_matrix_vectorised(m, m);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct, A, myPointStruct->cloud_index[i*n]);
        grady_matrix_vectorised(myPointStruct, grad, myPointStruct->cloud_index[i*n]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m);
        multiply_matrices_vectorised(grad, A_inv, B1, m, m, m);
        for (int j = 0; j < myPointStruct->num_cloud_points; j++) 
            myPointStruct->Dy[i*n + j] = B1[j];
    }
    free(grad);
    free(A);
    free(A_inv);
    free(B1);
}

void create_full_gradz_matrix_vectorised(PointStructure* myPointStruct) {
    double *grad, *A, *A_inv, *B1;
    int n = myPointStruct->num_cloud_points;
    int m = n + myPointStruct->num_poly_terms;
    
    // Allocate memory for matrices
    grad = create_matrix_vectorised(m, m);
    A = create_matrix_vectorised(m, m);
    A_inv = create_matrix_vectorised(m, m);
    B1 = create_matrix_vectorised(m, m);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct, A, myPointStruct->cloud_index[i*n]);
        gradz_matrix_vectorised(myPointStruct, grad, myPointStruct->cloud_index[i*n]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m);
        multiply_matrices_vectorised(grad, A_inv, B1, m, m, m);
        for (int j = 0; j < myPointStruct->num_cloud_points; j++)
            myPointStruct->Dz[i*n +j] = B1[j];
    }
    free(grad);
    free(A);
    free(A_inv);
    free(B1);
}

void create_full_laplacian_matrix_vectorised(PointStructure* myPointStruct) {
    double *lap, *A, *A_inv, *B1;
    int n = myPointStruct->num_cloud_points;
    int m = n + myPointStruct->num_poly_terms;
    
    // Allocate memory for matrices
    lap = create_matrix_vectorised(m, m);
    A = create_matrix_vectorised(m, m);
    A_inv = create_matrix_vectorised(m, m);
    B1 = create_matrix_vectorised(m, m);
   
    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct, A, myPointStruct->cloud_index[i*n]);
        laplacian_matrix_vectorised(myPointStruct, lap, myPointStruct->cloud_index[i*n]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m);
        multiply_matrices_vectorised(lap, A_inv, B1, m,m,m);
        for (int j = 0; j < myPointStruct->num_cloud_points; j++) 
            myPointStruct->lap[i*n +j] = B1[j];
    }
    free(lap);
    free(A);
    free(A_inv);
    free(B1);
}

void create_laplacian_Poisson_vectorised(PointStructure* myPointStruct) {
    // Parallelize the outer loop with OpenACC
    int n = myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++){ 
        if (myPointStruct->boundary_tag[i]){
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET || myPointStruct->node_bc[i].type == BC_WALL) // Neumann BC
                for (int j = 0; j < n; j++) {
                    myPointStruct->lap_Poison[i*n +j] = myPointStruct->Dx[i*n +j] * myPointStruct->x_normal[i] + myPointStruct->Dy[i*n +j] * myPointStruct->y_normal[i];
                    if (parameters.dimension == 3)
                        myPointStruct->lap_Poison[i*n +j] += myPointStruct->Dz[i*n +j] * myPointStruct->z_normal[i];
                }
            else{ // Dirichlet BC for pressure
                for (int j = 1; j < n; j++) 
                    myPointStruct->lap_Poison[i*n +j] = 0;
                myPointStruct->lap_Poison[i*n +0] = 1;
            }
        }
        else{
            for (int j = 0; j < n; j++) 
                myPointStruct->lap_Poison[i*n +j] = myPointStruct->lap[i*n +j];
        }
    }
}

void create_derivative_matrices_vectorised(PointStructure* myPointStruct){
    int m = myPointStruct->num_nodes;
    int n = myPointStruct->num_cloud_points;
    myPointStruct->Dx = create_matrix_vectorised(m,n);
    myPointStruct->Dy = create_matrix_vectorised(m,n);
    myPointStruct->lap = create_matrix_vectorised(m,n);
    if (parameters.dimension == 3){
        myPointStruct->Dz = create_matrix_vectorised(m,n);
        create_full_gradx_matrix_vectorised(myPointStruct);
        create_full_grady_matrix_vectorised(myPointStruct);
        create_full_gradz_matrix_vectorised(myPointStruct);
        create_full_laplacian_matrix_vectorised(myPointStruct);
    }
    if (parameters.dimension == 2){
        create_full_gradx_matrix_vectorised(myPointStruct);
        create_full_grady_matrix_vectorised(myPointStruct);
        create_full_laplacian_matrix_vectorised(myPointStruct);
    }
    create_laplacian_for_Poisson_equation_vectorised(myPointStruct);
}

void create_laplacian_for_Poisson_equation_vectorised(PointStructure* myPointStruct){
    int m = myPointStruct->num_nodes;
    int n = myPointStruct->num_cloud_points;
    myPointStruct->lap_Poison = create_matrix_vectorised(m,n);
    create_laplacian_Poisson_vectorised(myPointStruct);
}