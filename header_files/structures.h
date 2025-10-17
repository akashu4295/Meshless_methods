// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////////
// Structures
//////////////////////////////////////////////////////////////////////////////

// Structure to represent the common parameters for all grids
struct parameters
    {
    short dimension; //dimension of the problem, default is 3
    short poly_degree; //degree of the polynomial basis functions
    short phs_degree; //degree of the PHS basis functions
    short cloud_size_multiplier; //multiplier for the number of cloud points
    short test;
    double courant_number; // Courant number for the time step
    double steady_state_tolerance; // Tolerance for steady state convergence
    double poisson_solver_tolerance; // Tolerance for Poisson solver
    short num_vcycles; // Number of V-cycles
    short num_relax; // Number of relaxation steps
    int num_time_steps; // Number of time steps
    short num_levels; // Number of levels in the multigrid
    int write_interval; // Interval for writing the data
    bool neumann_flag_boundary; // Neumann boundary condition flag
    float omega; // relaxation parameter
    double dt; // time step
    short iter_momentum; // number of iterations for momentum equation
    short iter_timple; // number of iterations for timple time stepping
    double rho; // density
    double mu; // dynamic viscosity
    double Re; // Reynolds number
    double facRe;//Reynolds number factor for defect correction
    double facdt;//dual time step factor
    double nu; // kinematic viscosity
    bool fractional_step; // fractional step flag
    }parameters; // child created only once and used globally

// Structure to represent the mesh data
typedef struct PointStructure {   
    char mesh_filename[50]; // name of the mesh file
    int num_nodes; // number of nodes
    int num_corners; // number of corner nodes
    int num_boundary_nodes; // number of boundary nodes
    int num_elem; // number of elements
    double d_avg; // average distance between nodes
    short num_poly_terms; //number of polynomial terms  //////NUM_POLY_TERMS
    short num_cloud_points; //number of cloud points in the domain
    short poly_degree; //degree of the polynomial basis functions
    double* x;  // x coordinates of the nodes
    double* y;  // y coordinates of the nodes
    double* z;  // z coordinates of the nodes
    int* point_index; // index of the point in the original dataset
    double* x_normal; // x component of the normal vector
    double* y_normal; // y component of the normal vector
    double* z_normal; // z component of the normal vector
    bool* boundary_tag; // Boolean tag for boundary nodes
    bool* corner_tag; // Boolean tag for corner nodes
    short* pow_x; //power of x in the polynomial basis functions
    short* pow_y; //power of y in the polynomial basis functions
    short* pow_z; //power of z in the polynomial basis functions
    int** cloud_index; // indices of the INTERPOLATION cloud points 
    int* prolongation_points; // indices of the PROLONGATION cloud points
    int* restriction_points; // indices of the RESTRICTION cloud points
    double** Dx; // x-derivative matrix
    double** Dy; // y-derivative matrix
    double** Dz; // z-derivative matrix
    double** lap; // laplacian matrix
    double** lap_Poison; // laplacian matrix for Poisson solver
    double** restr_mat;// restriction matrix
    double** prol_mat;// prolongation matrix
    int** prolongation_points2; // indices of the PROLONGATION cloud points
    int** restriction_points2; // indices of the RESTRICTION cloud points
}PointStructure;

typedef struct FieldVariables {
    double* u; // x-component of the velocity field
    double* v; // y-component of the velocity field
    double* w; // z-component of the velocity field
    double* p; // pressure
    double* u_new; // x-component of the intermediate velocity field
    double* v_new; // y-component of the intermediate velocity field
    double* w_new; // z-component of the intermediate velocity field
    double* u_old; // x-component of the old velocity field
    double* v_old; // y-component of the old velocity field
    double* w_old; // z-component of the old velocity field
    double* p_old; // pressure field
    double* pprime; // pressure field
    double* res; // residual field
    double* source; // source field
    double* T; // temperature field
    double* dpdn; // dpdn field
    double* dpdx; // dpdx field
    double* dpdy; // dpdy field
    double* dpdz; // dpdz field
    double rho; // density field
    double mu; // dynamic viscosity field
    
    double nu; // kinematic viscosity field
    double Re; // Reynolds number field
    //double facRe; // Re factor for defect correction
} FieldVariables;

// Structure to represent a point in 3-dimensional space
typedef struct Point {
    double coords[3];
    int index; // Index of the point in the original dataset
} Point;

void AllocateMemoryPointStructure(PointStructure* myPointStruct, int nodes) {
    myPointStruct->x_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->x = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z = (double*)malloc(nodes * sizeof(double));
    myPointStruct->point_index = (int*)malloc(nodes * sizeof(int));
    myPointStruct->boundary_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->corner_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->num_nodes = nodes;
    myPointStruct->num_elem = 0;
}

void AllocateMemoryFieldVariables(FieldVariables** field, PointStructure* myPointStruct, int num_levels) {
    *field = (FieldVariables*) malloc(num_levels * sizeof(FieldVariables));
    for (int ii = 0; ii < num_levels; ii++) {
        (*field)[ii].u = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].u_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].pprime = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].p = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].u_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].p_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].res = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].source = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdn = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdx = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdy = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdz = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].T = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
            (*field)[ii].u[i] = 0.0;
            (*field)[ii].v[i] = 0.0;
            (*field)[ii].w[i] = 0.0;
            (*field)[ii].u_new[i] = 0.0;
            (*field)[ii].v_new[i] = 0.0;
            (*field)[ii].w_new[i] = 0.0;
            (*field)[ii].u_old[i] = 0.0;
            (*field)[ii].v_old[i] = 0.0;
            (*field)[ii].w_old[i] = 0.0;
            (*field)[ii].pprime[i] = 0.0;
            (*field)[ii].p_old[i] = 0.0;
            (*field)[ii].p[i] = 0.0;
            (*field)[ii].res[i] = 0.0;
            (*field)[ii].source[i] = 0.0;
            (*field)[ii].dpdn[i] = 0.0;
            (*field)[ii].dpdx[i] = 0.0;
            (*field)[ii].dpdy[i] = 0.0;
            (*field)[ii].dpdz[i] = 0.0;
            (*field)[ii].T[i] = 0.0;
        }
    }
}

void free_PointStructure(PointStructure* myPointStruct, int num_levels) {
    for (int i = 0; i<num_levels; i++){
        free(myPointStruct[i].x);
        free(myPointStruct[i].y);
        free(myPointStruct[i].z);
        free(myPointStruct[i].point_index);
        free(myPointStruct[i].x_normal);
        free(myPointStruct[i].y_normal);
        free(myPointStruct[i].z_normal);
        free(myPointStruct[i].boundary_tag);
        free(myPointStruct[i].corner_tag);
        free(myPointStruct[i].pow_x);
        free(myPointStruct[i].pow_y);
        free(myPointStruct[i].pow_z);
        for (int j = 0; j < myPointStruct[i].num_nodes; j++){
            free(myPointStruct[i].cloud_index[j]);
            free(myPointStruct[i].Dx[j]);
            free(myPointStruct[i].Dy[j]);
            if (parameters.dimension == 3)
                free(myPointStruct[i].Dz[j]);
            free(myPointStruct[i].lap[j]);
            if (parameters.num_levels > 1){
                free(myPointStruct[i].restr_mat[j]);
                free(myPointStruct[i].prol_mat[j]);
            }
        }
        free(myPointStruct[i].cloud_index);
        free(myPointStruct[i].Dx);
        free(myPointStruct[i].Dy);
        if (parameters.dimension == 3)
            free(myPointStruct[i].Dz);
        free(myPointStruct[i].lap);            
        free(myPointStruct[i].lap_Poison);
        if (parameters.num_levels > 1){
            free(myPointStruct[i].restr_mat);
            free(myPointStruct[i].prol_mat);
            free(myPointStruct[i].prolongation_points);
            free(myPointStruct[i].restriction_points);
        }
    }
    free(myPointStruct);
}

void free_field(FieldVariables* field, int num_levels) {
    for (int i = 0; i < num_levels; i++){
        free(field[i].u);
        free(field[i].v);
        free(field[i].w);
        free(field[i].pprime);
        free(field[i].T);
        free(field[i].res);
        free(field[i].source);
        free(field[i].dpdn);
        free(field[i].u_new);
        free(field[i].v_new);
        free(field[i].w_new);
    }
    free(field);
}

#endif
