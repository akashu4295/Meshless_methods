// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef FRACTIONAL_STEP_SOLVERS_H
#define FRACTIONAL_STEP_SOLVERS_H

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "structures.h"
#include "general_functions.h"
#include "mesh_functions.h"
#include "mat_lib.h"
#include "kdtree_functions.h"
#include "rbf.h"

/////////////////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////////////////

// Fractional Step Explicit Solver Modules

double fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_intermediate_velocity(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_mass_residual(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_boundary_dpdn(PointStructure* myPointStruct, FieldVariables* field);
void FS_multigrid_Poisson_solver(PointStructure* myPointStruct, FieldVariables* field);
void FS_relaxation(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_residuals(PointStructure* myPointStruct, FieldVariables* field);
void FS_restrict_residuals(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_prolongate_corrections(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_update_velocity(PointStructure* myPointStruct, FieldVariables* field);
void FS_update_boundary_pressure(PointStructure* myPointStruct, FieldVariables* field);

/////////////////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Fractional Step Explicit Solver Modules
/////////////////////////////////////////////////////////////////////////////
double fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* field){   
    double steady_state_error = 0.0;
    # pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        if (parameters.dimension == 3)
            field[0].w_old[i] = field[0].w[i];
    }

    FS_calculate_intermediate_velocity(&myPointStruct[0], &field[0]);
    FS_calculate_mass_residual(&myPointStruct[0], &field[0]);
    FS_calculate_boundary_dpdn(&myPointStruct[0], &field[0]);
    FS_multigrid_Poisson_solver(myPointStruct, field);
    FS_update_velocity(&myPointStruct[0], &field[0]);

    # pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        steady_state_error += (field[0].u[i] - field[0].u_old[i])*(field[0].u[i] - field[0].u_old[i]) + (field[0].v[i] - field[0].v_old[i])*(field[0].v[i] - field[0].v_old[i]);
        if (parameters.dimension == 3)
            steady_state_error += (field[0].w[i] - field[0].w_old[i])*(field[0].w[i] - field[0].w_old[i]);
    }
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}

void FS_calculate_intermediate_velocity(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->u, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->u, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->lap, field->u, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu *field->dpdn[i]);

// y-momentum
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->v, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->v, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->lap, field->v, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);

//  z-momentum
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->w, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->w, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_gpu(myPointStruct->lap, field->w, field->dpdn, myPointStruct->cloud_index, myPointStruct->num_nodes,myPointStruct->num_cloud_points);
        
        # pragma acc parallel loop present(field, parameters, myPointStruct)
        for (int i = 0; i < myPointStruct->num_nodes; i++)
            field->w_new[i] = field->w[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);
    }
    
    # pragma acc update host(field[0].u_new[:num_nodes], field[0].v_new[:num_nodes])
}

void FS_calculate_mass_residual(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->u_new, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->v_new, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->w_new, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct, temp1, temp2, temp3)
    for (int i = 0; i < num_nodes; i++)
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
}

void FS_calculate_boundary_dpdn(PointStructure* myPointStruct, FieldVariables* field){
    double dpdx, dpdy, dpdz;
    // Copy only boundary normals to gpu
    # pragma acc parallel loop present(field, parameters, myPointStruct, zeros[:myPointStruct->num_nodes])
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i]){
            dpdx = (field->u_new[i] - field->u[i]) * parameters.rho/parameters.dt; 
            dpdy = (field->v_new[i] - field->v[i]) * parameters.rho/parameters.dt;
            field->dpdn[i] = dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i];
            if (parameters.dimension == 3){
                dpdz = (field->w_new[i] - field->w[i]) * parameters.rho/parameters.dt; 
                field->dpdn[i] += dpdz*myPointStruct->z_normal[i];
            }
        }
    }
}

void FS_multigrid_Poisson_solver(PointStructure* myPointStruct, FieldVariables* field){
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            FS_relaxation(&myPointStruct[ilev], &field[ilev]);
            FS_calculate_residuals(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                FS_restrict_residuals(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
                }
            }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            FS_prolongate_corrections(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1) 
                FS_relaxation(&myPointStruct[ilev-1], &field[ilev-1]);
        }   
        // update_boundary_pressure(&myPointStruct[0], &field[0]);
    } 
}

void FS_relaxation(PointStructure* mypointstruct, FieldVariables* field){
    # pragma acc loop
    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            double sum = 0.0;
            # pragma acc loop reduction(+:sum)
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                sum += mypointstruct->lap_Poison[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            }
            field->p[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[i][0]) + (1-parameters.omega)*field->p[i]; 
        }
        double pref = field->p[0];
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->p[i] = field->p[i] - pref;
        }
        
        # pragma acc parallel loop present(field, mypointstruct)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            field->res[i]=field->source[i];
            for (int j = 0; j < mypointstruct->num_cloud_points; j++){
                field->res[i] = field->res[i] - mypointstruct->lap[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            }
        }
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));
    
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->p, zeros, mypointstruct->num_nodes));
    }
}

void FS_restrict_residuals(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            double results = 0.0;
            int i_restr = mypointStruct_c->restriction_points[i];
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_f->num_cloud_points; j++){
                results += mypointStruct_c->restr_mat[i][j] * field_f->res[mypointStruct_f->cloud_index[i_restr][j]];
            }
            field_c->source[i] = results;
        }
    }
}

void FS_prolongate_corrections(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol = mypointStruct_f->prolongation_points[i];
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_c->num_cloud_points; j++){
                results += mypointStruct_f->prol_mat[i][j] * field_c->p[mypointStruct_c->cloud_index[i_prol][j]];
            }
        field_f->p[i] = field_f->p[i] + results;
        }
    }
    # pragma acc parallel loop present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            field_c->p[i] = 0.0;
        }
    }
}

void FS_calculate_residuals(PointStructure* mypointStruct, FieldVariables* field){
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            for (int j = 0; j < mypointStruct->num_cloud_points; j++){
                sum += mypointStruct->lap[i][j]*field->p[mypointStruct->cloud_index[i][j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}

void FS_update_velocity(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    }
    // Update Interior nodes
    # pragma acc parallel loop present(field, parameters, myPointStruct, dpdx, dpdy, dpdz)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            if (parameters.dimension == 3){
                field->w[i] = field->w_new[i] - parameters.dt * field->dpdz[i]/parameters.rho;
            }
        }
    }
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    double sum;
    sum = 0.0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}


void FS_update_boundary_pressure(PointStructure* mypointstruct, FieldVariables* field){
    double sumx, sumy, sumz, Ap, term;

    if (parameters.neumann_flag_boundary){
    # pragma acc parallel loop present(field, parameters, mypointstruct)
    for (int i = 0; i < mypointstruct->num_boundary_nodes; i++)
    {
        sumx = 0; sumy = 0; sumz = 0; Ap = 0;
        # pragma acc loop reduction(+:sumx, sumy, sumz, Ap)
        for (int j = 1; j < mypointstruct->num_cloud_points; j++){
            sumx += mypointstruct->Dx[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            sumy += mypointstruct->Dy[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            if (parameters.dimension == 3){
                sumz += mypointstruct->Dz[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            }
        }
        Ap += mypointstruct->Dx[i][0]*mypointstruct->x_normal[i];
        Ap += mypointstruct->Dy[i][0]*mypointstruct->y_normal[i];
        if (parameters.dimension == 3){
            Ap += mypointstruct->Dz[i][0]*mypointstruct->z_normal[i];
        }
        term = (field->dpdn[i]-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i] -sumz*mypointstruct->z_normal[i])/Ap;
        field->p[i] = field->p[i] * 0.5 + 0.5 * term;
    }
    }
}

#endif
