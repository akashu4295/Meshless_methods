// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "structures.h"
#include "functions.h"


/////////////////////////////////////////////////////////////////////////////
// Fractional Step Explicit Solver Modules
/////////////////////////////////////////////////////////////////////////////

double fractional_step_explicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){   
    # pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        field[0].w_old[i] = field[0].w[i];
        field[0].p_old[i] = field[0].p[i];
    }

    #pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        FS_calculate_intermediate_velocity_vectorised(&myPointStruct[0], &field[0]);
        FS_calculate_mass_residual_vectorised(&myPointStruct[0], &field[0]);
        FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
        FS_update_velocity_vectorised(&myPointStruct[0], &field[0]);
    }
    
    double steady_state_error_par = 0.0;
    # pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error_par)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        double dw = field[0].w[i] - field[0].w_old[i];
        steady_state_error_par += fabs(du) + fabs(dv) + fabs(dw);
    }
    return (steady_state_error_par/(parameters.dimension * myPointStruct[0].num_nodes * parameters.dt));
}

double fractional_step_explicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{   
    #pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        field[0].p_old[i] = field[0].p[i];
    }

    #pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        FS_calculate_intermediate_velocity_vectorised_2d(myPointStruct, field);
        FS_calculate_mass_residual_vectorised_2d(myPointStruct, field);
        FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
        FS_update_velocity_vectorised_2d(myPointStruct, field);
    }

    double steady_state_error_par = 0.0;    
    #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error_par) 
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        steady_state_error_par += fabs(du) + fabs(dv);
    }
    
    return (steady_state_error_par/(parameters.dimension * myPointStruct[0].num_nodes * parameters.dt));
}

void FS_calculate_intermediate_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    # pragma acc data present(field, myPointStruct, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->u, field->dudy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->u, field->dudz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->lap, field->u, field->lapu, myPointStruct->cloud_index, num_nodes, num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->v, field->dvdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 5);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 6);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->v, field->dvdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 7);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->lap, field->v, field->lapv, myPointStruct->cloud_index, num_nodes, num_cloud_points, 8);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->w, field->dwdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 9);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->w, field->dwdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 10);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w, field->dwdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 11);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->lap, field->w, field->lapw, myPointStruct->cloud_index, num_nodes, num_cloud_points, 12);    
        # pragma acc wait(1,2,3,4,5,6,7,8,9,10,11,12)
    }
    # pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dudx[i] + field->v[i] * field->dudy[i] + field->w[i] * field->dudz[i] - parameters.nu * field->lapu[i]);
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dvdx[i] + field->v[i] * field->dvdy[i] + field->w[i] * field->dvdz[i] - parameters.nu * field->lapv[i]);
        field->w_new[i] = field->w[i] - parameters.dt * (field->u[i] * field->dwdx[i] + field->v[i] * field->dwdy[i] + field->w[i] * field->dwdz[i] - parameters.nu * field->lapw[i]);
    }
    // /* ---- ENFORCE BOUNDARY CONDITIONS ON u* ---- */
    // # pragma acc parallel loop gang vector default(present)
    // for (int i = 0; i < myPointStruct->num_nodes; i++){
    //     if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]){
    //         if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET 
    //                 || myPointStruct->node_bc[i].type == BC_WALL 
    //                 || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET){
    //             field->u_new[i] = myPointStruct->node_bc[i].u;
    //             field->v_new[i] = myPointStruct->node_bc[i].v;
    //             field->w_new[i] = myPointStruct->node_bc[i].w;
    //         }
    //     }
    // }
}

// # pragma acc routine
void FS_calculate_intermediate_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    #pragma acc data present(field, myPointStruct, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->u, field->dudy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->lap, field->u, field->lapu, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->v, field->dvdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 5);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->lap, field->v, field->lapv, myPointStruct->cloud_index, num_nodes, num_cloud_points, 6);
        # pragma acc wait(1,2,3,4,5,6)
    }
    # pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dudx[i] + field->v[i] * field->dudy[i] - parameters.nu *field->lapu[i]);
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dvdx[i] + field->v[i] * field->dvdy[i] - parameters.nu *field->lapv[i]);
    }
    #pragma acc data present(field, myPointStruct, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p_old, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p_old, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        # pragma acc wait(1,2)
    }
    # pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i] 
                        && !myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            field->u_new[i] = field->u[i] + parameters.dt * (field->dpdx[i]);
            field->v_new[i] = field->v[i] + parameters.dt * (field->dpdy[i]);
        }
    }
}

void FS_calculate_mass_residual_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u_new, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v_new, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w_new, field->dwdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
        # pragma acc wait(1,2,3)
    }
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) // skip corners
            continue;
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dudx[i]+field->dvdy[i]+field->dwdz[i])/parameters.dt;
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) // pressure outlet boundary nodes
            field->source[i] = myPointStruct->node_bc[i].p;
        else
            field->source[i] = parameters.rho*((field->u_new[i] - field->u[i])*myPointStruct->x_normal[i] + (field->v_new[i] - field->v[i])*myPointStruct->y_normal[i] + (field->w_new[i] - field->w[i])*myPointStruct->z_normal[i])/parameters.dt;
            // field->source[i] = 0;
    }
}

void FS_calculate_mass_residual_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u_new, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v_new, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        # pragma acc wait(1,2)
    }
    # pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) // skip corners
            continue;
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dudx[i]+field->dvdy[i])/parameters.dt;
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) // pressure outlet boundary nodes
            field->source[i] = myPointStruct->node_bc[i].p;
        else
            field->source[i] = parameters.rho*((field->u_new[i] - field->u[i])*myPointStruct->x_normal[i] + (field->v_new[i] - field->v[i])*myPointStruct->y_normal[i])/parameters.dt;
    }
}

void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            if (parameters.poisson_solver_type == 1)    
                FS_relaxation_vectorised_Jacobi(&myPointStruct[ilev], &field[ilev]);
            else if (parameters.poisson_solver_type == 2)
                FS_relaxation_vectorised_Gauss_Seidel(&myPointStruct[ilev], &field[ilev]);
            else if (parameters.poisson_solver_type == 3)
                FS_relaxation_vectorised_BiCGStab(&myPointStruct[ilev], field[ilev].source, field[ilev].p, parameters.num_relax, parameters.poisson_solver_tolerance);
            FS_calculate_residuals_vectorised(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                FS_restrict_residuals_vectorised(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }
        }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            FS_prolongate_corrections_vectorised(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1) {
                if (parameters.poisson_solver_type == 1)    
                    FS_relaxation_vectorised_Jacobi(&myPointStruct[ilev-1], &field[ilev-1]);
                else if (parameters.poisson_solver_type == 2)
                    FS_relaxation_vectorised_Gauss_Seidel(&myPointStruct[ilev-1], &field[ilev-1]);
                else if (parameters.poisson_solver_type == 3)
                    FS_relaxation_vectorised_BiCGStab(&myPointStruct[ilev-1], field[ilev-1].source, field[ilev-1].p, parameters.num_relax, parameters.poisson_solver_tolerance);   
            } 
        }
    } 
}

void FS_relaxation_vectorised_Jacobi(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;
    #pragma acc parallel loop present(field->p[:N], field->p_old[:N])
    for (int i = 0; i < N; i++) {
        field->p_old[i] = field->p[i];
    }
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector present(field->p[:N], field->p_old[:N], mypointstruct->cloud_index[:N*n], mypointstruct->lap_Poison[:N*n], parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 1; j < n; j++) {
                int idx = i*n + j;
                sum += mypointstruct->lap_Poison[idx] * field->p_old[mypointstruct->cloud_index[idx]];
            }
            field->p[i] = parameters.omega * (field->source[i] - sum) / mypointstruct->lap_Poison[i*n] 
                            + (1.0 - parameters.omega) * field->p_old[i];
        }
        // double mean = 0;
        // #pragma acc parallel loop present(field->p[:N], field->p_old[:N])
        // for (int i = 0; i < N; i++) {
        //     mean += field->p[i];
        // }
        // mean = mean/N;
        // #pragma acc parallel loop present(field->p[:N], field->p_old[:N])
        // for (int i = 0; i < N; i++) {
        //     field->p_old[i] = field->p[i]-mean;
        // }
    }
}

void FS_relaxation_vectorised_Gauss_Seidel(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;
    
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector_length(128) present(field, mypointstruct, parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            #pragma acc loop vector reduction(+:sum)
            for (int j = 1; j < n; j++) {
                sum += mypointstruct->lap_Poison[i*n + j] *
                       field->p[mypointstruct->cloud_index[i*n + j]];
            }
            field->p[i] = parameters.omega * ((field->source[i] - sum) / mypointstruct->lap_Poison[i*n])
                + (1 - parameters.omega) * field->p[i];
        }
    }
}

void FS_restrict_residuals_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            double results = 0.0;
            int i_restr = mypointStruct_c->restriction_points[i];
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_c->restr_mat[i*n + j] * field_f->res[mypointStruct_f->cloud_index[i_restr*n +j]];
            }
            field_c->source[i] = results;
        }
    }
}

void FS_prolongate_corrections_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_c->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol = mypointStruct_f->prolongation_points[i];
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_f->prol_mat[i*n + j] * field_c->p[mypointStruct_c->cloud_index[i_prol*n +j]];
            }
        field_f->p[i] = field_f->p[i] + results;
        }
    }
    # pragma acc parallel loop gang vector present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            field_c->p[i] = 0.0;
        }
    }
}

void FS_calculate_residuals_vectorised(PointStructure* mypointStruct, FieldVariables* field)
{
    int n = mypointStruct->num_cloud_points;
    double sum_res = 0.0;
    # pragma acc parallel loop gang vector present(field, mypointStruct) reduction(+:sum_res)
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap_Poison[i*n + j]*field->p[mypointStruct->cloud_index[i*n + j]];
            }
            field->res[i] = field->source[i] - sum;
            sum_res += fabs(field->res[i]);
        }
    }
    // printf("Poisson residual: %e\n", sum_res);
}

void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 3);
    }
    # pragma acc wait(1,2,3)
    // Update Interior nodes
    # pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            double sumux = 0.0, sumuy = 0.0, sumuz = 0.0;
            double sumvx = 0.0, sumvy = 0.0, sumvz = 0.0;
            double sumwx = 0.0, sumwy = 0.0, sumwz = 0.0;

            for (int j = 1; j < myPointStruct->num_cloud_points; j++){
                int k = i * myPointStruct->num_cloud_points + j;
                int neighbor = myPointStruct->cloud_index[k];
                
                // Collect all neighbor contributions for u, v, and w
                sumux -= myPointStruct->Dx[k] * field->u[neighbor];
                sumuy -= myPointStruct->Dy[k] * field->u[neighbor];
                sumuz -= myPointStruct->Dz[k] * field->u[neighbor];

                sumvx -= myPointStruct->Dx[k] * field->v[neighbor];
                sumvy -= myPointStruct->Dy[k] * field->v[neighbor];
                sumvz -= myPointStruct->Dz[k] * field->v[neighbor];

                sumwx -= myPointStruct->Dx[k] * field->w[neighbor];
                sumwy -= myPointStruct->Dy[k] * field->w[neighbor];
                sumwz -= myPointStruct->Dz[k] * field->w[neighbor];
            }

            double Ap = myPointStruct->Dx[i*myPointStruct->num_cloud_points] * myPointStruct->x_normal[i]
                    + myPointStruct->Dy[i*myPointStruct->num_cloud_points] * myPointStruct->y_normal[i]
                    + myPointStruct->Dz[i*myPointStruct->num_cloud_points] * myPointStruct->z_normal[i];

            // Ensure the velocity update projects correctly along the normals
            field->u[i] = (sumux * myPointStruct->x_normal[i] + sumuy * myPointStruct->y_normal[i] + sumuz * myPointStruct->z_normal[i]) / Ap;
            field->v[i] = (sumvx * myPointStruct->x_normal[i] + sumvy * myPointStruct->y_normal[i] + sumvz * myPointStruct->z_normal[i]) / Ap;
            field->w[i] = (sumwx * myPointStruct->x_normal[i] + sumwy * myPointStruct->y_normal[i] + sumwz * myPointStruct->z_normal[i]) / Ap;
        }
        else if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]){
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET 
                    || myPointStruct->node_bc[i].type == BC_WALL 
                    || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET){
                field->u[i] = myPointStruct->node_bc[i].u;
                field->v[i] = myPointStruct->node_bc[i].v;
                field->w[i] = myPointStruct->node_bc[i].w;
            }
        }
        else {//if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            field->w[i] = field->w_new[i] - parameters.dt * field->dpdz[i]/parameters.rho;
        }
    }
    
    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 5);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w, field->dwdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 6);
    }
    # pragma acc wait(4,5,6)
    double sum = 0.0;
    # pragma acc parallel loop gang vector default(present) reduction(+:sum)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        if (!myPointStruct->corner_tag[i])
            sum += parameters.rho*fabs(field->dudx[i]+field->dvdy[i]+field->dwdz[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}

void FS_update_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
    }
    # pragma acc wait(1,2)
    #pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            double sumux = 0.0, sumuy = 0.0, sumvx = 0.0, sumvy = 0.0;
            for (int j = 1; j < num_cloud_points; j++){
                int k = i*num_cloud_points + j;
                sumux -= myPointStruct->Dx[k] * field->u[myPointStruct->cloud_index[k]];
                sumuy -= myPointStruct->Dy[k] * field->u[myPointStruct->cloud_index[k]];
                sumvx -= myPointStruct->Dx[k] * field->v[myPointStruct->cloud_index[k]];
                sumvy -= myPointStruct->Dy[k] * field->v[myPointStruct->cloud_index[k]];
            }
            double Ap = myPointStruct->Dx[i*num_cloud_points] * myPointStruct->x_normal[i]
                        + myPointStruct->Dy[i*num_cloud_points] * myPointStruct->y_normal[i];
            field->u[i] = (sumux * myPointStruct->x_normal[i] + sumuy * myPointStruct->y_normal[i]) / Ap;
            field->v[i] = (sumvx * myPointStruct->x_normal[i] + sumvy * myPointStruct->y_normal[i]) / Ap;
        }
        else if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET 
                    || myPointStruct->node_bc[i].type == BC_WALL 
                    || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET){
            field->u[i] = myPointStruct->node_bc[i].u;
            field->v[i] = myPointStruct->node_bc[i].v;
        }
        else {//if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
        }
    }
    
    # pragma acc data present(myPointStruct, field, parameters)
    {   
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy,  myPointStruct->cloud_index, num_nodes, num_cloud_points, 4);
    }
    # pragma acc wait(3,4)
    
    double sum = 0.0;

    #pragma acc parallel loop gang vector reduction(+:sum) default(present)
    for (int i = 0; i < num_nodes; i++)
        if (!myPointStruct->corner_tag[i])
            sum += parameters.rho*fabs(field->dudx[i]+field->dvdy[i]);

    printf("Mass residual: %e\n", sum/num_nodes);
}

/* ============================================================
   Meshless Laplacian operator identical to Jacobi/GS
   ============================================================ */
static void apply_poisson_operator(PointStructure *ps, const double *x, double *Ax)
{
    int N = ps->num_nodes;
    int n = ps->num_cloud_points;

    #pragma acc parallel loop present(ps->lap_Poison[:N*n], ps->cloud_index[:N*n], x[:N], Ax[:N])
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        /* OFF-DIAGONALS ONLY (match Jacobi/GS) */
        for (int j = 1; j < n; j++) {
            int idx = i*n + j;
            sum += ps->lap_Poison[idx] * x[ps->cloud_index[idx]];
        }
        /* FULL operator */
        Ax[i] = ps->lap_Poison[i*n] * x[i] + sum;
    }
}

#define BICGSTAB_DEBUG 0


void FS_relaxation_vectorised_BiCGStab(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
{
    int N = ps->num_nodes;
    int n = ps->num_cloud_points;

    double *diag_inv = malloc(N * sizeof(double));
    double *r  = malloc(N * sizeof(double));
    double *r0 = malloc(N * sizeof(double));
    double *p  = malloc(N * sizeof(double));
    double *v  = malloc(N * sizeof(double));
    double *s  = malloc(N * sizeof(double));
    double *t  = malloc(N * sizeof(double));
    double *z  = malloc(N * sizeof(double));
    double *y  = malloc(N * sizeof(double));

    // Check if a Dirichlet pressure BC exists — if so, nullspace is already fixed
    int has_pressure_bc = 0;
    for (int i = 0; i < N; i++) {
        if (ps->node_bc[i].type == BC_PRESSURE_OUTLET) {
            has_pressure_bc = 1;
            break;
        }
    }

    #pragma acc data present(ps->lap_Poison[0:N*n], ps->cloud_index[0:N*n], \
                             ps->boundary_tag[0:N]) \
                     copyin(b[0:N]) \
                     copy(x[0:N]) \
                     create(diag_inv[0:N], r[0:N], r0[0:N], p[0:N], v[0:N], \
                            s[0:N], t[0:N], z[0:N], y[0:N])
    {
        /* Diagonal inverse (Jacobi preconditioner) */
        #pragma acc parallel loop present(diag_inv, ps->lap_Poison)
        for (int i = 0; i < N; i++) {
            double diag = ps->lap_Poison[i*n];
            diag_inv[i] = (fabs(diag) < 1e-14) ? 1.0 : 1.0 / diag;
        }

        /* Initial residual r = b - A*x */
        apply_poisson_operator(ps, x, v);

        #pragma acc parallel loop present(r, r0, p, v, b)
        for (int i = 0; i < N; i++) {
            r[i]  = b[i] - v[i];
            r0[i] = r[i];
            p[i]  = 0.0;
            v[i]  = 0.0;
        }

        double norm_b = 0.0;
        #pragma acc parallel loop reduction(+:norm_b) present(b)
        for (int i = 0; i < N; i++)
            norm_b += b[i] * b[i];
        norm_b = sqrt(norm_b);
        if (norm_b < 1e-14) norm_b = 1.0; // guard against zero RHS

        double norm_r0 = 0.0;
        #pragma acc parallel loop reduction(+:norm_r0) present(r)
        for (int i = 0; i < N; i++)
            norm_r0 += r[i] * r[i];
        norm_r0 = sqrt(norm_r0);

        double rho = 1.0, alpha = 1.0, omega_bicg = 1.0;

        for (int iter = 0; iter < max_iter; iter++) {

            double rho_new = 0.0;
            #pragma acc parallel loop reduction(+:rho_new) present(r0, r)
            for (int i = 0; i < N; i++)
                rho_new += r0[i] * r[i];

            if (fabs(rho_new) < 1e-30) {
                // Breakdown — restart with current residual as r0
                #pragma acc parallel loop present(r0, r)
                for (int i = 0; i < N; i++)
                    r0[i] = r[i];
                rho_new = norm_r0 * norm_r0;
                if (fabs(rho_new) < 1e-30) break;
            }

            double beta = (rho_new / rho) * (alpha / omega_bicg);

            #pragma acc parallel loop present(p, r, v)
            for (int i = 0; i < N; i++)
                p[i] = r[i] + beta * (p[i] - omega_bicg * v[i]);

            /* Precondition: z = M^-1 * p */
            #pragma acc parallel loop present(z, diag_inv, p)
            for (int i = 0; i < N; i++)
                z[i] = diag_inv[i] * p[i];

            /* v = A*z */
            apply_poisson_operator(ps, z, v);

            double r0v = 0.0;
            #pragma acc parallel loop reduction(+:r0v) present(r0, v)
            for (int i = 0; i < N; i++)
                r0v += r0[i] * v[i];

            if (fabs(r0v) < 1e-30) break;

            alpha = rho_new / r0v;

            #pragma acc parallel loop present(s, r, v)
            for (int i = 0; i < N; i++)
                s[i] = r[i] - alpha * v[i];

            /* Early exit if s is small enough */
            double norm_s = 0.0;
            #pragma acc parallel loop reduction(+:norm_s) present(s)
            for (int i = 0; i < N; i++)
                norm_s += s[i] * s[i];
            norm_s = sqrt(norm_s);

            if (norm_s / norm_r0 < tol) {
                #pragma acc parallel loop present(x, z)
                for (int i = 0; i < N; i++)
                    x[i] += alpha * z[i];
                break;
            }

            /* Precondition: y = M^-1 * s */
            #pragma acc parallel loop present(y, diag_inv, s)
            for (int i = 0; i < N; i++)
                y[i] = diag_inv[i] * s[i];

            /* t = A*y */
            apply_poisson_operator(ps, y, t);

            double ts = 0.0, tt = 0.0;
            #pragma acc parallel loop reduction(+:ts, tt) present(t, s)
            for (int i = 0; i < N; i++) {
                ts += t[i] * s[i];
                tt += t[i] * t[i];
            }

            if (fabs(tt) < 1e-30) break;
            omega_bicg = ts / tt;

            #pragma acc parallel loop present(x, z, y)
            for (int i = 0; i < N; i++)
                x[i] += alpha * z[i] + omega_bicg * y[i];

            #pragma acc parallel loop present(r, s, t)
            for (int i = 0; i < N; i++)
                r[i] = s[i] - omega_bicg * t[i];

            double norm_r = 0.0;
            #pragma acc parallel loop reduction(+:norm_r) present(r)
            for (int i = 0; i < N; i++)
                norm_r += r[i] * r[i];
            norm_r = sqrt(norm_r);

            #if BICGSTAB_DEBUG
            if (iter % 40 == 0)
                printf("BiCGStab iter %d  rel_res = %.6e  alpha=%.4e  omega=%.4e\n",
                        iter, norm_r/norm_b, alpha, omega_bicg);
            if (!isfinite(norm_r)) {
                printf("BiCGStab: NaN detected at iter %d, aborting.\n", iter);
                break;
            }
            #endif

            if (norm_r / norm_r0 < tol) break;
            if (fabs(omega_bicg) < 1e-30) break; // stagnation

            rho = rho_new;
        }

        if (!has_pressure_bc) {
            double mean = 0.0;
            #pragma acc parallel loop reduction(+:mean) present(x)
            for (int i = 0; i < N; i++)
                mean += x[i];
            
            mean /= (double)N;
            
            #pragma acc parallel loop present(x)
            for (int i = 0; i < N; i++)
                x[i] -= mean;
        }
    }

    free(r); free(r0); free(p); free(v);
    free(s); free(t); free(z); free(y);
    free(diag_inv);
}