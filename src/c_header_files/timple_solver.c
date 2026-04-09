// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules 2D
////////////////////////////////////////////////////////////////////////////////////

double time_implicit_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    double steady_state_error = 0.0;

    # pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field[0].u_old[i] = field[0].u[i]; 
        field[0].v_old[i] = field[0].v[i]; 
    }
    
    # pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        for (int iter = 0; iter<parameters.iter_timple; iter++){ 
            #pragma acc parallel loop
            for (int i = 0; i < myPointStruct->num_nodes; i++){
                field[0].u_new[i] = field[0].u[i];
                field[0].v_new[i] = field[0].v[i];
                field[0].pprime[i] = 0; 
            }
            
            calculate_intermediate_velocity_implicit_vectorised_2d(myPointStruct, field);
            calculate_mass_residual_implicit_vectorised_2d(myPointStruct, field);
            multigrid_Poisson_solver_vectorised(myPointStruct, field);
            update_velocity_implicit_vectorised_2d(myPointStruct, field);
            update_boundary_pressure_vectorised_2d(myPointStruct, field);
        }
    }

    # pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i=0; i<myPointStruct[0].num_nodes; i++){
        steady_state_error += pow(field[0].u[i]-field[0].u_old[i],2) + pow(field[0].v[i]-field[0].v_old[i],2);
    }
    return sqrt(steady_state_error/myPointStruct[0].num_nodes);
}

void calculate_intermediate_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    // Gradient calculation
    multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
    multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
    #pragma acc wait(1,2,3)
    
    for (int iter = 0; iter < parameters.iter_momentum; iter++) {
        #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
        for (int i = 0; i < num_nodes; i++) {
            if (!myPointStruct->boundary_tag[i]) {
                double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0;
                int base = i * num_cloud_points;

                double adv_x = field->u_new[i] * myPointStruct->Dx[base];
                double adv_y = field->v_new[i] * myPointStruct->Dy[base];
                double diff  = -parameters.mu * myPointStruct->lap[base];
                double unst  = 1.0 / (parameters.dt);
                double denom = parameters.rho * (adv_x + adv_y + unst) + diff;

                #pragma acc loop seq
                for (int j = 1; j < num_cloud_points; j++) {
                    int k = base + j;
                    int idx = myPointStruct->cloud_index[k];
                    t1 += myPointStruct->Dx[k] * field->u[idx];
                    t2 += myPointStruct->Dy[k] * field->u[idx];
                    t3 += myPointStruct->lap[k] * field->u[idx];
                    t4 += myPointStruct->Dx[k] * field->v[idx];
                    t5 += myPointStruct->Dy[k] * field->v[idx];
                    t6 += myPointStruct->lap[k] * field->v[idx];
                }

                double advection_u = field->u_new[i] * t1 + field->v_new[i] * t2;
                field->u[i] = (parameters.rho * (field->u_old[i]*unst - advection_u)  + parameters.mu*t3 - field->dpdx[i]) / denom;

                double advection_v = field->u_new[i] * t4 + field->v_new[i] * t5;
                field->v[i] = (parameters.rho * (field->v_old[i]*unst - advection_v)  + parameters.mu*t6 - field->dpdy[i]) / denom;
            }
        }

        #pragma acc parallel loop gang vector present(field, myPointStruct)
        for (int i = 0; i < num_nodes; i++) {
            if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
                if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET || myPointStruct->node_bc[i].type == BC_WALL) {
                    field->u[i] = myPointStruct->node_bc[i].u;
                    field->v[i] = myPointStruct->node_bc[i].v;
                }
            }
        }
    }
}

void calculate_mass_residual_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    # pragma acc data present(field, parameters, myPointStruct)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
    }
    # pragma acc wait(1,2)
    
    double sum = 0.0;
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) continue;

        if(!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho * (field->dudx[i] + field->dvdy[i]) / parameters.dt;
            sum += fabs(field->source[i]);
        }
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            field->source[i] = 0;
        }
        else {
            field->source[i] = parameters.rho * ((field->u[i] - field->u_old[i]) * myPointStruct->x_normal[i] + 
                (field->v[i] - field->v_old[i]) * myPointStruct->y_normal[i]) / parameters.dt;
        }
    }
}

void update_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    #pragma acc data present(field, parameters, myPointStruct)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->pprime, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->pprime, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
    }
    #pragma acc wait(1,2,3)
    #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (myPointStruct->corner_tag[i]) continue;

        if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
            double sumux = 0, sumuy = 0;
            double sumvx = 0, sumvy = 0;
            int base = i * num_cloud_points;
            
            #pragma acc loop seq
            for (int j = 1; j < num_cloud_points; j++) {
                int k = base + j;
                int idx = myPointStruct->cloud_index[k];
                sumux -= myPointStruct->Dx[k] * field->u[idx];
                sumuy -= myPointStruct->Dy[k] * field->u[idx];
                sumvx -= myPointStruct->Dx[k] * field->v[idx];
                sumvy -= myPointStruct->Dy[k] * field->v[idx];
            }
            
            double Ap = myPointStruct->Dx[base] * myPointStruct->x_normal[i] +
                        myPointStruct->Dy[base] * myPointStruct->y_normal[i];
            
            if (fabs(Ap) > 1e-12) {
                field->u[i] = (sumux * myPointStruct->x_normal[i] + sumuy * myPointStruct->y_normal[i]) / Ap;
                field->v[i] = (sumvx * myPointStruct->x_normal[i] + sumvy * myPointStruct->y_normal[i]) / Ap;
            }
        }
        else if(myPointStruct->node_bc[i].type == BC_WALL){
            field->u[i] = myPointStruct->node_bc[i].u;
            field->v[i] = myPointStruct->node_bc[i].v;
        }
        else if (!myPointStruct->boundary_tag[i]) {
            // Interior nodes: standard pressure correction
            field->u[i] -= (parameters.dt / parameters.rho) * field->dpdx[i];
            field->v[i] -= (parameters.dt / parameters.rho) * field->dpdy[i];
            field->p[i] += parameters.omega * field->pprime[i];
        }
    }
}

void update_boundary_pressure_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field) {
    double nu = parameters.mu / parameters.rho;
    int n = myPointStruct->num_cloud_points;
    int num_nodes = myPointStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for(int i = 0; i < num_nodes; i++) {
        if(myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
            if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
                field->p[i] = myPointStruct->node_bc[i].p;
                continue;
            }

            double t1 = 0, t2 = 0;
            int base = i * n;

            #pragma acc loop seq
            for(int j = 0; j < n; j++) {
                int idx = myPointStruct->cloud_index[base + j];
                
                t1 += nu * myPointStruct->lap[base + j] * field->u[idx];
                t1 -= field->u[i] * field->u[idx] * myPointStruct->Dx[base + j];
                t1 -= field->v[i] * field->u[idx] * myPointStruct->Dy[base + j];

                t2 += nu * myPointStruct->lap[base + j] * field->v[idx];
                t2 -= field->u[i] * field->v[idx] * myPointStruct->Dx[base + j];
                t2 -= field->v[i] * field->v[idx] * myPointStruct->Dy[base + j];
            }

            field->dpdn[i] = parameters.rho * (t1 * myPointStruct->x_normal[i] + 
                                               t2 * myPointStruct->y_normal[i]);

            double sumx = 0, sumy = 0;
            #pragma acc loop seq
            for (int j = 1; j < n; j++) {
                int idx = myPointStruct->cloud_index[base + j];
                sumx += myPointStruct->Dx[base + j] * field->p[idx];
                sumy += myPointStruct->Dy[base + j] * field->p[idx];
            }

            double Ap = myPointStruct->Dx[base] * myPointStruct->x_normal[i] +
                        myPointStruct->Dy[base] * myPointStruct->y_normal[i];

            if (fabs(Ap) > 1e-12) {
                field->p[i] = (field->dpdn[i] - sumx * myPointStruct->x_normal[i] 
                                             - sumy * myPointStruct->y_normal[i]) / Ap;
            }
        }
    }
     # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 5);
    }
    # pragma acc wait(4,5)
    double sum = 0.0;
    # pragma acc parallel loop gang vector default(present) reduction(+:sum)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        if (!myPointStruct->corner_tag[i])
            sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}


////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules 3D
////////////////////////////////////////////////////////////////////////////////////


double time_implicit_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    double steady_state_error = 0.0;

    # pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field[0].u_old[i] = field[0].u[i]; 
        field[0].v_old[i] = field[0].v[i]; 
        field[0].w_old[i] = field[0].w[i];
    }
    
    # pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        for (int iter = 0; iter<parameters.iter_timple; iter++){ 
            #pragma acc parallel loop
            for (int i = 0; i < myPointStruct->num_nodes; i++){
                field[0].u_new[i] = field[0].u[i];
                field[0].v_new[i] = field[0].v[i];
                field[0].w_new[i] = field[0].w[i];
                field[0].pprime[i] = 0; 
            }
            
            calculate_intermediate_velocity_implicit_vectorised(myPointStruct, field);
            calculate_mass_residual_implicit_vectorised(myPointStruct, field);
            multigrid_Poisson_solver_vectorised(myPointStruct, field);
            update_velocity_implicit_vectorised(myPointStruct, field);
            update_boundary_pressure_vectorised(myPointStruct, field);
        }
    }

    # pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i=0; i<myPointStruct[0].num_nodes; i++){
        steady_state_error += pow(field[0].u[i]-field[0].u_old[i],2) + pow(field[0].v[i]-field[0].v_old[i],2) + pow(field[0].w[i]-field[0].w_old[i],2);
    }
    return sqrt(steady_state_error/myPointStruct[0].num_nodes);
}

void calculate_intermediate_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    // Gradient calculation
    multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
    multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
    multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
    #pragma acc wait(1,2,3)
    
    for (int iter = 0; iter < parameters.iter_momentum; iter++) {
        #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
        for (int i = 0; i < num_nodes; i++) {
            if (!myPointStruct->boundary_tag[i]) {
                double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0, t9=0, t10=0, t11=0, t12=0;
                int base = i * num_cloud_points;

                double adv_x = field->u_new[i] * myPointStruct->Dx[base];
                double adv_y = field->v_new[i] * myPointStruct->Dy[base];
                double adv_z = field->w_new[i] * myPointStruct->Dz[base];
                double diff  = -parameters.mu * myPointStruct->lap[base];
                double unst  = 1.0 / (parameters.dt);
                double denom = parameters.rho * (adv_x + adv_y + adv_z + unst) + diff;

                #pragma acc loop seq
                for (int j = 1; j < num_cloud_points; j++) {
                    int k = base + j;
                    int idx = myPointStruct->cloud_index[k];
                    t1 += myPointStruct->Dx[k] * field->u[idx];
                    t2 += myPointStruct->Dy[k] * field->u[idx];
                    t3 += myPointStruct->Dz[k] * field->u[idx];
                    t4 += myPointStruct->lap[k] * field->u[idx];
                    t5 += myPointStruct->Dx[k] * field->v[idx];
                    t6 += myPointStruct->Dy[k] * field->v[idx];
                    t7 += myPointStruct->Dz[k] * field->v[idx];
                    t8 += myPointStruct->lap[k] * field->v[idx];
                    t9 += myPointStruct->Dx[k] * field->w[idx];
                    t10 += myPointStruct->Dy[k] * field->w[idx];
                    t11 += myPointStruct->Dz[k] * field->w[idx];
                    t12 += myPointStruct->lap[k] * field->w[idx];
                }

                double advection_u = field->u_new[i] * t1 + field->v_new[i] * t2 + field->w_new[i] * t3;
                field->u[i] = (parameters.rho * (field->u_old[i]*unst - advection_u) + parameters.mu*t4 - field->dpdx[i]) / denom;

                double advection_v = field->u_new[i] * t5 + field->v_new[i] * t6 + field->w_new[i] * t7;
                field->v[i] = (parameters.rho * (field->v_old[i]*unst - advection_v) + parameters.mu*t8 - field->dpdy[i]) / denom;

                double advection_w = field->u_new[i] * t9 + field->v_new[i] * t10 + field->w_new[i] * t11;
                field->w[i] = (parameters.rho * (field->w_old[i]*unst - advection_w) + parameters.mu*t12 - field->dpdz[i]) / denom;
            }
        }

        #pragma acc parallel loop gang vector present(field, myPointStruct)
        for (int i = 0; i < num_nodes; i++) {
            if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
                if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET || myPointStruct->node_bc[i].type == BC_WALL) {
                    field->u[i] = myPointStruct->node_bc[i].u;
                    field->v[i] = myPointStruct->node_bc[i].v;
                    field->w[i] = myPointStruct->node_bc[i].w;
                }
            }
        }
    }
}


void calculate_mass_residual_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    // Use velocity derivative fields instead of overwriting pressure gradients
    # pragma acc data present(field, parameters, myPointStruct)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w, field->dwdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
    }
    # pragma acc wait(1,2,3)
    
    double sum = 0.0;
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) continue;
        if(!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho * (field->dudx[i] + field->dvdy[i] + field->dwdz[i]) / parameters.dt;
            sum += fabs(field->source[i]);
        }
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            field->source[i] = 0.0;
        }
        else {
            field->source[i] = parameters.rho * (
                (field->u[i] - field->u_old[i]) * myPointStruct->x_normal[i] + 
                (field->v[i] - field->v_old[i]) * myPointStruct->y_normal[i] + 
                (field->w[i] - field->w_old[i]) * myPointStruct->z_normal[i]
            ) / parameters.dt;
        }
    }
    
}

void update_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    #pragma acc data present(field, parameters, myPointStruct)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->pprime, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->pprime, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->pprime, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
    }
    #pragma acc wait(1,2,3)

    #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (myPointStruct->corner_tag[i]) continue;

        if (myPointStruct->boundary_tag[i] && myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
            double sumux = 0, sumuy = 0, sumuz = 0;
            double sumvx = 0, sumvy = 0, sumvz = 0;
            double sumwx = 0, sumwy = 0, sumwz = 0;
            int base = i * num_cloud_points;
            
            #pragma acc loop seq
            for (int j = 1; j < num_cloud_points; j++) {
                int k = base + j;
                int idx = myPointStruct->cloud_index[k];
                sumux -= myPointStruct->Dx[k] * field->u[idx];
                sumuy -= myPointStruct->Dy[k] * field->u[idx];
                sumuz -= myPointStruct->Dz[k] * field->u[idx];
                sumvx -= myPointStruct->Dx[k] * field->v[idx];
                sumvy -= myPointStruct->Dy[k] * field->v[idx];
                sumvz -= myPointStruct->Dz[k] * field->v[idx];
                sumwx -= myPointStruct->Dx[k] * field->w[idx];
                sumwy -= myPointStruct->Dy[k] * field->w[idx];
                sumwz -= myPointStruct->Dz[k] * field->w[idx];
            }
            
            double Ap = myPointStruct->Dx[base] * myPointStruct->x_normal[i] +
                        myPointStruct->Dy[base] * myPointStruct->y_normal[i] +
                        myPointStruct->Dz[base] * myPointStruct->z_normal[i];
            
            if (fabs(Ap) > 1e-12) {
                field->u[i] = (sumux * myPointStruct->x_normal[i] + 
                               sumuy * myPointStruct->y_normal[i] + 
                               sumuz * myPointStruct->z_normal[i]) / Ap;
                field->v[i] = (sumvx * myPointStruct->x_normal[i] + 
                               sumvy * myPointStruct->y_normal[i] + 
                               sumvz * myPointStruct->z_normal[i]) / Ap;
                field->w[i] = (sumwx * myPointStruct->x_normal[i] + 
                               sumwy * myPointStruct->y_normal[i] + 
                               sumwz * myPointStruct->z_normal[i]) / Ap;
            }
        }
        else if (!myPointStruct->boundary_tag[i]) {
            // Interior nodes: standard pressure correction
            field->u[i] -= (parameters.dt / parameters.rho) * field->dpdx[i];
            field->v[i] -= (parameters.dt / parameters.rho) * field->dpdy[i];
            field->w[i] -= (parameters.dt / parameters.rho) * field->dpdz[i];
            field->p[i] += parameters.omega * field->pprime[i];
        }
        // Wall and inlet boundaries: velocity already set in intermediate step, don't touch
    }
    # pragma acc data present(myPointStruct, field, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 5);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points, 6);
    }
    # pragma acc wait(4,5,6)
    double sum = 0.0;
    # pragma acc parallel loop gang vector default(present) reduction(+:sum)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        if (!myPointStruct->corner_tag[i])
            sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]+field->dpdz[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}


void update_boundary_pressure_vectorised(PointStructure* myPointStruct, FieldVariables* field) {
    double nu = parameters.mu / parameters.rho;
    int n = myPointStruct->num_cloud_points;
    int num_nodes = myPointStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for(int i = 0; i < num_nodes; i++) {
        if(myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
            if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
                field->p[i] = myPointStruct->node_bc[i].p;
                continue;
            }

            double t1 = 0, t2 = 0, t3 = 0;
            int base = i * n;

            #pragma acc loop seq
            for(int j = 0; j < n; j++) {
                int idx = myPointStruct->cloud_index[base + j];
                
                t1 += nu * myPointStruct->lap[base + j] * field->u[idx];
                t1 -= field->u[i] * field->u[idx] * myPointStruct->Dx[base + j];
                t1 -= field->v[i] * field->u[idx] * myPointStruct->Dy[base + j];
                t1 -= field->w[i] * field->u[idx] * myPointStruct->Dz[base + j];

                t2 += nu * myPointStruct->lap[base + j] * field->v[idx];
                t2 -= field->u[i] * field->v[idx] * myPointStruct->Dx[base + j];
                t2 -= field->v[i] * field->v[idx] * myPointStruct->Dy[base + j];
                t2 -= field->w[i] * field->v[idx] * myPointStruct->Dz[base + j];

                t3 += nu * myPointStruct->lap[base + j] * field->w[idx];
                t3 -= field->u[i] * field->w[idx] * myPointStruct->Dx[base + j];
                t3 -= field->v[i] * field->w[idx] * myPointStruct->Dy[base + j];
                t3 -= field->w[i] * field->w[idx] * myPointStruct->Dz[base + j];
            }

            field->dpdn[i] = parameters.rho * (t1 * myPointStruct->x_normal[i] + 
                                               t2 * myPointStruct->y_normal[i] + 
                                               t3 * myPointStruct->z_normal[i]);

            double sumx = 0, sumy = 0, sumz = 0;
            #pragma acc loop seq
            for (int j = 1; j < n; j++) {
                int idx = myPointStruct->cloud_index[base + j];
                sumx += myPointStruct->Dx[base + j] * field->p[idx];
                sumy += myPointStruct->Dy[base + j] * field->p[idx];
                sumz += myPointStruct->Dz[base + j] * field->p[idx];
            }

            double Ap = myPointStruct->Dx[base] * myPointStruct->x_normal[i] +
                        myPointStruct->Dy[base] * myPointStruct->y_normal[i] +
                        myPointStruct->Dz[base] * myPointStruct->z_normal[i];

            if (fabs(Ap) > 1e-12) {
                field->p[i] = (field->dpdn[i] - sumx * myPointStruct->x_normal[i] 
                                             - sumy * myPointStruct->y_normal[i] 
                                             - sumz * myPointStruct->z_normal[i]) / Ap;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//  Poisson Solver
///////////////////////////////////////////////////////////////////////////////

void multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field) {
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++) {
        // Downward Path (Restriction)
        for (int ilev = 0; ilev < parameters.num_levels; ilev++) {
            if (parameters.poisson_solver_type == 1)    
                relaxation_vectorised_Jacobi(&myPointStruct[ilev], field[ilev].source, field[ilev].pprime, field[ilev].p_old);
            else if (parameters.poisson_solver_type == 3)
                relaxation_vectorised_BiCGStab(&myPointStruct[ilev], field[ilev].source, field[ilev].pprime, parameters.num_relax, parameters.poisson_solver_tolerance);
            else
                relaxation_vectorised_GaussSeidel(&myPointStruct[ilev], field[ilev].source, field[ilev].pprime);  
            calculate_residuals_vectorised(&myPointStruct[ilev], &field[ilev]);
            if (ilev < parameters.num_levels - 1) {
                restrict_residuals_vectorised(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
                #pragma acc parallel loop present(field[ilev+1])
                for(int i=0; i < myPointStruct[ilev+1].num_nodes; i++) 
                    field[ilev+1].pprime[i] = 0.0;
            }
        }
        
        // Upward Path (Prolongation)
        for (int ilev = parameters.num_levels - 1; ilev > 0; ilev--) {
            prolongate_corrections_vectorised(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (parameters.poisson_solver_type == 1)    
                relaxation_vectorised_Jacobi(&myPointStruct[ilev-1], field[ilev-1].source, field[ilev-1].pprime, field[ilev-1].p_old);
            else if (parameters.poisson_solver_type == 3)
                relaxation_vectorised_BiCGStab(&myPointStruct[ilev-1], field[ilev-1].source, field[ilev-1].pprime, parameters.num_relax, parameters.poisson_solver_tolerance);
            else
                relaxation_vectorised_GaussSeidel(&myPointStruct[ilev-1], field[ilev-1].source, field[ilev-1].pprime);  
        }
    }
}

void calculate_residuals_vectorised(PointStructure* mypointStruct, FieldVariables* field) {
    int n = mypointStruct->num_cloud_points;
    int num_nodes = mypointStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field, mypointStruct)
    for (int i = 0; i < num_nodes; i++) {
        field->res[i] = 0.0;
        if (!mypointStruct->boundary_tag[i] && !mypointStruct->corner_tag[i]) {
            double sum = 0;
            int base = i * n;
            #pragma acc loop seq 
            for (int j = 0; j < n; j++) {
                sum += mypointStruct->lap_Poison[base + j] * field->pprime[mypointStruct->cloud_index[base + j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}

void restrict_residuals_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                   FieldVariables* field_f, FieldVariables* field_c) {
    int n_f = mypointStruct_f->num_cloud_points;
    // int n_c = mypointStruct_c->num_cloud_points;
    int num_nodes_c = mypointStruct_c->num_nodes;

    #pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < num_nodes_c; i++) {
        double results = 0.0;
        if (mypointStruct_c->boundary_tag[i] == false) {
            int base_c = i * n_f; 
            int i_restr_node = mypointStruct_c->restriction_points[i];
            int base_f = i_restr_node * n_f;

            #pragma acc loop seq
            for (int j = 0; j < n_f; j++) {
                results += mypointStruct_c->restr_mat[base_c + j] * field_f->res[mypointStruct_f->cloud_index[base_f + j]];
            }
        }
        field_c->source[i] = results;
    }
}

void prolongate_corrections_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                       FieldVariables* field_f, FieldVariables* field_c) {
    int n_c = mypointStruct_c->num_cloud_points;
    int num_nodes_f = mypointStruct_f->num_nodes;
    int num_nodes_c = mypointStruct_c->num_nodes;

    #pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < num_nodes_f; i++) {
        // Only prolongate corrections for interior nodes
        if (mypointStruct_f->boundary_tag[i] == false) {
            int i_prol_node = mypointStruct_f->prolongation_points[i];
            int base_f = i * n_c;
            int base_c = i_prol_node * n_c;
            double results = 0.0;

            #pragma acc loop seq
            for (int j = 0; j < n_c; j++) {
                results += mypointStruct_f->prol_mat[base_f + j] * field_c->pprime[mypointStruct_c->cloud_index[base_c + j]];
            }
            field_f->pprime[i] += results;
        }
    }

    #pragma acc parallel loop gang vector present(field_c, mypointStruct_c)
    for (int i = 0; i < num_nodes_c; i++) {
        // Zero out only non-boundary nodes to maintain BC consistency
        if (mypointStruct_c->boundary_tag[i] == false) {
            field_c->pprime[i] = 0.0;
        }
    }
}


void relaxation_vectorised_Jacobi(PointStructure* mypointstruct, const double* source, double* pprime, double* p_old)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;
    
    #pragma acc parallel loop present(pprime[:N], p_old[:N])
    for (int i = 0; i < N; i++) {
        p_old[i] = pprime[i];
    }
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector present(pprime[:N], p_old[:N], mypointstruct->cloud_index[:N*n], mypointstruct->lap_Poison[:N*n], parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 1; j < n; j++) {
                int idx = i*n + j;
                sum += mypointstruct->lap_Poison[idx] * p_old[mypointstruct->cloud_index[idx]];
            }
            pprime[i] = parameters.omega * (source[i] - sum) / mypointstruct->lap_Poison[i*n] 
                            + (1.0 - parameters.omega) * p_old[i];
        }
        
        #pragma acc parallel loop present(pprime[:N], p_old[:N])
        for (int i = 0; i < N; i++) {
            p_old[i] = pprime[i];
        }
    }
}

void relaxation_vectorised_GaussSeidel(PointStructure* mypointstruct, const double* source, double* pprime) {
    int n = mypointstruct->num_cloud_points;
    int num_nodes = mypointstruct->num_nodes;

    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector present(pprime[:num_nodes], source[:num_nodes], parameters, mypointstruct)
        for (int i = 0; i < num_nodes; i++) {
            double sum = 0.0;
            int base = i * n;
            #pragma acc loop seq // Inner loops should be handled carefully for reduction
            for (int j = 1; j < n; j++) {
                sum += mypointstruct->lap_Poison[base + j] * pprime[mypointstruct->cloud_index[base + j]];
            }
            double diag = mypointstruct->lap_Poison[base];
            double new_val = (source[i] - sum) / diag;
            pprime[i] = parameters.omega * new_val + (1.0 - parameters.omega) * pprime[i];
        }
    }
}

void relaxation_vectorised_BiCGStab(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
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

    #pragma acc data copyin(b[0:N], ps->lap_Poison[0:N*n], ps->cloud_index[0:N*n]) \
                     copy(x[0:N]) \
                     create(diag_inv[0:N], r[0:N], r0[0:N], p[0:N], v[0:N], \
                            s[0:N], t[0:N], z[0:N], y[0:N])
    {
        /* Diagonal inverse (Jacobi preconditioner) */
        #pragma acc parallel loop
        for (int i = 0; i < N; i++) {
            double diag = ps->lap_Poison[i*n];
            if (fabs(diag) < 1e-14) diag = 1.0;
            diag_inv[i] = 1.0 / diag;
        }

        /* Initial residual r = b - A*x */
        #pragma acc parallel loop
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * x[col];
            }
            r[i]  = b[i] - sum;
            r0[i] = r[i];
            p[i]  = 0.0;
            v[i]  = 0.0;
        }

        double norm_r0 = 0.0;
        #pragma acc parallel loop reduction(+:norm_r0)
        for (int i = 0; i < N; i++) {
            norm_r0 += r[i] * r[i];
        }
        norm_r0 = sqrt(norm_r0);

        double rho = 1.0, alpha = 1.0, omega = 1.0;

        for (int iter = 0; iter < max_iter; iter++) {

            double rho_new = 0.0;
            #pragma acc parallel loop reduction(+:rho_new)
            for (int i = 0; i < N; i++) {
                rho_new += r0[i] * r[i];
            }

            if (fabs(rho_new) < 1e-30) break;

            double beta = (rho_new / rho) * (alpha / omega);

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                z[i] = diag_inv[i] * p[i];
            }

            /* v = A*z */
            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    int col = ps->cloud_index[i*n + j];
                    sum += ps->lap_Poison[i*n + j] * z[col];
                }
                v[i] = sum;
            }

            double r0v = 0.0;
            #pragma acc parallel loop reduction(+:r0v)
            for (int i = 0; i < N; i++) {
                r0v += r0[i] * v[i];
            }

            if (fabs(r0v) < 1e-30) break;

            alpha = rho_new / r0v;

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                s[i] = r[i] - alpha * v[i];
            }

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                y[i] = diag_inv[i] * s[i];
            }

            /* t = A*y */
            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    int col = ps->cloud_index[i*n + j];
                    sum += ps->lap_Poison[i*n + j] * y[col];
                }
                t[i] = sum;
            }

            double ts = 0.0, tt = 0.0;
            #pragma acc parallel loop reduction(+:ts,tt)
            for (int i = 0; i < N; i++) {
                ts += t[i] * s[i];
                tt += t[i] * t[i];
            }

            if (fabs(tt) < 1e-30) break;

            omega = ts / tt;

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                x[i] += alpha * z[i] + omega * y[i];
            }

            #pragma acc parallel loop
            for (int i = 0; i < N; i++) {
                r[i] = s[i] - omega * t[i];
            }

            double norm_r = 0.0;
            #pragma acc parallel loop reduction(+:norm_r)
            for (int i = 0; i < N; i++) {
                norm_r += r[i] * r[i];
            }
            norm_r = sqrt(norm_r);

            if (norm_r / norm_r0 < tol) break;

            rho = rho_new;
        }
    }

    free(r); free(r0); free(p); free(v);
    free(s); free(t); free(z); free(y);
    free(diag_inv);
}