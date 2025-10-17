// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef TIMPLE_SOLVERS_H
#define TIMPLE_SOLVERS_H

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

////////////////////////////////////////////////////////////////////////////////////
// Function Declarations
//////////////////////////////////////////////////////////////////////////////////
double time_implicit_solver(PointStructure* myPointStruct, FieldVariables* field);
void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, FieldVariables* field);
void calculate_mass_residual_implicit(PointStructure* myPointStruct, FieldVariables* field);
void multigrid_Poisson_solver(PointStructure* myPointStruct, FieldVariables* field);
void update_velocity_implicit(PointStructure* myPointStruct, FieldVariables* field);
void relaxation(PointStructure* mypointstruct, FieldVariables* field);
void calculate_residuals(PointStructure* myPointStruct, FieldVariables* field);
void restrict_residuals(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void prolongate_corrections(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void update_boundary_pressure(PointStructure* mypointstruct, FieldVariables* field);
void update_boundary_pprime(PointStructure* mypointstruct, FieldVariables* field);

////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules
////////////////////////////////////////////////////////////////////////////////////
double time_implicit_solver(PointStructure* myPointStruct, FieldVariables* field){
    double steady_state_error = 0.0;

    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field[0].u_old[i] = field[0].u_new[i]; field[0].u[i] = field[0].u_new[i];
        field[0].v_old[i] = field[0].v_new[i]; field[0].v[i] = field[0].v_new[i];
        field[0].w_old[i] = field[0].w_new[i]; field[0].w[i] = field[0].w_new[i];
        field[0].pprime[i] = 0;
    }
    for (int iter = 0; iter<parameters.iter_timple; iter++){ 
           for (int i = 0; i < myPointStruct->num_nodes; i++){
            field[0].u[i] = field[0].u_new[i];
            field[0].v[i] = field[0].v_new[i];
            field[0].w[i] = field[0].w_new[i];
            field[0].pprime[i] = 0; 
           }
        calculate_intermediate_velocity_implicit(myPointStruct, field);
        calculate_mass_residual_implicit(myPointStruct, field);
        multigrid_Poisson_solver(myPointStruct, field);
        update_velocity_implicit(myPointStruct, field);
        update_boundary_pressure(myPointStruct, field);
    }

    for (int i=0; i<myPointStruct[0].num_nodes; i++){
        steady_state_error += pow(field->u_new[i]-field->u_old[i],2) + pow(field->v_new[i]-field->v_old[i],2);
        if (parameters.dimension == 3)
        steady_state_error += pow(field->w_new[i]-field->w_old[i],2);
    }
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}

void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    double denom = 0;
    double temp1, temp2, temp3, temp4;
    double advection_x, advection_y, advection_z, advection, diffusion, unsteady;

    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->p_old, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->p_old, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->p_old, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    }   
    
    # pragma acc loop
    for (int iter = 0; iter<parameters.iter_momentum; iter++){
        # pragma acc loop
        for (int i = 0; i < num_nodes; i++){
            if(!myPointStruct->boundary_tag[i]){
        // u-velocity
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;			
                advection_x = field->u[i] * myPointStruct->Dx[i][0];                      
                advection_y = field->v[i] * myPointStruct->Dy[i][0]; 
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[i][0];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + unsteady) + diffusion; 
                # pragma acc loop
                for (int j = 1; j < num_cloud_points; j++){
                    temp1  += myPointStruct->Dx[i][j] * field->u_new[myPointStruct->cloud_index[i][j]];
                    temp2  += myPointStruct->Dy[i][j] * field->u_new[myPointStruct->cloud_index[i][j]];
                    temp4  += myPointStruct->lap[i][j]* field->u_new[myPointStruct->cloud_index[i][j]];
                }
                if (parameters.dimension == 3){
                    denom  += parameters.rho * field->w[i] * myPointStruct->Dz[i][0];
                # pragma acc loop
                for (int j = 1; j < num_cloud_points; j++)
                    temp3  += myPointStruct->Dz[i][j]* field->u_new[myPointStruct->cloud_index[i][j]];
                }
                advection = field->u[i] * temp1 + field->v[i] * temp2 + field->w[i] * temp3;
                unsteady  = field->u_old[i]/(parameters.dt*parameters.facdt);
                field->u_new[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) *   parameters.mu* myPointStruct->lap[i][0] * field->u_new[i] + 
                               parameters.mu*temp4 - field->dpdx[i]) / denom;
        // v-velocity
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;			
	            advection_x = field->u[i] * myPointStruct->Dx[i][0];                      
	            advection_y = field->v[i] * myPointStruct->Dy[i][0]; 
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[i][0];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + unsteady) + diffusion; 
                # pragma acc loop
                for (int j = 1; j < num_cloud_points; j++){
                    temp1 += myPointStruct->Dx[i][j] * field->v_new[myPointStruct->cloud_index[i][j]];
                    temp2 += myPointStruct->Dy[i][j] * field->v_new[myPointStruct->cloud_index[i][j]];
                    temp4 += myPointStruct->lap[i][j]* field->v_new[myPointStruct->cloud_index[i][j]];
                }
                if (parameters.dimension == 3){
                    denom += parameters.rho * field->w[i] * myPointStruct->Dz[i][0];
                    # pragma acc loop
                    for (int j = 1; j < num_cloud_points; j++)
                        temp3 += myPointStruct->Dz[i][j]* field->v_new[myPointStruct->cloud_index[i][j]];
                }
                advection = field->u[i] * temp1 + field->v[i] * temp2 + field->w[i] * temp3;
                unsteady  = field->v_old[i]/(parameters.dt*parameters.facdt);
                field->v_new[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i][0] * field->v_new[i] +
                                    parameters.mu*temp4 - field->dpdy[i]) / denom;

                if (parameters.dimension==3){
            // w-velocity
                    temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;			
                    advection_x = field->u[i] * myPointStruct->Dx[i][0];                      
                    advection_y = field->v[i] * myPointStruct->Dy[i][0];
                    advection_z = field->w[i] * myPointStruct->Dz[i][0];  
                    diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[i][0];
                    unsteady    = 1./(parameters.dt*parameters.facdt);			
                    denom       = parameters.rho* (advection_x + advection_y + advection_z + unsteady) + diffusion; 
                    # pragma acc loop
                    for (int j = 1; j < num_cloud_points; j++){
                        temp1 += myPointStruct->Dx[i][j] * field->w_new[myPointStruct->cloud_index[i][j]];
                        temp2 += myPointStruct->Dy[i][j] * field->w_new[myPointStruct->cloud_index[i][j]];
                        temp3 += myPointStruct->Dz[i][j] * field->w_new[myPointStruct->cloud_index[i][j]];
                        temp4 += myPointStruct->lap[i][j]* field->w_new[myPointStruct->cloud_index[i][j]];
                    }
                    advection = field->u[i] * temp1 + field->v[i] * temp2 + field->w[i] * temp3;
                    unsteady  = field->w_old[i]/(parameters.dt*parameters.facdt);
                    field->w_new[i] = (parameters.rho * (unsteady - advection)- (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i][0] * field->w_new[i] + 
                                    parameters.mu*temp4 - field->dpdz[i]) / denom;
                }                   
            }
        }
    }
}

void calculate_mass_residual_implicit(PointStructure* myPointStruct, FieldVariables* field){

    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    double sum = 0.0;
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->u_new, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->v_new, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->w_new, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);

    # pragma acc parallel loop present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
            sum += (field->source[i]*parameters.dt);
      	}
      	else{
      	    field->source[i] = 0;
    	}
    }
    printf("Mass residual: %e\n", fabs(sum)/num_nodes);         
    # pragma acc update host(field[0].source[:num_nodes])
}
    
void multigrid_Poisson_solver(PointStructure* myPointStruct, FieldVariables* field){
    # pragma acc loop
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        # pragma acc loop
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            relaxation(&myPointStruct[ilev], &field[ilev]);  
            calculate_residuals(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                restrict_residuals(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }
        }
        # pragma acc loop
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            prolongate_corrections(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1)
                relaxation(&myPointStruct[ilev-1], &field[ilev-1]);
        }
    }        
}

void update_boundary_pprime(PointStructure* mypointstruct, FieldVariables* field){
    double sumx, sumy, sumz, Ap;
    if (parameters.neumann_flag_boundary){
        # pragma acc parallel loop present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_boundary_nodes; i++){
            sumx = 0; sumy = 0; sumz = 0; Ap = 0;
            # pragma acc loop reduction(+:sumx, sumy, sumz, Ap)
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                sumx += mypointstruct->Dx[i][j]*field->pprime[mypointstruct->cloud_index[i][j]];
                sumy += mypointstruct->Dy[i][j]*field->pprime[mypointstruct->cloud_index[i][j]];
                if (parameters.dimension == 3){
                    sumz += mypointstruct->Dz[i][j]*field->pprime[mypointstruct->cloud_index[i][j]];
                }
            }
            Ap += mypointstruct->Dx[i][0]*mypointstruct->x_normal[i];
            Ap += mypointstruct->Dy[i][0]*mypointstruct->y_normal[i];
            if (parameters.dimension == 3){
                Ap += mypointstruct->Dz[i][0]*mypointstruct->z_normal[i];
            }
            field->pprime[i] = (-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i] -sumz*mypointstruct->z_normal[i])/Ap;
        }
    }
}

void relaxation(PointStructure* mypointstruct, FieldVariables* field){
     //double* zeros=create_vector(mypointstruct->num_nodes);
    # pragma acc loop
    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            double sum = 0.0;
            # pragma acc loop reduction(+:sum)
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                sum += mypointstruct->lap_Poison[i][j]*field->pprime[mypointstruct->cloud_index[i][j]];
            }
            field->pprime[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[i][0]) + (1-parameters.omega)*field->pprime[i]; 
        }
        double pref = field->pprime[0];
        # pragma acc parallel loop present(field, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->pprime[i] = field->pprime[i] - pref;
        }
        
        # pragma acc parallel loop present(field, mypointstruct)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            if(!mypointstruct->boundary_tag[i]){
                field->res[i]=field->source[i];
                # pragma acc loop
                for (int j = 0; j < mypointstruct->num_cloud_points; j++){
                    field->res[i] = field->res[i] - mypointstruct->lap_Poison[i][j]*field->pprime[mypointstruct->cloud_index[i][j]];
                }
            }
            else{
                field->res[i]=0;
            }
        }
           //printf("Iteration: %d, Pressure Residual: %e\n", iter, l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));    
    }
           // printf(" Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));    
           //printf("Pressure corrections: %e\n", l2_norm_gen(mypointstruct, field->pprime, zeros, mypointstruct->num_nodes));
    //free(zeros);
}

void calculate_residuals(PointStructure* mypointStruct, FieldVariables* field){
    double* zeros=create_vector(mypointStruct->num_nodes);
    # pragma acc parallel loop present(field, mypointStruct)
    for (int i = 0; i < mypointStruct->num_nodes; i++){
           field->res[i] = 0.0;
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            # pragma acc loop reduction(+:sum)
            for (int j = 0; j < mypointStruct->num_cloud_points; j++){
                sum += mypointStruct->lap_Poison[i][j]*field->pprime[mypointStruct->cloud_index[i][j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
     //printf("Pressure Residual: %e\n", l2_norm_gen(mypointStruct, field->res, zeros, mypointStruct->num_nodes));
     free(zeros);
}

void restrict_residuals(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        double results = 0.0;
        if (mypointStruct_c->boundary_tag[i]==false){
            int i_restr = mypointStruct_c->restriction_points[i];
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_f->num_cloud_points; j++){
                results += mypointStruct_c->restr_mat[i][j] * field_f->res[mypointStruct_f->cloud_index[i_restr][j]];
            }
        }
                field_c->source[i] = results;
    }
}

void prolongate_corrections(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol = mypointStruct_f->prolongation_points[i];
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_c->num_cloud_points; j++){
                results += mypointStruct_f->prol_mat[i][j] * field_c->pprime[mypointStruct_c->cloud_index[i_prol][j]];
            }
            field_f->pprime[i] = field_f->pprime[i] + results;
        }
    }
    # pragma acc parallel loop present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
       // if (mypointStruct_c->boundary_tag[i]==false){
            field_c->pprime[i] = 0.0;
       // }
    }
}

void update_velocity_implicit(PointStructure* myPointStruct, FieldVariables* field){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    multiply_sparse_matrix_vector_gpu(myPointStruct->Dx, field->pprime, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dy, field->pprime, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
    multiply_sparse_matrix_vector_gpu(myPointStruct->Dz, field->pprime, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    }

    // Update Interior nodes
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->u_new[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v_new[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            if (parameters.dimension == 3){
                field->w_new[i] = field->w_new[i] - parameters.dt * field->dpdz[i]/parameters.rho;
            }
	        field->p_old[i] = field->p_old[i] + 1.0*field->pprime[i];
        }
    }
    double pref = field->p_old[0];
    for (int i = 0; i < num_nodes; i++){
        field->p_old[i] = field->p_old[i] - pref;
    }
}

void update_boundary_pressure(PointStructure* myPointStruct, FieldVariables* field){
    double temp1, temp2, temp3;
    double sumx, sumy, sumz, Ap;
    double nu = parameters.mu / parameters.rho;

    for(int i = 0; i < myPointStruct->num_nodes; i++){
    if(myPointStruct->boundary_tag[i]){
        temp1 = 0; temp2 = 0; temp3 = 0;
        for(int j = 0; j<myPointStruct->num_cloud_points; j++){
            temp1 = nu * myPointStruct->lap[i][j] * field->u_new[myPointStruct->cloud_index[i][j]];
            temp1 -= field->u[i]* field->u_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
            temp1 -= field->v[i]* field->u_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];
        //  temp1 -= (field->u[i] - field->u_old[i])/parameters.dt;
            if (parameters.dimension == 3){
                temp1 -= field->w[i]* field->u_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
            temp2 = nu* myPointStruct->lap[i][j] * field->v_new[myPointStruct->cloud_index[i][j]];
            temp2 -= field->u[i]* field->v_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
            temp2 -= field->v[i]* field->v_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];     
            //   temp2 -= (field->v[i] - field->v_old[i])/parameters.dt;
            if (parameters.dimension == 3){
                temp2 -= field->w[i]* field->v_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
            
            if (parameters.dimension == 3){
                temp3 = nu* myPointStruct->lap[i][j] * field->w_new[myPointStruct->cloud_index[i][j]];
                temp3 -= field->u[i]* field->w_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
                temp3 -= field->v[i]* field->w_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];
                temp3 -= field->w[i]* field->w_new[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
        }
        field->dpdn[i] = field->rho * (temp1*myPointStruct->x_normal[i] + temp2*myPointStruct->y_normal[i] + temp3*myPointStruct->z_normal[i]);
		
        sumx = 0; sumy = 0; sumz = 0; Ap = 0.0;
        # pragma acc loop reduction(+:sumx, sumy, sumz, Ap)

        for (int j = 1; j < myPointStruct->num_cloud_points; j++){
            sumx += myPointStruct->Dx[i][j]*field->p_old[myPointStruct->cloud_index[i][j]];
            sumy += myPointStruct->Dy[i][j]*field->p_old[myPointStruct->cloud_index[i][j]];
            if (parameters.dimension == 3){
                sumz += myPointStruct->Dz[i][j]*field->p_old[myPointStruct->cloud_index[i][j]];
            }
        }
            
        Ap += myPointStruct->Dx[i][0]*myPointStruct->x_normal[i];
        Ap += myPointStruct->Dy[i][0]*myPointStruct->y_normal[i];
        if (parameters.dimension == 3){
            Ap += myPointStruct->Dz[i][0]*myPointStruct->z_normal[i];
        }
        field->p_old[i] =(field->dpdn[i] - sumx*myPointStruct->x_normal[i] -sumy*myPointStruct->y_normal[i] - sumz*myPointStruct->z_normal[i])/Ap;
    	}
    }
}
	
#endif

