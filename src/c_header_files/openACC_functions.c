// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "structures.h"
#include "functions.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Function definitions
void copyin_parameters_to_gpu(){
    # pragma acc enter data copyin(parameters)
}

void copyin_pointstructure_to_gpu(PointStructure *myPointStruct)
{
    int L = parameters.num_levels;

    /* Create struct shells */
    #pragma acc enter data copyin(myPointStruct[:L])

    for (int l = 0; l < L; l++) {

        int N  = myPointStruct[l].num_nodes;
        int Nc = myPointStruct[l].num_cloud_points;

        /* ---------------- coordinates ---------------- */
        #pragma acc enter data copyin( \
            myPointStruct[l].x[:N], \
            myPointStruct[l].y[:N], \
            myPointStruct[l].z[:N], \
            myPointStruct[l].x_normal[:N], \
            myPointStruct[l].y_normal[:N], \
            myPointStruct[l].z_normal[:N], \
            myPointStruct[l].corner_tag[:N], \
            myPointStruct[l].boundary_tag[:N], \
            myPointStruct[l].point_index[:N], \
            myPointStruct[l].node_bc[:N]  )

        #pragma acc enter data attach( \
            myPointStruct[l].x, \
            myPointStruct[l].y, \
            myPointStruct[l].z, \
            myPointStruct[l].x_normal, \
            myPointStruct[l].y_normal, \
            myPointStruct[l].z_normal, \
            myPointStruct[l].corner_tag, \
            myPointStruct[l].boundary_tag, \
            myPointStruct[l].point_index, \
            myPointStruct[l].node_bc )

        /* ---------------- cloud connectivity ---------------- */
        #pragma acc enter data copyin( \
            myPointStruct[l].cloud_index[:N * Nc], \
            myPointStruct[l].Dx[:N * Nc], \
            myPointStruct[l].Dy[:N * Nc], \
            myPointStruct[l].lap[:N * Nc], \
            myPointStruct[l].lap_Poison[:N * Nc] )

        #pragma acc enter data attach( \
            myPointStruct[l].cloud_index, \
            myPointStruct[l].Dx, \
            myPointStruct[l].Dy, \
            myPointStruct[l].lap, \
            myPointStruct[l].lap_Poison )

        if (parameters.dimension == 3) {
            #pragma acc enter data copyin(myPointStruct[l].Dz[:N * Nc])
            #pragma acc enter data attach(myPointStruct[l].Dz)
        }
        
        /* ---------------- Poisson / multigrid operators ---------------- */
        if (myPointStruct[l].prol_mat != NULL) {
            #pragma acc enter data copyin(myPointStruct[l].prol_mat[:N * Nc])
            #pragma acc enter data attach(myPointStruct[l].prol_mat)
        }

        if (myPointStruct[l].restr_mat != NULL) {
            #pragma acc enter data copyin(myPointStruct[l].restr_mat[:N * Nc])
            #pragma acc enter data attach(myPointStruct[l].restr_mat)
        }
    }
}


void copyin_field_to_gpu(FieldVariables *field,
                         PointStructure *myPointStruct)
{
    int L = parameters.num_levels;

    #pragma acc enter data copyin(field[:L])

    for (int l = 0; l < L; l++) {

        int N = myPointStruct[l].num_nodes;

        /* ---------------- velocity ---------------- */
        #pragma acc enter data copyin( \
            field[l].u[:N], \
            field[l].v[:N], \
            field[l].u_old[:N], \
            field[l].v_old[:N], \
            field[l].u_new[:N], \
            field[l].v_new[:N], \
            field[l].dudx[:N], \
            field[l].dudy[:N], \
            field[l].dvdx[:N], \
            field[l].dvdy[:N], \
            field[l].lapu[:N], \
            field[l].lapv[:N], \
            field[l].p[:N], \
            field[l].p_old[:N], \
            field[l].pprime[:N], \
            field[l].dpdn[:N], \
            field[l].dpdx[:N], \
            field[l].dpdy[:N]  )

        #pragma acc enter data attach( \
            field[l].u, \
            field[l].v, \
            field[l].u_old, \
            field[l].v_old, \
            field[l].u_new, \
            field[l].v_new, \
            field[l].dudx, \
            field[l].dudy, \
            field[l].dvdx, \
            field[l].dvdy, \
            field[l].lapu, \
            field[l].lapv, \
            field[l].p, \
            field[l].p_old, \
            field[l].pprime, \
            field[l].dpdn, \
            field[l].dpdx, \
            field[l].dpdy )

        if (parameters.dimension == 3) {
            #pragma acc enter data copyin( \
                field[l].w[:N], \
                field[l].w_old[:N], \
                field[l].w_new[:N], \
                field[l].dwdx[:N], \
                field[l].dwdy[:N], \
                field[l].dwdz[:N], \
                field[l].dvdz[:N], \
                field[l].dudz[:N], \
                field[l].lapw[:N], \
                field[l].dpdz[:N] )
            #pragma acc enter data attach( \
                field[l].w, \
                field[l].w_old, \
                field[l].w_new, \
                field[l].dwdx, \
                field[l].dwdy, \
                field[l].dwdz, \
                field[l].dvdz, \
                field[l].dudz, \
                field[l].lapw, \
                field[l].dpdz )
        }

        /* ---------------- res, source ---------------- */
        #pragma acc enter data copyin( \
            field[l].res[:N], \
            field[l].source[:N] )

        #pragma acc enter data attach( \
            field[l].res, \
            field[l].source )
        /* ---------------- compressible variables ---------------- */
        if (parameters.compressible_flow) {

            #pragma acc enter data copyin( \
                field[l].rho[:N], field[l].rho_old[:N], field[l].rho_new[:N], \
                field[l].T_new[:N], field[l].T[:N], field[l].T_old[:N], \
                field[l].e[:N], field[l].e_old[:N], \
                field[l].drhodx[:N], field[l].drhody[:N], field[l].drhodz[:N], \
                field[l].dTdx[:N], field[l].dTdy[:N], field[l].dTdz[:N], \
                field[l].dedx[:N], field[l].dedy[:N], field[l].dedz[:N], \
                field[l].mu[:N], field[l].kappa[:N], \
                field[l].tau_xx[:N], field[l].tau_yy[:N], field[l].tau_zz[:N], \
                field[l].tau_xy[:N], field[l].tau_xz[:N], field[l].tau_yz[:N], \
                field[l].div_tau_x[:N], field[l].div_tau_y[:N], field[l].div_tau_z[:N], \
                field[l].Q_visc[:N], field[l].Q_source[:N] )

            #pragma acc enter data attach( \
                field[l].rho, field[l].rho_old, field[l].rho_new, \
                field[l].T_new, field[l].T, field[l].T_old, \
                field[l].e, field[l].e_old, \
                field[l].drhodx, field[l].drhody, field[l].drhodz, \
                field[l].dTdx, field[l].dTdy, field[l].dTdz, \
                field[l].dedx, field[l].dedy, field[l].dedz, \
                field[l].mu, field[l].kappa, \
                field[l].tau_xx, field[l].tau_yy, field[l].tau_zz, \
                field[l].tau_xy, field[l].tau_xz, field[l].tau_yz, \
                field[l].div_tau_x, field[l].div_tau_y, field[l].div_tau_z, \
                field[l].Q_visc, field[l].Q_source )
        }
    }
}

void copypout_pointstructure_from_gpu(PointStructure* myPointStruct){
    # pragma acc exit data copyout(myPointStruct[:parameters.num_levels])
}

void copypout_field_from_gpu(FieldVariables* field, PointStructure* myPointStruct){
    # pragma acc exit data copyout(field[:parameters.num_levels])
}