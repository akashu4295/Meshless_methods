// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef INIT_C
#define INIT_C

#include "../header_files/structures.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);
void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);

void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    double x,y,z,r;
    double r_o = 2.0, omega_i = -1.0;
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            myfieldvariables[ii].u[i] = 0;
            myfieldvariables[ii].v[i] = 0;
            myfieldvariables[ii].w[i] = 0;
            myfieldvariables[ii].u_new[i] = 0;
            myfieldvariables[ii].v_new[i] = 0;
            myfieldvariables[ii].w_new[i] = 0;
            myfieldvariables[ii].p[i] = 0;
            myfieldvariables[ii].p_old[i] = 0;
            myfieldvariables[ii].pprime[i] = 0;
            
        }
    }
}

void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    double x, y, r;
    double r_o = 2, omega_i = 1;
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            x = myPointStruct[ii].x[i]; y = myPointStruct[ii].y[i];
            r = sqrt(x*x + y*y);   
            if (fabs(r-r_o)<1e-9){
                myfieldvariables[ii].u[i] = -y/r;
                myfieldvariables[ii].v[i] = x/r;
                myfieldvariables[ii].p[i] = 0.0;
                myfieldvariables[ii].u_new[i] = -y/r;
                myfieldvariables[ii].v_new[i] = x/r;
                myfieldvariables[ii].pprime[i] = 0.0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
        }
    }
}


#endif
