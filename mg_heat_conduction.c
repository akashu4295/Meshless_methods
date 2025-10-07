// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// compile command: gcc mg_heat_conduction.c -lm
// Run command: ./a.out

////////////////////////////////////////////////////////////////
// Multigrid solver for Heat conduction equation
// using Radial Basis Function Finite Difference (RBF-FD) method
////////////////////////////////////////////////////////////////

////////////// Header files

#include "header_files/heat_conduction.h"
#include "header_files/openACC_functions.h"
#include "init/init_heat_conduction.c"

////////////// Main Program

int main()
{
    clock_t clock_start = clock(), clock_program_begin = clock();
    PointStructure* myPointStruct;
    FieldVariables* field;
    double steady_state_error = 1e10;
    FILE *file1;
    FILE *file2;

///// Read Input File Name and Parameters
    read_flow_parameters("flow_parameters.csv");     // Read the parameters from the file
    read_grid_filenames(&myPointStruct, "grid_filenames.csv", &parameters.num_levels);     // Read the parameters from the file

///// Read Mesh data on all levels and compute derivative matrices
    read_complete_mesh_data(myPointStruct, parameters.num_levels);
    double time_read_mesh = (double)(clock()-clock_start)/CLOCKS_PER_SEC;

///// Create and/or test derivative matrices
    clock_start = clock();    // Start the clock
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_derivative_matrices(&myPointStruct[ii]);
    double time_create_derivatives = (double)(clock()-clock_start)/CLOCKS_PER_SEC;
    clock_start = clock();    // Start the clock
    if(parameters.test>0)
        test_derivatives (myPointStruct, parameters.num_levels, parameters.dimension);
    double time_test_derivatives = (double)(clock()-clock_start)/CLOCKS_PER_SEC;

///// Setting up the boudary condition (Lid driven cavity)
    clock_start = clock();    // Start the clock    
    AllocateMemoryFieldVariables(&field, myPointStruct, parameters.num_levels);

    int k = 1;  // wave number
    initial_conditions(myPointStruct, field, parameters.num_levels, k);
    boundary_conditions(myPointStruct, field, parameters.num_levels, k);

    double time_presolve = (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC;
    double* temp = create_vector(myPointStruct[0].num_nodes);

///// Time stepping 
    file2 = fopen("Convergence.csv", "w"); // Write data to a file

    for (int it =0; it < parameters.num_time_steps; it++) 
    {   
        for (int i = 0; i < myPointStruct[0].num_nodes; i++){
            temp[i] = field[0].p[i];
        }
        for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
            set_boundary_rhs(&myPointStruct[0], field, parameters.num_levels, k);
            for (int ilev = 0; ilev < parameters.num_levels; ilev++){
                relaxation(&myPointStruct[ilev], &field[ilev]);
                calculate_residuals(&myPointStruct[ilev], &field[ilev]);
                if (ilev != parameters.num_levels-1){
                    restrict_residuals(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
                    }
                }
            for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
                prolongate_corrections(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
                if (ilev != 1) 
                    relaxation(&myPointStruct[ilev-1], &field[ilev-1]);
            } 
        }
        fprintf(file2, "%d, %lf\n", it, steady_state_error);  
        steady_state_error = 0.0;  
        for (int i = 0; i < myPointStruct[0].num_nodes; i++){
            steady_state_error += (temp[i] - field[0].p[i])*(temp[i] - field[0].p[i]);
        }
        printf("Time step: %d, Steady state error: %e\n", it, sqrt(steady_state_error/myPointStruct[0].num_nodes));

        if (steady_state_error < parameters.steady_state_tolerance){
            break;
        }
    }

    file1 = fopen("Solution.csv", "w");
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        fprintf(file1, "%lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], field[0].p[i]);
        fflush(file1);
    }

    printf("Time_step, dt : %lf\n",parameters.dt);
    printf("Courant number: %lf\n",parameters.courant_number);
    printf("Average point spacing: %lf\n",myPointStruct[0].d_avg);
    printf("Polynomial degree: %d\n",myPointStruct[0].poly_degree);
    printf("Time for reading mesh data: %lf\n", time_read_mesh);
    printf("Time for creating derivative matrices: %lf\n", time_create_derivatives);
    printf("Time for testing derivatives: %lf\n", time_test_derivatives);
    printf("Time for presolve: %lf\n", time_presolve);
    // printf("Time for copying data to GPU: %lf\n", time_copy_gpu);
    printf("Time for execution (total): %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);
    
    
    // Clear read data
    free_PointStructure(myPointStruct, parameters.num_levels);
    free_field(field, parameters.num_levels);

    fclose(file1);
    fclose(file2); 
    return 0;
} 
