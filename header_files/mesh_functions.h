// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H

#include "structures.h"
#include "kdtree_functions.h"
#include "general_functions.h"
#include "write_functions.h"
#include "rbf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double dtemp; int itemp; char temp[50]; // temporary variables used multiple times in the code

///////////////////////////////////////////////////////////////////////////////
// Function definitions
///////////////////////////////////////////////////////////////////////////////

void read_PointStructure(PointStructure* myPointStruct);
void read_flow_parameters(char *filename);
void calculate_parameters(PointStructure* myPointStruct);
void correct_normal_directions(PointStructure* myPointStruct);
void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels);
void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels);   
void create_restriction_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);
void create_prolongation_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);
void rcm_reordering(PointStructure* myPointStruct);

///////////////////////////////////////////////////////////////////////////////
// Function definitions
///////////////////////////////////////////////////////////////////////////////

void read_flow_parameters(char *filename) {
    FILE *file;
    char ctemp[100];
    int temp, temp1;

    file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.dimension);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.dimension);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.poly_degree);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.poly_degree);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.phs_degree);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.phs_degree);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.cloud_size_multiplier);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.cloud_size_multiplier);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.test);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.test);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.courant_number);
    printf("PARAMETERS: %s = %lf\n", ctemp, parameters.courant_number);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.steady_state_tolerance);
    printf("PARAMETERS: %s = %e\n", ctemp, parameters.steady_state_tolerance);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.poisson_solver_tolerance);
    printf("PARAMETERS: %s = %e\n", ctemp, parameters.poisson_solver_tolerance);
    temp = fscanf(file, "%[^,],%f\n", ctemp, &parameters.omega);
    printf("PARAMETERS: %s = %f\n", ctemp, parameters.omega);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.dt);
    printf("PARAMETERS: %s = %lf\n", ctemp, parameters.dt);
    temp = fscanf(file, "%[^,],%d\n", ctemp, &parameters.num_time_steps);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.num_time_steps);
    temp = fscanf(file, "%[^,],%d\n", ctemp, &parameters.write_interval);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.write_interval);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.Re);
    printf("PARAMETERS: %s = %lf\n", ctemp, parameters.Re);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.iter_momentum);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.iter_momentum);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.iter_timple);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.iter_timple);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.num_vcycles);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.num_vcycles);
    temp = fscanf(file, "%[^,],%hd\n", ctemp, &parameters.num_relax);
    printf("PARAMETERS: %s = %hd\n", ctemp, parameters.num_relax);
    temp = fscanf(file, "%[^,],%d\n", ctemp, &temp1);
    if (temp1 == 1)
        parameters.neumann_flag_boundary = true;
    else
        parameters.neumann_flag_boundary = false;
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.neumann_flag_boundary);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.facRe);
    printf("PARAMETERS: %s = %lf\n", ctemp, parameters.facRe);
    temp = fscanf(file, "%[^,],%lf\n", ctemp, &parameters.facdt);
    printf("PARAMETERS: %s = %lf\n", ctemp, parameters.facdt);
    temp = fscanf(file, "%[^,],%d\n", ctemp, &temp1);
    if (temp1 == 1)
        parameters.fractional_step = true;
    else
        parameters.fractional_step = false;
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.fractional_step);
    if (temp == EOF){
        printf("Error: Unable to read the file\n");
        exit(1);
    }
    fclose(file);
    parameters.rho = 1;
    parameters.mu = parameters.rho / parameters.Re;
    parameters.nu = parameters.mu / parameters.rho;
}

void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels)
{   
    FILE *file;
    file = fopen(filename, "r");
    int fscan_temp;
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    fscan_temp = fscanf(file, "%[^,],%hd\n", temp, num_levels);
    if (fscan_temp == EOF){
        printf("Error: Unable to read the file\n");
        exit(1);
    }
    else
        printf("PARAMETERS: %s = %hd\n", temp, *num_levels);
    // Allocate memory for point structure for all levels and read meshfile names
    *myPointStruct = (PointStructure*)malloc((*num_levels) * sizeof(PointStructure));
    for (short ii = 0; ii<*num_levels ; ii = ii +1){
        fscan_temp = fscanf(file, "%s\n", (*myPointStruct)[ii].mesh_filename);
        if (fscan_temp == EOF){
            printf("Error: Unable to read the meshfile names\n");
            exit(1);
        }
    }
    fclose(file);
}

void calculate_parameters(PointStructure* myPointStruct) {
    if (parameters.dimension == 3)
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2) * (myPointStruct->poly_degree + 3)/ 6;
    else
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2)/ 2;
    myPointStruct->num_cloud_points = parameters.cloud_size_multiplier * myPointStruct->num_poly_terms;
    printf("PARAMETERS: num_poly_terms = %d\n", myPointStruct->num_poly_terms);
    printf("PARAMETERS: num_cloud_points = %d\n", myPointStruct->num_cloud_points);
    printf("PARAMETERS: dimension = %d\n", parameters.dimension);
    myPointStruct->pow_x = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_y = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_z = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    int count = 0;
    
    if (parameters.dimension==2)
        for (int i = 0; i <= myPointStruct->poly_degree; i++)
        {
            for (int j = 0; j <= myPointStruct->poly_degree - i; j++)
            {
                myPointStruct->pow_x[count] = i;
                myPointStruct->pow_y[count] = j;
                myPointStruct->pow_z[count] = 0;
                count++;
            }
        }
    else if (parameters.dimension==3)
        for (int i = 0; i <= myPointStruct->poly_degree; i++)
        {
            for (int j = 0; j <= myPointStruct->poly_degree - i; j++)
            {
                for (int k = 0; k <= myPointStruct->poly_degree - i - j; k++)
                {
                    myPointStruct->pow_x[count] = i;
                    myPointStruct->pow_y[count] = j;
                    myPointStruct->pow_z[count] = k;
                    count++;
                }
            }
        }
}

void read_PointStructure(PointStructure* myPointStruct)
{   
    FILE *file;
    int fscan_temp;
    char filename[50];
    strcpy(filename, myPointStruct->mesh_filename);

    file = fopen(filename, "r");
    double x_min = 1e6, y_min = 1e6, z_min = 1e6;
    double x_max = -1e6, y_max = -1e6, z_max = -1e6;
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }

    // Read the mesh file nodes
    while (true)
    {
        fscan_temp = fscanf(file, "%s ", temp);
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            if (strcmp(temp, "$Nodes") == 0)
                break;
    }
    fscan_temp = fscanf(file, "%i ", &itemp);

    // Allocate memory for the members of point structure
    AllocateMemoryPointStructure(myPointStruct, itemp);
    
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {
        fscan_temp = fscanf(file, "%i ", &itemp); //node number
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        fscan_temp =fscanf(file, "%lf ", &dtemp);
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            myPointStruct->x[i]=dtemp;
        fscan_temp = fscanf(file, "%lf ", &dtemp); 
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            myPointStruct->y[i]=dtemp;
        fscan_temp = fscanf(file, "%lf ", &dtemp); 
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            myPointStruct->z[i]=dtemp;
        
        myPointStruct->point_index[i] = i;
        myPointStruct->boundary_tag[i]=false;
        myPointStruct->corner_tag[i]=false;
        myPointStruct->x_normal[i]=0;
        myPointStruct->y_normal[i]=0;
        myPointStruct->z_normal[i]=0;

        if (myPointStruct->x[i] < x_min)
            x_min = myPointStruct->x[i];
        if (myPointStruct->y[i] < y_min)
            y_min = myPointStruct->y[i];
        if (myPointStruct->z[i] < z_min)
            z_min = myPointStruct->z[i];
        if (myPointStruct->x[i] > x_max)
            x_max = myPointStruct->x[i];
        if (myPointStruct->y[i] > y_max)
            y_max = myPointStruct->y[i];
        if (myPointStruct->z[i] > z_max)
            z_max = myPointStruct->z[i];
    }

    while (true)
    {
        fscan_temp = fscanf (file, "%s ", temp);
        if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            if (strcmp(temp, "$Elements") == 0)
                break;
    }
    fscan_temp=fscanf(file, "%i ", &itemp);
    if (fscan_temp == EOF){
            printf("Error: Unable to read the file\n");
            exit(1);
        }
        else
            myPointStruct->num_elem = itemp;

    if (parameters.dimension == 2){ //if all z coordinates are zero, then set dimension of the problem to 2
        myPointStruct->d_avg = sqrt((x_max - x_min)*(y_max - y_min) / myPointStruct->num_elem);
        printf("Average distance between nodes = %lf\n", myPointStruct->d_avg);
        }
    else {
        myPointStruct->d_avg = cbrt((x_max - x_min)*(y_max - y_min)*(z_max - z_min)/myPointStruct->num_elem);
        printf("Average distance between nodes = %lf\n", myPointStruct->d_avg);
    }

    // Read the mesh file elements and boundary nodes
    int e_type, tag_int, e_node1, e_node2, e_node3, e_node4; //variables for reading elements
    double dx,dy,dz,dx1,dy1,dz1; // vectors for calculating normal to the surface

    if (parameters.dimension==2) //if 2D problem, then read only line elements as boundary elements
    {
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscan_temp = fscanf(file, "%i ", &itemp);   //element number
            fscan_temp = fscanf(file, "%i ", &e_type); //type of the element
            if (e_type == 1) //reading line elements
            { 
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); // edge/boundary tag, if a square, elements on one edge will have same tag
                fscan_temp = fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1]; //calculate normal to the line
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1]; 
                myPointStruct->x_normal[e_node1 - 1] += -dy; myPointStruct->y_normal[e_node1 - 1] += dx; //add the normals to each node
                myPointStruct->x_normal[e_node2 - 1] += -dy; myPointStruct->y_normal[e_node2 - 1] += dx; //add the normals to each node
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscan_temp = fscanf(file, "%i ", &itemp);   //node number
                myPointStruct->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of elements
                fscan_temp = fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    else //if 3D problem, then read only surface elements as boundary elements
    {   
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscan_temp = fscanf(file, "%i ", &itemp);   //element number
            fscan_temp = fscanf(file, "%i ", &e_type); //element type
            if (e_type == 2)
            { //important 3 node triangle CV on the boundary
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); //element tag
                fscan_temp = fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node3);     //node number
                myPointStruct->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                //calculate normal to the surface of the triangle
                // by calculating cross product of the two vectors
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1]; 
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                dz = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
                dx1 = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                dy1 = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                dz1 = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];
                dtemp = dy*dz1 - dz*dy1;
                myPointStruct->x_normal[e_node1 - 1] += dtemp;
                myPointStruct->x_normal[e_node2 - 1] += dtemp;
                myPointStruct->x_normal[e_node3 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                myPointStruct->y_normal[e_node1 - 1] += dtemp;
                myPointStruct->y_normal[e_node2 - 1] += dtemp;
                myPointStruct->y_normal[e_node3 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                myPointStruct->z_normal[e_node1 - 1] += dtemp;
                myPointStruct->z_normal[e_node2 - 1] += dtemp;
                myPointStruct->z_normal[e_node3 - 1] += dtemp;
            }
            else if (e_type == 3)
            { //important 4 node quad element on the boundary
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); //element tag
                fscan_temp = fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node3);     //node number
                myPointStruct->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                fscan_temp = fscanf(file, "%i ", &e_node4);     //node number
                myPointStruct->boundary_tag[e_node4 - 1] = true; //set boundary tag to true
                // calculate normal to the surface of the quad element
                // by calculating cross product of any two vectors/lines on the surface
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1];
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                dz = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
                dx1 = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                dy1 = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                dz1 = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];
                dtemp = dy*dz1 - dz*dy1;
                myPointStruct->x_normal[e_node1 - 1] += dtemp;
                myPointStruct->x_normal[e_node2 - 1] += dtemp;
                myPointStruct->x_normal[e_node3 - 1] += dtemp;
                myPointStruct->x_normal[e_node4 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                myPointStruct->y_normal[e_node1 - 1] += dtemp;
                myPointStruct->y_normal[e_node2 - 1] += dtemp;
                myPointStruct->y_normal[e_node3 - 1] += dtemp;
                myPointStruct->y_normal[e_node4 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                myPointStruct->z_normal[e_node1 - 1] += dtemp;
                myPointStruct->z_normal[e_node2 - 1] += dtemp;
                myPointStruct->z_normal[e_node3 - 1] += dtemp;
                myPointStruct->z_normal[e_node4 - 1] += dtemp;
            }
            else if (e_type == 1)
            { //2 node line: log these 2 vertex numbers to delete them later
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); //element tag
                fscan_temp = fscanf(file, "%i ", &itemp);   //first point
                myPointStruct->corner_tag[itemp - 1] = true;
                fscan_temp = fscanf(file, "%i ", &itemp); //second point
                myPointStruct->corner_tag[itemp - 1] = true;
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &itemp);
                fscan_temp = fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscan_temp = fscanf(file, "%i ", &itemp);   //node number
                myPointStruct->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of element
                fscan_temp = fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    fclose(file);
    myPointStruct->num_corners = 0;
    myPointStruct->num_boundary_nodes = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {
        if (myPointStruct->boundary_tag[i] == true)
        {
            myPointStruct->num_boundary_nodes++;
        }
        if (myPointStruct->corner_tag[i] == true)
        {
            myPointStruct->num_corners++;
        }
    }
    // Remove first num_corners nodes
    int remove_count = myPointStruct->num_corners;

    memmove(myPointStruct->x, myPointStruct->x + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->y, myPointStruct->y + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->z, myPointStruct->z + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->x_normal, myPointStruct->x_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->y_normal, myPointStruct->y_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->z_normal, myPointStruct->z_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
    memmove(myPointStruct->point_index, myPointStruct->point_index + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(int));
    memmove(myPointStruct->boundary_tag, myPointStruct->boundary_tag + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(bool));
    memmove(myPointStruct->corner_tag, myPointStruct->corner_tag + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(bool));


    // Reallocate memory to shrink the array
    int new_size = myPointStruct->num_nodes - remove_count;
    myPointStruct->x = (double*)realloc(myPointStruct->x, new_size * sizeof(double));
    myPointStruct->y = (double*)realloc(myPointStruct->y, new_size * sizeof(double));
    myPointStruct->z = (double*)realloc(myPointStruct->z, new_size * sizeof(double));
    myPointStruct->x_normal = (double*)realloc(myPointStruct->x_normal, new_size * sizeof(double));
    myPointStruct->y_normal = (double*)realloc(myPointStruct->y_normal, new_size * sizeof(double));
    myPointStruct->z_normal = (double*)realloc(myPointStruct->z_normal, new_size * sizeof(double));
    myPointStruct->point_index = (int*)realloc(myPointStruct->point_index, new_size * sizeof(int));
    myPointStruct->boundary_tag = (bool*)realloc(myPointStruct->boundary_tag, new_size * sizeof(bool));
    myPointStruct->corner_tag = (bool*)realloc(myPointStruct->corner_tag, new_size * sizeof(bool));

    if (myPointStruct->x == NULL || myPointStruct->y == NULL || myPointStruct->z == NULL 
        || myPointStruct->x_normal == NULL || myPointStruct->y_normal == NULL || 
        myPointStruct->z_normal == NULL || myPointStruct->point_index == NULL || 
        myPointStruct->boundary_tag == NULL|| myPointStruct->corner_tag == NULL) {
        perror("realloc failed");
        exit(1);
    }
    myPointStruct->num_nodes = new_size;
    myPointStruct->num_boundary_nodes = myPointStruct->num_boundary_nodes - myPointStruct->num_corners;
    myPointStruct->num_corners = 0;

    for(int i = 0; i<new_size; i++){
        myPointStruct->point_index[i] = myPointStruct->point_index[i]-remove_count;
    }
    printf("No of nodes = %d \nNo of elements = %d \n", myPointStruct->num_nodes, myPointStruct->num_elem);
    printf("No of boundary nodes = %d \nNo of corner nodes = %d \n", myPointStruct->num_boundary_nodes, myPointStruct->num_corners);
}

void correct_normal_directions(PointStructure *myPointStruct)
{
    double dx, dy, dz, mag, dot_product, x_avg, y_avg, z_avg;
    short count = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {   
        if(myPointStruct->boundary_tag[i] == true)
        {   
            x_avg = 0; y_avg = 0; z_avg = 0; count = 0;
            // for (int j = 0; j < myPointStruct->num_cloud_points; j++){
            //     if ((myPointStruct->boundary_tag[myPointStruct->cloud_index[i][j]] == false)&&(myPointStruct->corner_tag[myPointStruct->cloud_index[i][j]] == false))
            //     {
            //         x_avg += myPointStruct->x[myPointStruct->cloud_index[i][j]];
            //         y_avg += myPointStruct->y[myPointStruct->cloud_index[i][j]];
            //         z_avg += myPointStruct->z[myPointStruct->cloud_index[i][j]];
            //         count++;
            //         // if (count == 10)
            //         //     break;
            //     }
            // }
            for (int j = myPointStruct->num_cloud_points-1; j > 0; j--){
                if (myPointStruct->boundary_tag[myPointStruct->cloud_index[i][j]] == false){
                    x_avg += myPointStruct->x[myPointStruct->cloud_index[i][j]];
                    y_avg += myPointStruct->y[myPointStruct->cloud_index[i][j]];
                    z_avg += myPointStruct->z[myPointStruct->cloud_index[i][j]];
                    count++;
                    if (count == 1)
                        break;
                }
            }
            x_avg = x_avg/count; y_avg = y_avg/count; z_avg = z_avg/count;

            dx = x_avg - myPointStruct->x[i];
            dy = y_avg - myPointStruct->y[i];
            dz = z_avg - myPointStruct->z[i];
            dot_product = dx * myPointStruct->x_normal[i] + dy * myPointStruct->y_normal[i] + dz * myPointStruct->z_normal[i];
            if (dot_product > 0)
            {
                myPointStruct->x_normal[i] = -myPointStruct->x_normal[i];
                myPointStruct->y_normal[i] = -myPointStruct->y_normal[i];
                myPointStruct->z_normal[i] = -myPointStruct->z_normal[i];
            }

            mag = sqrt(myPointStruct->x_normal[i] * myPointStruct->x_normal[i] + myPointStruct->y_normal[i] * myPointStruct->y_normal[i] + myPointStruct->z_normal[i] * myPointStruct->z_normal[i]);
            myPointStruct->x_normal[i] = myPointStruct->x_normal[i] / mag;
            myPointStruct->y_normal[i] = myPointStruct->y_normal[i] / mag;
            myPointStruct->z_normal[i] = myPointStruct->z_normal[i] / mag;
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
//////// Restriction and Prolongation matrix creation 

void create_restriction_matrix(PointStructure* myPointStruct_f, 
                                    PointStructure* myPointStruct_c)
{
    // allocate memory for restriction matrix
    
    short m = myPointStruct_f->num_cloud_points;
    short n = myPointStruct_f->num_poly_terms;
    short mpn = m+n;

    myPointStruct_c->restr_mat = (double**)malloc(myPointStruct_c->num_nodes * sizeof(double*));

    for (int i = 0; i < myPointStruct_c->num_nodes; i++)
        myPointStruct_c->restr_mat[i] = (double*)malloc(m * sizeof(double));

    // Initialise to zeros
    for (int i = 0; i < myPointStruct_c->num_nodes; i++)
    {
        for (short j = 0; j < m; j++)
        {
            myPointStruct_c->restr_mat[i][j] = 0;
        }
    }
    
    double **A_inv = create_matrix1(mpn,mpn);
    double **A = create_matrix1(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m); // coefficients of restriction matrix
    int i_restr;

    double point1[3], point2[3];
    for (int i = myPointStruct_c->num_boundary_nodes; 
                    i < myPointStruct_c->num_nodes; i++)
    {   
        i_restr = myPointStruct_c->restriction_points[i];
        create_A_matrix_from_cloud_indices(myPointStruct_f, A, myPointStruct_f->cloud_index[i_restr]);
        matrixInverse_Gauss_Jordan(A, A_inv, m+n);
        
        point1[0] = myPointStruct_c->x[i];
        point1[1] = myPointStruct_c->y[i];
        point1[2] = myPointStruct_c->z[i];
        
        for (short j = 0; j < m; j++) {
            point2[0] = myPointStruct_f->x[myPointStruct_f->cloud_index[i_restr][j]];
            point2[1] = myPointStruct_f->y[myPointStruct_f->cloud_index[i_restr][j]];
            point2[2] = myPointStruct_f->z[myPointStruct_f->cloud_index[i_restr][j]];
            temp[j] = calculate_phs_rbf(point1, point2, parameters.phs_degree, 
                                                        parameters.dimension);
        }
        double seed_pt[3] = {myPointStruct_f->x[i_restr],
                             myPointStruct_f->y[i_restr],
                             myPointStruct_f->z[i_restr]};
        for (short j = 0; j < n; j++) 
            temp[j+m] = pow(point1[0]-seed_pt[0], myPointStruct_f->pow_x[j])
                                * pow(point1[1]-seed_pt[1], myPointStruct_f->pow_y[j])
                                * pow(point1[2]-seed_pt[2], myPointStruct_f->pow_z[j]);
        
        multiply_vector_matrix_columnwise(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < m; j++)
        {
            myPointStruct_c->restr_mat[i][j] = temp1[j];
        }
    }
    free(temp);
    free(temp1);
    free_matrix(A, mpn);
    free_matrix(A_inv, mpn);
}

void create_prolongation_matrix(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c)
{
    // allocate memory for restriction matrix
    
    short m = myPointStruct_c->num_cloud_points;
    short n = myPointStruct_c->num_poly_terms;
    short mpn = m+n;
    
    myPointStruct_f->prol_mat = (double**)malloc(myPointStruct_f->num_nodes * sizeof(double*));
    
    for (int i = 0;
                     i < myPointStruct_f->num_nodes; i++)
        myPointStruct_f->prol_mat[i] = (double*)malloc(m * sizeof(double));

    // Initialise to zeros
    for (int i = 0;
                     i < myPointStruct_f->num_nodes; i++)
    {
        for (short j = 0; j < m; j++)
        {
            myPointStruct_f->prol_mat[i][j] = 0;
        }
    }
    
    double **A_inv = create_matrix1(mpn,mpn);
    double **A = create_matrix1(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m);
    int i_prol;

    double point1[3], point2[3], seed_pt[3];
    for (int i = myPointStruct_f->num_boundary_nodes;
                     i < myPointStruct_f->num_nodes; i++)
    {   
        i_prol = myPointStruct_f->prolongation_points[i];
        create_A_matrix_from_cloud_indices(myPointStruct_c, A, myPointStruct_c->cloud_index[i_prol]);
        matrixInverse_Gauss_Jordan(A, A_inv, m+n);
        
        point1[0] = myPointStruct_f->x[i];
        point1[1] = myPointStruct_f->y[i];
        point1[2] = myPointStruct_f->z[i];
        
        for (short j = 0; j < m; j++) {
            point2[0] = myPointStruct_c->x[myPointStruct_c->cloud_index[i_prol][j]];
            point2[1] = myPointStruct_c->y[myPointStruct_c->cloud_index[i_prol][j]];
            point2[2] = myPointStruct_c->z[myPointStruct_c->cloud_index[i_prol][j]];
            temp[j] = calculate_phs_rbf(point1, point2, parameters.phs_degree, 
                                                        parameters.dimension);
        }
        seed_pt[0]= myPointStruct_c->x[i_prol];
        seed_pt[1]= myPointStruct_c->y[i_prol];
        seed_pt[2]= myPointStruct_c->z[i_prol];
        for (short j = 0; j < n; j++) 
            temp[j+m] = pow(point1[0]-seed_pt[0], myPointStruct_c->pow_x[j])
                                * pow(point1[1]-seed_pt[1], myPointStruct_c->pow_y[j])
                                * pow(point1[2]-seed_pt[2], myPointStruct_c->pow_z[j]);
        
        multiply_vector_matrix_columnwise(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < (m); j++)
        {
            myPointStruct_f->prol_mat[i][j] = temp1[j];
        }
    }
    free(temp);
    free(temp1);
    free_matrix(A, mpn);
    free_matrix(A_inv, mpn);
}

void rcm_reordering(PointStructure* myPointstruct) {
    int *queue1 = (int*)malloc(myPointstruct->num_nodes * sizeof(int));  // RCM ordered indices
    int *queue2 = (int*)malloc(myPointstruct->num_nodes * sizeof(int));  // Maps from original to RCM
    bool *visited = (bool*)malloc(myPointstruct->num_nodes * sizeof(bool));

    // Temporary arrays to hold reordered data
    double *temp_x = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_y = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_z = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_nx = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_ny = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_nz = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    
    int **temp_cloud_index = (int**)malloc(myPointstruct->num_nodes * sizeof(int*));
    for (int i = 0; i < myPointstruct->num_nodes; i++)
        temp_cloud_index[i] = (int*)malloc(myPointstruct->num_cloud_points * sizeof(int));

    // Initialize visited array and queue counters
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        visited[i] = false;
    }

    int count = 0;        // Counter for filling queue1 (RCM order)
    int count_queue = 0;   // Pointer to track current node in BFS

    // Start BFS with node 0 (or a node of your choice)
    for (int i = 0; i<myPointstruct->num_boundary_nodes; i++){
        queue1[count++] = i;   // Starting with node 0 (can change to another node)
        count_queue++;
        visited[i] = true;
    }
    queue1[count++] = myPointstruct->point_index[myPointstruct->num_boundary_nodes];   // Starting with node 0 (can change to another node)
    visited[myPointstruct->num_boundary_nodes] = true;

    // Perform BFS-like traversal to generate RCM order in queue1
    while (count_queue < count) {
        int current = queue1[count_queue++];  // Dequeue node

        // For each neighbor (cloud point) of the current node
        for (int i = 1; i < myPointstruct->num_cloud_points; i++) {
            int neighbor = myPointstruct->cloud_index[current][i];  // Get neighboring node

            // Enqueue neighbor if not visited
            if (!visited[neighbor]) {
                queue1[count++] = neighbor;
                visited[neighbor] = true;
            }
        }
    }

    // Now queue1 contains nodes in the RCM order. Generate queue2 for reverse mapping
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        queue2[queue1[i]] = i;  // Mapping original node to RCM index
    }

    // Reorder the data based on RCM order (queue1) and apply the reverse mapping (queue2)
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        temp_x[i] = myPointstruct->x[queue1[i]];
        temp_y[i] = myPointstruct->y[queue1[i]];
        temp_z[i] = myPointstruct->z[queue1[i]];

        // Reorder the cloud index array based on queue2 mapping
        for (int j = 0; j < myPointstruct->num_cloud_points; j++) {
            temp_cloud_index[i][j] = queue2[myPointstruct->cloud_index[queue1[i]][j]];
        }

        temp_nx[i] = myPointstruct->x_normal[queue1[i]];
        temp_ny[i] = myPointstruct->y_normal[queue1[i]];
        temp_nz[i] = myPointstruct->z_normal[queue1[i]];
    }

    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        myPointstruct->x[i] = temp_x[i];
        myPointstruct->y[i] = temp_y[i];
        myPointstruct->z[i] = temp_z[i];

        for (int j = 0; j < myPointstruct->num_cloud_points; j++) {
            myPointstruct->cloud_index[i][j] = temp_cloud_index[i][j];
        }

        myPointstruct->x_normal[i] = temp_nx[i];
        myPointstruct->y_normal[i] = temp_ny[i];
        myPointstruct->z_normal[i] = temp_nz[i];
    }

    // Free allocated memory
    free(queue1);
    free(queue2);
    free(visited);
    free(temp_x);
    free(temp_y);
    free(temp_z);
    free(temp_nx);
    free(temp_ny);
    free(temp_nz);

    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        free(temp_cloud_index[i]);
    }
    free(temp_cloud_index);
}

void rcm_reordering_with_boundarynodes(PointStructure* myPointstruct) {
    int *queue1 = (int*)malloc(myPointstruct->num_nodes * sizeof(int));  // RCM ordered indices
    int *queue2 = (int*)malloc(myPointstruct->num_nodes * sizeof(int));  // Maps from original to RCM
    bool *visited = (bool*)malloc(myPointstruct->num_nodes * sizeof(bool));

    // Temporary arrays to hold reordered data
    double *temp_x = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_y = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_z = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_nx = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_ny = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    double *temp_nz = (double*)malloc(myPointstruct->num_nodes * sizeof(double));
    bool *temp_boundary_tag = (bool*)malloc(myPointstruct->num_nodes * sizeof(bool));
    
    int **temp_cloud_index = (int**)malloc(myPointstruct->num_nodes * sizeof(int*));
    for (int i = 0; i < myPointstruct->num_nodes; i++)
        temp_cloud_index[i] = (int*)malloc(myPointstruct->num_cloud_points * sizeof(int));

    // Initialize visited array and queue counters
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        visited[i] = false;
    }

    int count = 0;        // Counter for filling queue1 (RCM order)
    int count_queue = 0;   // Pointer to track current node in BFS

    queue1[0] = 0;   // Starting with node 0 (can change to another node)
    visited[0] = true;
    count++;
    
    // Perform BFS-like traversal to generate RCM order in queue1
    while (count_queue < count) {
        int current = queue1[count_queue++];  // Dequeue node

        // For each neighbor (cloud point) of the current node
        for (int i = 1; i < myPointstruct->num_cloud_points; i++) {
            int neighbor = myPointstruct->cloud_index[current][i];  // Get neighboring node

            // Enqueue neighbor if not visited
            if (!visited[neighbor]) {
                queue1[count++] = neighbor;
                visited[neighbor] = true;
            }
        }
    }

    // Now queue1 contains nodes in the RCM order. Generate queue2 for reverse mapping
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        queue2[queue1[i]] = i;  // Mapping original node to RCM index
    }

    // Reorder the data based on RCM order (queue1) and apply the reverse mapping (queue2)
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        temp_x[i] = myPointstruct->x[queue1[i]];
        temp_y[i] = myPointstruct->y[queue1[i]];
        temp_z[i] = myPointstruct->z[queue1[i]];

        // Reorder the cloud index array based on queue2 mapping
        for (int j = 0; j < myPointstruct->num_cloud_points; j++) {
            temp_cloud_index[i][j] = queue2[myPointstruct->cloud_index[queue1[i]][j]];
        }

        temp_nx[i] = myPointstruct->x_normal[queue1[i]];
        temp_ny[i] = myPointstruct->y_normal[queue1[i]];
        temp_nz[i] = myPointstruct->z_normal[queue1[i]];
        temp_boundary_tag[i] = myPointstruct->boundary_tag[queue1[i]];
    }

    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        myPointstruct->x[i] = temp_x[i];
        myPointstruct->y[i] = temp_y[i];
        myPointstruct->z[i] = temp_z[i];

        for (int j = 0; j < myPointstruct->num_cloud_points; j++) {
            myPointstruct->cloud_index[i][j] = temp_cloud_index[i][j];
        }

        myPointstruct->x_normal[i] = temp_nx[i];
        myPointstruct->y_normal[i] = temp_ny[i];
        myPointstruct->z_normal[i] = temp_nz[i];
        myPointstruct->boundary_tag[i] = temp_boundary_tag[i];
    }

    // Free allocated memory
    free(queue1);
    free(queue2);
    free(visited);
    free(temp_x);
    free(temp_y);
    free(temp_z);
    free(temp_nx);
    free(temp_ny);
    free(temp_nz);
    free(temp_boundary_tag);
    for (int i = 0; i < myPointstruct->num_nodes; i++) {
        free(temp_cloud_index[i]);
    }
    free(temp_cloud_index);
}

void create_prolongation_and_restriction_matrices(PointStructure* myPointStruct, short num_levels){
    
    printf("\nIdentifying prolongation and restriction points\n");
    for (short ii = 0; ii<num_levels ; ii = ii +1){
        printf("Level %d\n", ii+1);
        if (ii == num_levels-1)
            myPointStruct[ii].prolongation_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii], 1);
        else{
            myPointStruct[ii].prolongation_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii+1], 1);
            create_prolongation_matrix(&myPointStruct[ii], &myPointStruct[ii+1]);
        }
        if (ii == 0)
            myPointStruct[ii].restriction_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii], 1);  
        else{
            myPointStruct[ii].restriction_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii-1], 1);
            create_restriction_matrix(&myPointStruct[ii-1], &myPointStruct[ii]);
        }
    }
    printf("Prolongation and restriction points identified\n");
}

void calculate_avg_dx(PointStructure* myPointStruct){
    double dx, dy, dz;
    myPointStruct->d_avg = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            dx = myPointStruct->x[i] - myPointStruct->x[myPointStruct->cloud_index[i][1]];
            dy = myPointStruct->y[i] - myPointStruct->y[myPointStruct->cloud_index[i][1]];
            if (parameters.dimension == 3)
                dz = myPointStruct->z[i] - myPointStruct->z[myPointStruct->cloud_index[i][1]];
            else
                dz = 0;
            myPointStruct->d_avg += sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    myPointStruct->d_avg = myPointStruct->d_avg/myPointStruct->num_nodes;
}

///////////////////////////////////////////////////////////////////////////////

void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels)
{   
    for (short ii = 0; ii<num_levels ; ii = ii +1){
        myPointStruct[ii].poly_degree = 3; 
    	myPointStruct[0].poly_degree = parameters.poly_degree;
        printf("\nReading mesh data for level %d\n", ii+1);
        calculate_parameters(&myPointStruct[ii]);
        read_PointStructure(&myPointStruct[ii]);
        find_cloud_index(&myPointStruct[ii]);
        rcm_reordering_with_boundarynodes(&myPointStruct[ii]);
        correct_normal_directions(&myPointStruct[ii]);
        calculate_avg_dx(&myPointStruct[ii]);
    }
    create_prolongation_and_restriction_matrices(myPointStruct, num_levels);
}

#endif
