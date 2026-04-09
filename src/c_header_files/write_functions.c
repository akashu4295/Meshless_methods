// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// Functions used to write the output to files

#include "functions.h"

//////////////////////////////////////////////////////////////////////
// Function Definitions
//////////////////////////////////////////////////////////////////////

void write_normals(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing normals to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %lf %lf %lf\n", i, myPointStruct->x_normal[i], myPointStruct->y_normal[i], myPointStruct->z_normal[i]);
    fclose(file);
}

void write_boundary_tags(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing boundary tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %d\n", i, myPointStruct->boundary_tag[i]);
    fclose(file);
}

void write_corner_tags(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing corner tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %d\n", i, myPointStruct->corner_tag[i]);
    fclose(file);
}

void write_coordinates(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing coordinates to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%lf %lf %lf\n", myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i]);
    fclose(file);
}

void write_cloud_index(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing cloud index to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }

    int n = myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fprintf(file, "%d ", i);
        for (int j = 0; j < n; j++)
            fprintf(file, "%d ", myPointStruct->cloud_index[i*n +j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_prolongation_and_restriction_points(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing prolongation and restriction points to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fprintf(file, "%d %d %d\n", i, myPointStruct->prolongation_points[i], myPointStruct->restriction_points[i]);
    }
    fclose(file);
}

void write_test_files(double* f, double* fx, double* fy, double* fz, double* lapf, double* fxx, double* fyy, double* fzz, int num_nodes, char* folder1)
{
    FILE *file;
    char temp[100];
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"f.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", f[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fx.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fx[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fy.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fy[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fz.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fz[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"lapf.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", lapf[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fxx.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fxx[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fyy.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fyy[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fzz.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fzz[i]);
    }
    fclose(file);
    printf("Files written\n");
}

void write_processed_grid_data(PointStructure* myPointStruct, int ii)
{   
    char filename[50];
    sprintf(filename, "normals_%d.csv", ii);
    write_normals(myPointStruct, filename); // Write normals of all points
    sprintf(filename, "boundary_tags_%d.csv", ii);
    write_boundary_tags(myPointStruct, filename); // Write boundary tags of all points
    sprintf(filename, "corner_tags_%d.csv", ii);
    write_corner_tags(myPointStruct, filename); // Write corner tags of all points
    sprintf(filename, "coordinates_%d.csv", ii);
    write_coordinates(myPointStruct, filename); // Write coordinates of all points
    sprintf(filename, "cloud_index_%d.csv", ii);
    write_cloud_index(myPointStruct, filename); // Write coordinates of all points
    sprintf(filename, "prolongation_and_restriction_%d.csv", ii);
    write_prolongation_and_restriction_points(myPointStruct, filename);
    printf("\n\n");
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Structure to hold field data (velocities and pressure)
typedef struct {
    double *u;  // x-velocity
    double *v;  // y-velocity
    double *w;  // z-velocity
    double *p;  // pressure
} Field;

int write_vtk(char *gmsh_filename, FieldVariables *field, PointStructure* myPS, char *sol_filename)
{
    FILE *fp_in, *fp_out;
    char line[256];
    char vtk_filename[256];
    sprintf(vtk_filename, "%s.vtk", sol_filename);

    int num_nodes = 0, num_elements = 0;
    int i, node_id;
    double x, y, z;

    fp_in = fopen(gmsh_filename, "r");
    if (!fp_in) {
        fprintf(stderr, "Error: Cannot open %s\n", gmsh_filename);
        return -1;
    }

    /* ---------------- READ NODES ---------------- */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Nodes")) {
            fscanf(fp_in, "%d", &num_nodes);
            break;
        }
    }

    double *nodes_x = malloc(num_nodes * sizeof(double));
    double *nodes_y = malloc(num_nodes * sizeof(double));
    double *nodes_z = malloc(num_nodes * sizeof(double));

    for (i = 0; i < num_nodes; i++) {
        fscanf(fp_in, "%d %lf %lf %lf", &node_id, &x, &y, &z);
        nodes_x[node_id - 1] = x;
        nodes_y[node_id - 1] = y;
        nodes_z[node_id - 1] = z;
    }

    /* ---------------- READ ELEMENTS ---------------- */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Elements")) {
            fscanf(fp_in, "%d", &num_elements);
            break;
        }
    }

    int max_conn = (parameters.dimension == 2) ? 3 : 4;
    int vtk_type  = (parameters.dimension == 2) ? 5 : 10;

    int *conn = malloc(num_elements * max_conn * sizeof(int));
    int cell_count = 0;

    for (i = 0; i < num_elements; i++) {
        int elem_id, elem_type, num_tags;

        fscanf(fp_in, "%d %d %d", &elem_id, &elem_type, &num_tags);

        for (int j = 0; j < num_tags; j++) {
            int tmp;
            fscanf(fp_in, "%d", &tmp);
        }

        if (parameters.dimension == 2 && elem_type == 2) {
            int n1, n2, n3;
            fscanf(fp_in, "%d %d %d", &n1, &n2, &n3);

            conn[cell_count*3+0] = n1-1;
            conn[cell_count*3+1] = n2-1;
            conn[cell_count*3+2] = n3-1;
            cell_count++;
        }
        else if (parameters.dimension == 3 && elem_type == 4) {
            int n1, n2, n3, n4;
            fscanf(fp_in, "%d %d %d %d", &n1, &n2, &n3, &n4);

            conn[cell_count*4+0] = n1-1;
            conn[cell_count*4+1] = n2-1;
            conn[cell_count*4+2] = n3-1;
            conn[cell_count*4+3] = n4-1;
            cell_count++;
        }
        else {
            fgets(line, sizeof(line), fp_in);
        }
    }

    fclose(fp_in);

    /* ---------------- WRITE VTK ---------------- */
    fp_out = fopen(vtk_filename, "w");

    fprintf(fp_out, "# vtk DataFile Version 3.0\n");
    fprintf(fp_out, "Mesh\n");
    fprintf(fp_out, "ASCII\n");
    fprintf(fp_out, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp_out, "POINTS %d double\n", num_nodes);
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%lf %lf %lf\n",
                nodes_x[i], nodes_y[i], nodes_z[i]);
    }

    int vtk_cell_size = max_conn + 1;

    fprintf(fp_out, "\nCELLS %d %d\n",
            cell_count,
            cell_count * vtk_cell_size);

    for (i = 0; i < cell_count; i++) {
        fprintf(fp_out, "%d ", max_conn);
        for (int j = 0; j < max_conn; j++) {
            fprintf(fp_out, "%d ", conn[i*max_conn + j]);
        }
        fprintf(fp_out, "\n");
    }

    fprintf(fp_out, "\nCELL_TYPES %d\n", cell_count);
    for (i = 0; i < cell_count; i++) {
        fprintf(fp_out, "%d\n", vtk_type);
    }

    /* ---------------- POINT DATA ---------------- */
    fprintf(fp_out, "\nPOINT_DATA %d\n", num_nodes);

    /* Velocity */
    fprintf(fp_out, "VECTORS velocity double\n");

    int num_corner = myPS->num_corners;

    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "0.0 0.0 0.0\n");
    }

    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i - num_corner];

        if (parameters.dimension == 2) {
            fprintf(fp_out, "%.16e %.16e %.16e\n",
                    field->u[k], field->v[k], 0.0);
        } else {
            fprintf(fp_out, "%.16e %.16e %.16e\n",
                    field->u[k], field->v[k], field->w[k]);
        }
    }

    /* Pressure */
    fprintf(fp_out, "\nSCALARS pressure double 1\n");
    fprintf(fp_out, "LOOKUP_TABLE default\n");

    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "0.0\n");
    }

    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i - num_corner];
        fprintf(fp_out, "%.16e\n", field->p[k]);
    }

    fclose(fp_out);

    free(nodes_x);
    free(nodes_y);
    free(nodes_z);
    free(conn);

    printf("VTK written: %d cells\n", cell_count);

    return 0;
}
