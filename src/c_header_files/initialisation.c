#include "structures.h"
#include "functions.h"
#include "kdtree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* trim whitespace in-place */
static char* trim(char* s){
    while (isspace((unsigned char)*s)) s++;
    if (*s == 0) return s;

    char* end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return s;
}

BCType parse_bc_type(const char* s){
    if (!strcmp(s, "wall")) return BC_WALL;
    if (!strcmp(s, "velocity_inlet")) return BC_VELOCITY_INLET;
    if (!strcmp(s, "velocity_outlet")) return BC_VELOCITY_OUTLET;
    if (!strcmp(s, "pressure_outlet")) return BC_PRESSURE_OUTLET;
    if (!strcmp(s, "isothermal_wall")) return BC_ISOTHERMAL_WALL;
    if (!strcmp(s, "adiabatic_wall")) return BC_ADIABATIC_WALL;
    if (!strcmp(s, "symmetry")) return BC_SYMMETRY;
    if (!strcmp(s, "supersonic_inlet")) return BC_SUPERSONIC_INLET;
    if (!strcmp(s, "supersonic_outlet")) return BC_SUPERSONIC_OUTLET;
    if (!strcmp(s, "subsonic_inlet")) return BC_SUBSONIC_INLET;
    if (!strcmp(s, "subsonic_outlet")) return BC_SUBSONIC_OUTLET;
    if (!strcmp(s, "interior")) return BC_INTERIOR;
    return BC_DEFAULT;
}

void read_physical_names(char* meshfile, PointStructure* ps)
{
    FILE* file_pn = fopen(meshfile, "r");
    if (file_pn == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    char temp[128];
    int n;
    while (fscanf(file_pn, "%s", temp) != EOF) {
        if (strcmp(temp, "$PhysicalNames") == 0)
            break;
    }
    fscanf(file_pn, "%d", &n);
    printf("Number of physical names: %d\n", n);
    ps->num_boundary_types = n;
    if (n == 0) {
        ps->boundary_map[0].physical_id = 0;
        strcpy(ps->boundary_map[0].name, "default");
        ps->boundary_map[0].bc.type = BC_INTERIOR;
        ps->boundary_map[0].bc.u = 0.0;
        ps->boundary_map[0].bc.v = 0.0;
        ps->boundary_map[0].bc.w = 0.0;
        ps->boundary_map[0].bc.p = 0.0;
        ps->boundary_map[0].bc.T = parameters.T_ref;
        ps->boundary_map[0].bc.rho = parameters.rho_ref;
        printf("No physical names found in mesh file.\n");
        fclose(file_pn);
        return;
    }
    for (int i = 0; i < n; i++) {
        int dim;
        char raw[64];
        fscanf(file_pn, "%d %hd %63s", &dim, &ps->boundary_map[i].physical_id, raw);
        int len = strlen(raw);
        if (raw[0] == '"' && raw[len-1] == '"') {
            raw[len-1] = '\0';
            strcpy(ps->boundary_map[i].name, raw + 1);
        } else {
            strcpy(ps->boundary_map[i].name, raw);
        }
        printf("\tPhysical ID: %d, Name: %s\n", ps->boundary_map[i].physical_id, ps->boundary_map[i].name);
        
        // Initialize with default values
        ps->boundary_map[i].bc.type = BC_INTERIOR;
        ps->boundary_map[i].bc.u = 0.0;
        ps->boundary_map[i].bc.v = 0.0;
        ps->boundary_map[i].bc.w = 0.0;
        ps->boundary_map[i].bc.p = parameters.p_ref;
        ps->boundary_map[i].bc.T = parameters.T_ref;
        ps->boundary_map[i].bc.rho = parameters.rho_ref;
        ps->boundary_map[i].bc.p_total = 0.0;
        ps->boundary_map[i].bc.T_total = 0.0;
    }
    fclose(file_pn);
}

void read_boundary_conditions_file(char* bcfile, PointStructure* ps){
    FILE* file = fopen(bcfile, "r");
    char line[512];

    if (!file) {
        perror("Cannot open BC file");
        exit(1);
    }

    int lineno = 0;

    while (fgets(line, sizeof(line), file)) {
        lineno++;

        /* strip newline */
        line[strcspn(line, "\r\n")] = 0;

        char* s = trim(line);

        /* skip empty lines and comments */
        if (*s == '\0' || *s == '#')
            continue;

        BCValue bc;
        bc.type = BC_DEFAULT;
        bc.u = 0.0; bc.v = 0.0; bc.w = 0.0; bc.p = 0.0; bc.v_n = 0; bc.v_t = 0;
        if (parameters.compressible_flow){
            bc.p_total = parameters.p_ref;
            bc.rho = parameters.rho_ref;
            bc.T = parameters.T_ref;
            bc.T_total = parameters.T_ref;
        }

        // Set defaults for compressible flow
        bc.T = parameters.T_ref;
        bc.rho = parameters.rho_ref;
        bc.p = parameters.p_ref;
        bc.p_total = 0.0;
        bc.T_total = 0.0;

        char* tokens[16];
        int ntok = 0;

        /* tokenize by commas */
        char* tok = strtok(s, ",");
        while (tok && ntok < 16) {
            tokens[ntok++] = trim(tok);
            tok = strtok(NULL, ",");
        }

        if (ntok < 2) {
            printf("BC CSV error (line %d): need at least name,type\n", lineno);
            continue;
        }

        const char* name = tokens[0];
        const char* type = tokens[1];

        bc.type = parse_bc_type(type);
        if (bc.type == BC_DEFAULT) {
            printf("BC CSV error (line %d): unknown BC type '%s'\n",
                   lineno, type);
            continue;
        }

        int have_u = 0, have_v = 0, have_w = 0, have_p = 0, have_vn = 0, have_vt = 0;
        int have_T = 0, have_rho = 0, have_p_total = 0, have_T_total = 0;

        /* parse key=value pairs */
        for (int i = 2; i < ntok; i++) {
            char* eq = strchr(tokens[i], '=');
            if (!eq) {
                printf("BC CSV warning (line %d): ignoring token '%s'\n",
                       lineno, tokens[i]);
                continue;
            }

            *eq = 0;
            char* key = trim(tokens[i]);
            char* val = trim(eq + 1);

            if (!strcmp(key, "u")) {
                bc.u = atof(val); have_u = 1;
            }
            else if (!strcmp(key, "v")) {
                bc.v = atof(val); have_v = 1;
            }
            else if (!strcmp(key, "w")) {
                bc.w = atof(val); have_w = 1;
            }
            else if (!strcmp(key, "v_n")) {
                bc.v_n = atof(val); have_vn = 1;
            }
            else if (!strcmp(key, "v_t")) {
                bc.v_t = atof(val); have_vt = 1;
            }
            else if (!strcmp(key, "p")) {
                bc.p = atof(val); have_p = 1;
            }
            else if (!strcmp(key, "T")) {
                bc.T = atof(val); have_T = 1;
            }
            else if (!strcmp(key, "rho")) {
                bc.rho = atof(val); have_rho = 1;
            }
            else if (!strcmp(key, "p_total") || !strcmp(key, "ptotal")) {
                bc.p_total = atof(val); have_p_total = 1;
            }
            else if (!strcmp(key, "T_total") || !strcmp(key, "Ttotal")) {
                bc.T_total = atof(val); have_T_total = 1;
            }
            else {
                printf("BC CSV warning (line %d): unknown key '%s'\n",
                       lineno, key);
            }
        }

        /* semantic validation */
        if ((bc.type == BC_VELOCITY_INLET ||
             bc.type == BC_VELOCITY_OUTLET ||
             bc.type == BC_WALL) &&
            (!(have_u && have_v && have_w) &&
            !(have_vn && have_vt)))
        {
            printf("BC CSV error (line %d): velocity BC requires u,v,w or v_n,v_t\n", lineno);
            continue;
        }

        if (bc.type == BC_PRESSURE_OUTLET && !have_p) {
            printf("BC CSV error (line %d): pressure BC requires p\n", lineno);
            continue;
        }
        
        // Compressible flow validations
        if (parameters.compressible_flow) {
            // For velocity inlets in compressible flow, need temperature
            if (bc.type == BC_VELOCITY_INLET && !have_T && !have_rho) {
                printf("BC CSV warning (line %d): velocity inlet in compressible flow should specify T or rho (using T_ref)\n", lineno);
                bc.T = parameters.T_ref;
            }
            
            // For isothermal walls, need temperature
            if (bc.type == BC_ISOTHERMAL_WALL && !have_T) {
                printf("BC CSV error (line %d): isothermal wall requires T\n", lineno);
                continue;
            }
            
            // For subsonic inlet, need total conditions
            if (bc.type == BC_SUBSONIC_INLET && (!have_p_total || !have_T_total)) {
                printf("BC CSV error (line %d): subsonic inlet requires p_total and T_total\n", lineno);
                continue;
            }
            
            // Calculate density from ideal gas if not specified
            if (have_T && have_p && !have_rho) {
                bc.rho = bc.p / (parameters.R_gas * bc.T);
            } else if (have_rho && have_T && !have_p) {
                bc.p = bc.rho * parameters.R_gas * bc.T;
            }
        }

        /* map BC to physical name */
        int found = 0;
        for (int i = 0; i < ps->num_boundary_types; i++) {
            if (!strcmp(ps->boundary_map[i].name, name)) {
                ps->boundary_map[i].bc = bc;
                found = 1;
                break;
            }
        }

        if (!found) {
            printf("BC CSV warning (line %d): physical name '%s' not found in mesh\n",
                   lineno, name);
        }
    }
    fclose(file);
}

int bc_priority(BCType t){
    if (t == BC_VELOCITY_INLET) return 3;
    if (t == BC_WALL) return 2;
    if (t == BC_PRESSURE_OUTLET) return 1;
    if (t == BC_ISOTHERMAL_WALL) return 2;
    if (t == BC_ADIABATIC_WALL) return 2;
    return 0;
}

void assign_node_bc(PointStructure* ps, int node, BCValue new_bc){
    if (bc_priority(new_bc.type) > bc_priority(ps->node_bc[node].type)) 
        ps->node_bc[node] = new_bc;
}
static void decompose_velocity(double v_n, double v_t, double nx, double ny,
                                double* u_out, double* v_out)
{
    *u_out = v_n * nx - v_t * ny;
    *v_out = v_n * ny + v_t * nx;
}

void apply_boundary_conditions(PointStructure* myPointStruct, FieldVariables* field, int numlevels){
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            if (!myPointStruct[ii].boundary_tag[i]) continue;
            
            BCType bc_type = myPointStruct[ii].node_bc[i].type;
            
            switch (bc_type)
            {
                case BC_VELOCITY_INLET:
                case BC_VELOCITY_OUTLET:
                {
                    double vn = myPointStruct[ii].node_bc[i].v_n;
                    double vt = myPointStruct[ii].node_bc[i].v_t;

                    if (vn != 0.0 || vt != 0.0)
                    {
                        double nx = myPointStruct[ii].x_normal[i];
                        double ny = myPointStruct[ii].y_normal[i];
                        decompose_velocity(vn, vt, nx, ny,
                                           &field[ii].u[i], &field[ii].v[i]);
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    else
                    {
                        field[ii].u[i] = myPointStruct[ii].node_bc[i].u;
                        field[ii].v[i] = myPointStruct[ii].node_bc[i].v;
                        if (parameters.dimension == 3)
                            field[ii].w[i] = myPointStruct[ii].node_bc[i].w;
                    }

                    if (parameters.compressible_flow) {
                        field[ii].T[i]   = myPointStruct[ii].node_bc[i].T;
                        field[ii].rho[i] = myPointStruct[ii].node_bc[i].rho;
                        field[ii].p[i]   = myPointStruct[ii].node_bc[i].p;
                    }
                    break;
                }

                case BC_WALL:
                {
                    double vn = myPointStruct[ii].node_bc[i].v_n;
                    double vt = myPointStruct[ii].node_bc[i].v_t;

                    if (vn != 0.0 || vt != 0.0)
                    {
                        double nx = myPointStruct[ii].x_normal[i];
                        double ny = myPointStruct[ii].y_normal[i];
                        decompose_velocity(vn, vt, nx, ny,
                                           &field[ii].u[i], &field[ii].v[i]);
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    else
                    {
                        /* Stationary no-slip wall */
                        field[ii].u[i] = 0.0;
                        field[ii].v[i] = 0.0;
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    break;
                }

                case BC_ISOTHERMAL_WALL:
                {
                    double vn = myPointStruct[ii].node_bc[i].v_n;
                    double vt = myPointStruct[ii].node_bc[i].v_t;

                    if (vn != 0.0 || vt != 0.0)
                    {
                        double nx = myPointStruct[ii].x_normal[i];
                        double ny = myPointStruct[ii].y_normal[i];
                        decompose_velocity(vn, vt, nx, ny,
                                           &field[ii].u[i], &field[ii].v[i]);
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    else
                    {
                        field[ii].u[i] = 0.0;
                        field[ii].v[i] = 0.0;
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }

                    if (parameters.compressible_flow) {
                        field[ii].T[i]   = myPointStruct[ii].node_bc[i].T;
                        field[ii].rho[i] = parameters.p_ref /
                            (parameters.R_gas * field[ii].T[i]);
                    }
                    break;
                }

                case BC_ADIABATIC_WALL:
                {
                    double vn = myPointStruct[ii].node_bc[i].v_n;
                    double vt = myPointStruct[ii].node_bc[i].v_t;

                    if (vn != 0.0 || vt != 0.0)
                    {
                        double nx = myPointStruct[ii].x_normal[i];
                        double ny = myPointStruct[ii].y_normal[i];
                        decompose_velocity(vn, vt, nx, ny,
                                           &field[ii].u[i], &field[ii].v[i]);
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    else
                    {
                        field[ii].u[i] = 0.0;
                        field[ii].v[i] = 0.0;
                        if (parameters.dimension == 3)
                            field[ii].w[i] = 0.0;
                    }
                    /* Temperature gradient enforced in solver (∂T/∂n = 0) */
                    break;
                }

                case BC_PRESSURE_OUTLET:
                    field[ii].p[i] = myPointStruct[ii].node_bc[i].p;
                    myPointStruct[ii].flag_outlets = true;
                    break;

                case BC_SYMMETRY:
                    /* Zero normal velocity and gradients handled in solver */
                    break;

                case BC_SUPERSONIC_INLET:
                    if (parameters.compressible_flow) {
                        double vn = myPointStruct[ii].node_bc[i].v_n;
                        double vt = myPointStruct[ii].node_bc[i].v_t;

                        if (vn != 0.0 || vt != 0.0)
                        {
                            double nx = myPointStruct[ii].x_normal[i];
                            double ny = myPointStruct[ii].y_normal[i];
                            decompose_velocity(vn, vt, nx, ny,
                                               &field[ii].u[i], &field[ii].v[i]);
                            if (parameters.dimension == 3)
                                field[ii].w[i] = 0.0;
                        }
                        else
                        {
                            field[ii].u[i] = myPointStruct[ii].node_bc[i].u;
                            field[ii].v[i] = myPointStruct[ii].node_bc[i].v;
                            if (parameters.dimension == 3)
                                field[ii].w[i] = myPointStruct[ii].node_bc[i].w;
                        }

                        field[ii].p[i]   = myPointStruct[ii].node_bc[i].p;
                        field[ii].T[i]   = myPointStruct[ii].node_bc[i].T;
                        field[ii].rho[i] = myPointStruct[ii].node_bc[i].rho;
                    }
                    break;

                case BC_SUPERSONIC_OUTLET:
                    /* All variables extrapolated, handled in solver */
                    break;

                case BC_SUBSONIC_INLET:
                    if (parameters.compressible_flow) {
                        field[ii].T[i] = myPointStruct[ii].node_bc[i].T_total;
                        field[ii].p[i] = myPointStruct[ii].node_bc[i].p_total;
                    }
                    break;

                case BC_SUBSONIC_OUTLET:
                    if (parameters.compressible_flow) {
                        field[ii].p[i] = myPointStruct[ii].node_bc[i].p;
                    }
                    break;

                default:
                    field[ii].u[i] = 0.0;
                    field[ii].v[i] = 0.0;
                    if (parameters.dimension == 3)
                        field[ii].w[i] = 0.0;
                    field[ii].p[i] = 0.0;

                    if (parameters.compressible_flow) {
                        field[ii].T[i]   = parameters.T_ref;
                        field[ii].rho[i] = parameters.rho_ref;
                    }

                    printf("Applying default BC at node %d: u=v=w=p=0.0\n", i);
                    break;
            }
        }
    }
}

void check_restart_file(PointStructure* myPointStruct, FieldVariables *field) {
    if (parameters.restart) {
        printf("Checking restart file: %s\n", parameters.restart_filename);
        FILE *file = fopen(parameters.restart_filename, "r");
        if (file) {
            printf("Restart file found: %s\n", parameters.restart_filename);
            for (int i = 0; i < myPointStruct->num_nodes; i++) {
                double x, y, z, u, v, w, p;
                if (fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf", &x, &y, &z, &u, &v, &w, &p) != 7) {
                    fprintf(stderr, "Error reading restart file at line %d\n", i);
                    fclose(file);
                    exit(1);
                }
                field->u[i] = u;
                field->v[i] = v;
                field->w[i] = w;
                field->p[i] = p;
            }
            
        } 
        else {
            printf("Restart file not found: %s\n", parameters.restart_filename);
            exit(1);
        }
    }
}    
