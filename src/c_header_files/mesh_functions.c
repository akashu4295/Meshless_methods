// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "functions.h"
#include "kdtree.h"

// Safe guards for fscanf on reading 1 and 2 items
#define SAFE_SCAN1(x, msg) if ((x) != 1) { puts(msg); exit(1); }
#define SAFE_SCAN2(x, msg) if ((x) != 2) { puts(msg); exit(1); }
#define LINE_LEN 512

///////////////////////////////////////////////////////////////////////////////
// Function definitions
///////////////////////////////////////////////////////////////////////////////

static void trim(char *s)
{
    char *p = s;
    while (isspace(*p)) p++;
    memmove(s, p, strlen(p) + 1);

    for (int i = strlen(s) - 1; i >= 0 && isspace(s[i]); i--)
        s[i] = '\0';
}

// Remove // and block comments 
static void strip_comments(char *line, bool *in_block_comment)
{
    char *p = line;

    while (*p)
    {
        if (*in_block_comment)
        {
            if (p[0] == '*' && p[1] == '/')
            {
                *in_block_comment = false;
                memmove(p, p + 2, strlen(p + 2) + 1);
            }
            else
            {
                memmove(p, p + 1, strlen(p + 1) + 1);
            }
        }
        else
        {
            if (p[0] == '/' && p[1] == '*')
            {
                *in_block_comment = true;
                memmove(p, p + 2, strlen(p + 2) + 1);
            }
            else if (p[0] == '/' && p[1] == '/')
            {
                *p = '\0';
                break;
            }
            else
            {
                p++;
            }
        }
    }
}

static int parse_key_value(char *line, char *key, char *value)
{
    char *comma = strchr(line, ',');
    if (!comma) return 0;

    *comma = '\0';
    strcpy(key, line);
    strcpy(value, comma + 1);

    trim(key);
    trim(value);

    return 1;
}

void read_flow_parameters(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file){
        perror("Cannot open flow parameter file");
        exit(1);
    }

    char line[LINE_LEN], key[128], val[128];
    bool in_block_comment = false;

    while (fgets(line, sizeof(line), file)){
        strip_comments(line, &in_block_comment);
        trim(line);

        if (line[0] == '\0') continue;
        if (!parse_key_value(line, key, val)) continue;

        /* ---- dispatch table ---- */
        if (!strcmp(key, "domain_dimensions")){
            parameters.dimension = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.dimension);
        }
        else if (!strcmp(key, "poly_deg")){
            parameters.poly_degree = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.poly_degree);
        }
        else if (!strcmp(key, "phs_deg")){
            parameters.phs_degree = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.phs_degree);
        }
        else if (!strcmp(key, "cloud_size_multiplier")){
            parameters.cloud_size_multiplier = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.cloud_size_multiplier);
        }
        else if (!strcmp(key, "test_derivative")){
            parameters.test = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.test);
        }
        else if (!strcmp(key, "courant_number")){
            parameters.courant_number = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.courant_number);
        }
        else if (!strcmp(key, "steady_tolerance")){
            parameters.steady_state_tolerance = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.steady_state_tolerance);
        }
        else if (!strcmp(key, "poisson_solver_tolerance")){
            parameters.poisson_solver_tolerance = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.poisson_solver_tolerance);
        }
        else if (!strcmp(key, "sor_parameter")){
            parameters.omega = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.omega);
        }
        else if (!strcmp(key, "time_step")){
            parameters.dt = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.dt);
        }
        else if (!strcmp(key, "num_time_steps")){
            parameters.num_time_steps = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.num_time_steps);
        }
        else if (!strcmp(key, "write_interval")){
            parameters.write_interval = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.write_interval);
        }
        else if (!strcmp(key, "Re")){
            parameters.Re = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.Re);
        }
        else if (!strcmp(key, "iter_momentum")){
            parameters.iter_momentum = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.iter_momentum);
        }
        else if (!strcmp(key, "iter_timple")){
            parameters.iter_timple = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.iter_timple);
        }
        else if (!strcmp(key, "num_vcycles")){
            parameters.num_vcycles = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.num_vcycles);
        }
        else if (!strcmp(key, "num_relax")){
            parameters.num_relax = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.num_relax);
        }
        else if (!strcmp(key, "fractional_step")){
            parameters.fractional_step = atoi(val) != 0;
            printf("PARAMETERS: %s = %hd\n", key, parameters.fractional_step);
        }
        else if (!strcmp(key, "facRe")){
            parameters.facRe = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.facRe);
        }
        else if (!strcmp(key, "facdt")){
            parameters.facdt = atof(val);
            printf("PARAMETERS: %s = %f\n", key, parameters.facdt);
        }
        else if (!strcmp(key, "Poisson_solver_type")){
            parameters.poisson_solver_type = atoi(val);
            printf("PARAMETERS: %s = %hd\n", key, parameters.poisson_solver_type);
        }
        else if (!strcmp(key, "restart")){
            parameters.restart = atoi(val) != 0;
            printf("PARAMETERS: %s = %hd\n", key, parameters.restart);
        }
        else if (!strcmp(key, "restart_filename")){
            strncpy(parameters.restart_filename, val, 249);
            parameters.restart_filename[249] = '\0';
            printf("PARAMETERS: %s = %s\n", key, parameters.restart_filename);
        }
        /* ---- compressible block ---- */

        else if (!strcmp(key, "compressible_flow")){
            parameters.compressible_flow = atoi(val) != 0;
            printf("PARAMETERS: %s = %hd\n", key, parameters.compressible_flow);
        }
        else if (!strcmp(key, "gamma")){
            parameters.gamma = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.gamma);
        }
        else if (!strcmp(key, "R_gas")){
            parameters.R_gas = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.R_gas);
        }
        else if (!strcmp(key, "Pr")){
            parameters.Pr = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.Pr);
        }
        else if (!strcmp(key, "Cv")){
            parameters.cv = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.cv);
        }
        else if (!strcmp(key, "Cp")){
            parameters.cp = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.cp);
        }
        else if (!strcmp(key, "T_ref")){
            parameters.T_ref = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.T_ref);
        }
        else if (!strcmp(key, "rho_ref")){
            parameters.rho_ref = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.rho_ref);
        }
        else if (!strcmp(key, "p_ref")){
            parameters.p_ref = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.p_ref);
        }
        else if (!strcmp(key, "mu_ref")){
            parameters.mu_ref = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.mu_ref);
        }
        else if (!strcmp(key, "T_sutherland")){
            parameters.T_sutherland = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.T_sutherland);
        }
        else if (!strcmp(key, "viscosity_model")){
            parameters.viscosity_model = atoi(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %hd\n", key, parameters.viscosity_model);
        }
        else if (!strcmp(key, "Mach")){
            parameters.Mach = atof(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %f\n", key, parameters.Mach);
        }
        else if (!strcmp(key, "energy_equation")){
            parameters.energy_equation = atoi(val);
            if (parameters.compressible_flow)
                printf("PARAMETERS: %s = %hd\n", key, parameters.energy_equation);
        }
        else
            printf("!! Unknown parameter ignored: %s\n", key);
    }

    fclose(file);

    /* ---- derived quantities ---- */
    if (parameters.compressible_flow)
    {
        parameters.rho = parameters.rho_ref;
        parameters.mu  = parameters.mu_ref;
    }
    else
    {
        parameters.rho = 1.0;
        parameters.mu  = 1.0 / parameters.Re;
    }
    parameters.nu = parameters.mu / parameters.rho;
}

void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels)
{   
    FILE *file;
    file = fopen(filename, "r");
    char temp[250];
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    SAFE_SCAN2(fscanf(file, "%[^,],%hd\n", temp, num_levels), "Number of levels read error, Check grid_filenames.csv");
    printf("PARAMETERS: %s = %hd\n", temp, *num_levels);
    
    // Allocate memory for point structure for all levels and read meshfile names
    *myPointStruct = (PointStructure*)malloc((*num_levels) * sizeof(PointStructure));
    for (short ii = 0; ii<*num_levels ; ii = ii +1){
        SAFE_SCAN1(fscanf(file, "%s\n", (*myPointStruct)[ii].mesh_filename), "Mesh filename read error, Check grid_filenames.csv");
        printf("Mesh filename for level %d : %s\n", ii, (*myPointStruct)[ii].mesh_filename);
    }
    fclose(file);
}

void read_PointStructure(PointStructure* myPointStruct)
{   
    if (parameters.dimension == 3)
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2) * (myPointStruct->poly_degree + 3)/ 6;
    else
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2)/ 2;
    myPointStruct->num_cloud_points = parameters.cloud_size_multiplier * myPointStruct->num_poly_terms;
    printf("PARAMETERS: num_poly_terms = %d\n", myPointStruct->num_poly_terms);
    printf("PARAMETERS: num_cloud_points = %d\n", myPointStruct->num_cloud_points);
    printf("PARAMETERS: dimension = %d\n", parameters.dimension);

    FILE *file;
    char filename[250];
    double dtemp; int itemp; char temp[50]; // temporary variables used multiple times in the code
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
        SAFE_SCAN1(fscanf(file, "%s ", temp), "Error reading mesh file, Check mesh file format");
        if (strcmp(temp, "$Nodes") == 0)
            break;
    }
    SAFE_SCAN1(fscanf(file, "%i ", &itemp), "Number of nodes read error, Check mesh file format");

    // Allocate memory for the members of point structure
    AllocateMemoryPointStructure(myPointStruct, itemp);
    
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

    FILE *bcf = fopen("bc.csv", "r");
    if (bcf != NULL)
    {
        fclose(bcf);  // close immediately, we just tested existence
        read_physical_names(filename, myPointStruct); // Read physical names from the mesh file
        read_boundary_conditions_file("bc.csv", myPointStruct);
    }
    else
    {
        printf("bc.csv not found: using default / interior boundary conditions\n");
    }

    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {
        SAFE_SCAN1(fscanf(file, "%i ", &itemp), "Node number read error, Check mesh file format"); //node number
        SAFE_SCAN1(fscanf(file, "%lf ", &dtemp), "X coordinate read error, Check mesh file format");
        myPointStruct->x[i]=dtemp;
        SAFE_SCAN1(fscanf(file, "%lf ", &dtemp), "Y coordinate read error, Check mesh file format"); 
        myPointStruct->y[i]=dtemp;
        SAFE_SCAN1(fscanf(file, "%lf ", &dtemp), "Z coordinate read error, Check mesh file format"); 
        myPointStruct->z[i]=dtemp;
        
        myPointStruct->point_index[i] = i;
        myPointStruct->node_bc[i].type = BC_INTERIOR;
        myPointStruct->node_bc[i].u = 0;
        myPointStruct->node_bc[i].v = 0;
        myPointStruct->node_bc[i].w = 0;
        myPointStruct->node_bc[i].p = 0;
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
        SAFE_SCAN1(fscanf(file, "%s ", temp), "Error reading mesh file, Check mesh file format");
        if (strcmp(temp, "$Elements") == 0)
            break;
    }
    SAFE_SCAN1(fscanf(file, "%i ", &itemp), "Number of elements read error, Check mesh file format");
    myPointStruct->num_elem = itemp;

    if (parameters.dimension == 2) // Calculation done only for kdtree search
        myPointStruct->d_avg = sqrt((x_max - x_min)*(y_max - y_min) / myPointStruct->num_elem);
    else 
        myPointStruct->d_avg = cbrt((x_max - x_min)*(y_max - y_min)*(z_max - z_min)/myPointStruct->num_elem);

    // Read the mesh file elements and boundary nodes
    int e_type, e_node1, e_node2, e_node3, e_node4;
    int num_tags, geom_id;
    short phys_id;
    double dx, dy, dz, dx1, dy1, dz1;
    
    if (parameters.dimension == 2)
    {
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscanf(file, "%d %d", &itemp, &e_type);

            if (e_type == 1 || e_type == 15)   // line or point
            {
                fscanf(file, "%d", &num_tags);

                phys_id = 0;
                geom_id = 0;
                if (num_tags >= 1) fscanf(file, "%hd", &phys_id);
                if (num_tags >= 2) fscanf(file, "%d", &geom_id);
                for (int t = 2; t < num_tags; t++) fscanf(file, "%*d");
                
                BCValue bc;
                bc.type = BC_INTERIOR;
                bc.u = 0.0; bc.v = 0.0; bc.w = 0.0; bc.p = 0.0; 
                if (parameters.compressible_flow){
                    bc.p_total = parameters.p_ref;
                    bc.rho = parameters.rho_ref;
                    bc.T = parameters.T_ref;
                    bc.T_total = parameters.T_ref;
                }
                
                if ((phys_id != 0))
                {
                    for (int i = 0; i < myPointStruct->num_boundary_types; i++)
                    {
                        // printf("%d %d \n", myPointStruct->boundary_map[i].physical_id, phys_id);
                        if (myPointStruct->boundary_map[i].physical_id == phys_id)
                        {
                            bc = myPointStruct->boundary_map[i].bc;
                            break;
                        }
                    }
                }

                if (e_type == 1)  // 2-node line
                {
                    fscanf(file, "%d %d", &e_node1, &e_node2);

                    myPointStruct->boundary_tag[e_node1 - 1] = true;
                    myPointStruct->boundary_tag[e_node2 - 1] = true;

                    if (myPointStruct->node_bc[e_node1 - 1].type == BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node1 - 1, bc);
                    }
                    else if (myPointStruct->node_bc[e_node1 - 1].type != bc.type) {
                        myPointStruct->corner_tag[e_node1 - 1] = true;
                        assign_node_bc(myPointStruct, e_node1 - 1, bc);  // priority-based
                    }

                    if (myPointStruct->node_bc[e_node2 - 1].type == BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node2 - 1, bc);
                    }
                    else if (myPointStruct->node_bc[e_node2 - 1].type != bc.type) {
                        myPointStruct->corner_tag[e_node2 - 1] = true;
                        assign_node_bc(myPointStruct, e_node2 - 1, bc);  // priority-based
                    }

                    // printf(" %d %d \n", myPointStruct->node_bc[e_node1 - 1].type, myPointStruct->node_bc[e_node2 - 1].type);
                    dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1];
                    dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];

                    myPointStruct->x_normal[e_node1 - 1] += -dy;
                    myPointStruct->y_normal[e_node1 - 1] +=  dx;
                    myPointStruct->x_normal[e_node2 - 1] += -dy;
                    myPointStruct->y_normal[e_node2 - 1] +=  dx;
                }
                else   // e_type == 15 -> point (corner)
                {
                    fscanf(file, "%d", &itemp);
                    myPointStruct->corner_tag[itemp - 1]   = true;
                    myPointStruct->boundary_tag[itemp - 1] = true;

                    if (bc.type != BC_INTERIOR) {
                        assign_node_bc(myPointStruct, itemp - 1, bc);
                    }
                }
            }
            else
            {
                fscanf(file, "%*[^\n]\n");
            }
        }
    }
    else   // 3D 
    {
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscanf(file, "%d %d", &itemp, &e_type);

            if (e_type == 2 || e_type == 3 || e_type == 1 || e_type == 15)
            {
                fscanf(file, "%d", &num_tags);

                phys_id = 0;
                geom_id = 0;
                if (num_tags >= 1) fscanf(file, "%hd", &phys_id);
                if (num_tags >= 2) fscanf(file, "%d", &geom_id);
                for (int t = 2; t < num_tags; t++) fscanf(file, "%*d");

                BCValue bc;
                bc.type = BC_INTERIOR;
                bc.u = 0.0; bc.v = 0.0; bc.w = 0.0; bc.p = 0.0; 
                if (parameters.compressible_flow){
                    bc.p_total = parameters.p_ref;
                    bc.rho = parameters.rho_ref;
                    bc.T = parameters.T_ref;
                    bc.T_total = parameters.T_ref;
                }

                if ((phys_id != 0))
                {
                    for (int i = 0; i < myPointStruct->num_boundary_types; i++)
                    {
                        if (myPointStruct->boundary_map[i].physical_id == phys_id)
                        {
                            bc = myPointStruct->boundary_map[i].bc;
                            break;
                        }
                    }
                }

                if (e_type == 2 || e_type == 3)   // triangle or quad
                {
                    fscanf(file, "%d %d %d", &e_node1, &e_node2, &e_node3);

                    myPointStruct->boundary_tag[e_node1 - 1] = true;
                    myPointStruct->boundary_tag[e_node2 - 1] = true;
                    myPointStruct->boundary_tag[e_node3 - 1] = true;


                    if (myPointStruct->node_bc[e_node1 - 1].type == BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node1 - 1, bc);
                    }
                    else if (myPointStruct->node_bc[e_node1 - 1].type != bc.type) {
                        myPointStruct->corner_tag[e_node1 - 1] = true;
                        assign_node_bc(myPointStruct, e_node1 - 1, bc);
                    }
                    if (myPointStruct->node_bc[e_node2 - 1].type == BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node2 - 1, bc);
                    }
                    else if (myPointStruct->node_bc[e_node2 - 1].type != bc.type) {
                        myPointStruct->corner_tag[e_node2 - 1] = true;
                        assign_node_bc(myPointStruct, e_node2 - 1, bc);
                    }
                    if (myPointStruct->node_bc[e_node3 - 1].type == BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node3 - 1, bc);
                    }
                    else if (myPointStruct->node_bc[e_node3 - 1].type != bc.type) {
                        myPointStruct->corner_tag[e_node3 - 1] = true;
                        assign_node_bc(myPointStruct, e_node3 - 1, bc);
                    }

                    if (e_type == 3)
                    {
                        fscanf(file, "%d", &e_node4);
                        myPointStruct->boundary_tag[e_node4 - 1] = true;
                        if (myPointStruct->node_bc[e_node4 - 1].type == BC_INTERIOR) {
                            assign_node_bc(myPointStruct, e_node4 - 1, bc);
                        }
                        else if (myPointStruct->node_bc[e_node4 - 1].type != bc.type) {
                            myPointStruct->corner_tag[e_node4 - 1] = true;
                            assign_node_bc(myPointStruct, e_node4 - 1, bc);
                        }
                    }
                    if (e_type == 2)  // triangle
                    {
                        dx  = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1];
                        dy  = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                        dz  = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
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
                    else  // e_type == 3 (quad) - split into two triangles
                    {
                        // First triangle: nodes 1-2-3
                        dx  = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1];
                        dy  = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                        dz  = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
                        dx1 = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                        dy1 = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                        dz1 = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];

                        double nx1 = dy*dz1 - dz*dy1;
                        double ny1 = dz*dx1 - dx*dz1;
                        double nz1 = dx*dy1 - dy*dx1;

                        // Second triangle: nodes 1-3-4
                        dx  = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                        dy  = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                        dz  = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];
                        dx1 = myPointStruct->x[e_node4 - 1] - myPointStruct->x[e_node1 - 1];
                        dy1 = myPointStruct->y[e_node4 - 1] - myPointStruct->y[e_node1 - 1];
                        dz1 = myPointStruct->z[e_node4 - 1] - myPointStruct->z[e_node1 - 1];

                        double nx2 = dy*dz1 - dz*dy1;
                        double ny2 = dz*dx1 - dx*dz1;
                        double nz2 = dx*dy1 - dy*dx1;

                        // Average and accumulate
                        myPointStruct->x_normal[e_node1 - 1] += (nx1 + nx2) * 0.5;
                        myPointStruct->x_normal[e_node2 - 1] += nx1 * 0.5;
                        myPointStruct->x_normal[e_node3 - 1] += (nx1 + nx2) * 0.5;
                        myPointStruct->x_normal[e_node4 - 1] += nx2 * 0.5;

                        myPointStruct->y_normal[e_node1 - 1] += (ny1 + ny2) * 0.5;
                        myPointStruct->y_normal[e_node2 - 1] += ny1 * 0.5;
                        myPointStruct->y_normal[e_node3 - 1] += (ny1 + ny2) * 0.5;
                        myPointStruct->y_normal[e_node4 - 1] += ny2 * 0.5;

                        myPointStruct->z_normal[e_node1 - 1] += (nz1 + nz2) * 0.5;
                        myPointStruct->z_normal[e_node2 - 1] += nz1 * 0.5;
                        myPointStruct->z_normal[e_node3 - 1] += (nz1 + nz2) * 0.5;
                        myPointStruct->z_normal[e_node4 - 1] += nz2 * 0.5;
                    }
                }
                else if (e_type == 1)   // line -> corners
                {
                    fscanf(file, "%d %d", &e_node1, &e_node2);
                    myPointStruct->corner_tag[e_node1 - 1] = true;
                    myPointStruct->corner_tag[e_node2 - 1] = true;
                    myPointStruct->boundary_tag[e_node1 - 1] = true;    
                    myPointStruct->boundary_tag[e_node2 - 1] = true;
                    if (bc.type != BC_INTERIOR) {
                        assign_node_bc(myPointStruct, e_node1 - 1,  bc);
                        assign_node_bc(myPointStruct, e_node2 - 1,  bc);
                    }
                }
                else   // e_type == 15
                {
                    fscanf(file, "%d", &itemp);
                    myPointStruct->corner_tag[itemp - 1] = true;
                    myPointStruct->boundary_tag[itemp - 1] = true;

                    if (bc.type != BC_INTERIOR) {
                        assign_node_bc(myPointStruct, itemp - 1, bc);
                    }
                }
            }
            else
            {
                fscanf(file, "%*[^\n]\n");
            }
        }
    }

    fclose(file);
    myPointStruct->num_corners = 0;
    myPointStruct->num_boundary_nodes = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i] == true)
            myPointStruct->num_boundary_nodes++;
        if (myPointStruct->corner_tag[i] == true)
            myPointStruct->num_corners++;
    }
    // Remove first num_corners nodes
    int remove_count = myPointStruct->num_corners;
    int new_size = myPointStruct->num_nodes - myPointStruct->num_corners;

    if ((remove_count != 0) && (myPointStruct->num_boundary_types == 0)) {
        printf("Removing %d corner nodes from the mesh\n", remove_count);
        memmove(myPointStruct->x, myPointStruct->x + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->y, myPointStruct->y + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->z, myPointStruct->z + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->x_normal, myPointStruct->x_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->y_normal, myPointStruct->y_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->z_normal, myPointStruct->z_normal + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(double));
        memmove(myPointStruct->point_index, myPointStruct->point_index + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(int));
        memmove(myPointStruct->boundary_tag, myPointStruct->boundary_tag + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(bool));
        memmove(myPointStruct->corner_tag, myPointStruct->corner_tag + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(bool));
        memmove(myPointStruct->node_bc, myPointStruct->node_bc + remove_count, (myPointStruct->num_nodes - remove_count) * sizeof(BCValue));

        // Reallocate memory to shrink the array
        myPointStruct->x = (double*)realloc(myPointStruct->x, new_size * sizeof(double));
        myPointStruct->y = (double*)realloc(myPointStruct->y, new_size * sizeof(double));
        myPointStruct->z = (double*)realloc(myPointStruct->z, new_size * sizeof(double));
        myPointStruct->x_normal = (double*)realloc(myPointStruct->x_normal, new_size * sizeof(double));
        myPointStruct->y_normal = (double*)realloc(myPointStruct->y_normal, new_size * sizeof(double));
        myPointStruct->z_normal = (double*)realloc(myPointStruct->z_normal, new_size * sizeof(double));
        myPointStruct->point_index = (int*)realloc(myPointStruct->point_index, new_size * sizeof(int));
        myPointStruct->boundary_tag = (bool*)realloc(myPointStruct->boundary_tag, new_size * sizeof(bool));
        myPointStruct->corner_tag = (bool*)realloc(myPointStruct->corner_tag, new_size * sizeof(bool));
        myPointStruct->node_bc = (BCValue*)realloc(myPointStruct->node_bc, new_size * sizeof(BCValue));

        if (myPointStruct->x == NULL || myPointStruct->y == NULL || myPointStruct->z == NULL 
                || myPointStruct->x_normal == NULL || myPointStruct->y_normal == NULL || 
                myPointStruct->z_normal == NULL || myPointStruct->point_index == NULL || 
                myPointStruct->boundary_tag == NULL|| myPointStruct->corner_tag == NULL) {
            perror("realloc failed while deleting corner nodes");
            exit(1);
        }
        myPointStruct->num_nodes = new_size;
        myPointStruct->num_boundary_nodes = myPointStruct->num_boundary_nodes - myPointStruct->num_corners;
        // myPointStruct->num_corners = 0;

        for(int i = 0; i<new_size; i++){
            myPointStruct->point_index[i] = myPointStruct->point_index[i]-remove_count;
        }
    }
    else if ((myPointStruct->num_corners != 0) ) {

        // Allocate temporary arrays for non-corner nodes
        double *temp_x = (double*)malloc(new_size * sizeof(double));
        double *temp_y = (double*)malloc(new_size * sizeof(double));
        double *temp_z = (double*)malloc(new_size * sizeof(double));
        double *temp_x_normal = (double*)malloc(new_size * sizeof(double));
        double *temp_y_normal = (double*)malloc(new_size * sizeof(double));
        double *temp_z_normal = (double*)malloc(new_size * sizeof(double));
        int *temp_point_index = (int*)malloc(new_size * sizeof(int));
        bool *temp_boundary_tag = (bool*)malloc(new_size * sizeof(bool));
        bool *temp_corner_tag = (bool*)malloc(new_size * sizeof(bool));
        BCValue *temp_node_bc = (BCValue*)malloc(new_size * sizeof(BCValue));

        if (!temp_x || !temp_y || !temp_z || !temp_x_normal || !temp_y_normal || 
            !temp_z_normal || !temp_point_index || !temp_boundary_tag || 
            !temp_corner_tag || !temp_node_bc) {
            perror("malloc failed during corner node removal");
            // Free any successfully allocated memory
            free(temp_x); free(temp_y); free(temp_z);
            free(temp_x_normal); free(temp_y_normal); free(temp_z_normal);
            free(temp_point_index); free(temp_boundary_tag); 
            free(temp_corner_tag); free(temp_node_bc);
            exit(1);
        }

        // Copy non-corner nodes to temporary arrays
        int j = 0;
        for (int i = 0; i < myPointStruct->num_nodes; i++) {
            if (!myPointStruct->corner_tag[i]) {
                temp_x[j] = myPointStruct->x[i];
                temp_y[j] = myPointStruct->y[i];
                temp_z[j] = myPointStruct->z[i];
                temp_x_normal[j] = myPointStruct->x_normal[i];
                temp_y_normal[j] = myPointStruct->y_normal[i];
                temp_z_normal[j] = myPointStruct->z_normal[i];
                temp_point_index[j] = myPointStruct->point_index[i];
                temp_boundary_tag[j] = myPointStruct->boundary_tag[i];
                temp_corner_tag[j] = false; // All corners removed
                temp_node_bc[j] = myPointStruct->node_bc[i];
                j++;
            }
        }

        // Free old arrays and assign new ones
        free(myPointStruct->x);
        free(myPointStruct->y);
        free(myPointStruct->z);
        free(myPointStruct->x_normal);
        free(myPointStruct->y_normal);
        free(myPointStruct->z_normal);
        free(myPointStruct->point_index);
        free(myPointStruct->boundary_tag);
        free(myPointStruct->corner_tag);
        free(myPointStruct->node_bc);

        myPointStruct->x = temp_x;
        myPointStruct->y = temp_y;
        myPointStruct->z = temp_z;
        myPointStruct->x_normal = temp_x_normal;
        myPointStruct->y_normal = temp_y_normal;
        myPointStruct->z_normal = temp_z_normal;
        myPointStruct->point_index = temp_point_index;
        myPointStruct->boundary_tag = temp_boundary_tag;
        myPointStruct->corner_tag = temp_corner_tag;
        myPointStruct->node_bc = temp_node_bc;
        myPointStruct->num_boundary_nodes -= myPointStruct->num_corners;
        myPointStruct->num_nodes = new_size;
        for (int i = 0; i < new_size; i++) {
            myPointStruct->point_index[i] = i;
        }
    }
    for(int i = 0; i<new_size; i++){
        myPointStruct->point_index[i] = i;
    }
    printf("No of nodes = %d \nNo of elements = %d \n", myPointStruct->num_nodes, myPointStruct->num_elem);
}

void correct_normal_directions(PointStructure *myPointStruct)
{
    double dx, dy, dz, mag, dot_product, x_avg, y_avg, z_avg;
    short count = 0;
    int n = myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {   
        if(myPointStruct->boundary_tag[i] == true)
        {   
            x_avg = 0; y_avg = 0; z_avg = 0; count = 0;

            for (int j = myPointStruct->num_cloud_points-1; j > 0; j--){
                if (myPointStruct->boundary_tag[myPointStruct->cloud_index[i*n +j]] == false){
                    x_avg += myPointStruct->x[myPointStruct->cloud_index[i*n +j]];
                    y_avg += myPointStruct->y[myPointStruct->cloud_index[i*n +j]];
                    z_avg += myPointStruct->z[myPointStruct->cloud_index[i*n +j]];
                    count++;
                    if (count == 1)
                        break;
                }
            }
            //x_avg = x_avg/count; y_avg = y_avg/count; z_avg = z_avg/count;

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
            if (mag < 1e-14) continue; // or set a default normal
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

    myPointStruct_c->restr_mat = (double*)malloc(myPointStruct_c->num_nodes * m * sizeof(double));

    // for (int i = 0; i < myPointStruct_c->num_nodes; i++)
    //     myPointStruct_c->restr_mat[i] = (double*)malloc(m * sizeof(double));

    // Initialise to zeros
    for (int i = 0; i < myPointStruct_c->num_nodes * m; i++)
        myPointStruct_c->restr_mat[i] = 0;
    
    double *A_inv = create_matrix_vectorised(mpn,mpn);
    double *A = create_matrix_vectorised(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m); // coefficients of restriction matrix
    int i_restr;

    double point1[3], point2[3];
    for (int i = myPointStruct_c->num_boundary_nodes; i < myPointStruct_c->num_nodes; i++)
    {   
        i_restr = myPointStruct_c->restriction_points[i];
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct_f, A, myPointStruct_f->cloud_index[i_restr]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m+n);
        
        point1[0] = myPointStruct_c->x[i];
        point1[1] = myPointStruct_c->y[i];
        point1[2] = myPointStruct_c->z[i];
        
        for (short j = 0; j < m; j++) {
            int k = (i_restr - 1)* myPointStruct_f->num_cloud_points + j;
            point2[0] = myPointStruct_f->x[myPointStruct_f->cloud_index[k]];
            point2[1] = myPointStruct_f->y[myPointStruct_f->cloud_index[k]];
            point2[2] = myPointStruct_f->z[myPointStruct_f->cloud_index[k]];
            temp[j] = calculate_phs_rbf(point1, point2, parameters.phs_degree, parameters.dimension);
        }
        double seed_pt[3] = {myPointStruct_f->x[i_restr],
                             myPointStruct_f->y[i_restr],
                             myPointStruct_f->z[i_restr]};
        for (short j = 0; j < n; j++) 
            temp[j+m] = pow(point1[0]-seed_pt[0], myPointStruct_f->pow_x[j])
                                * pow(point1[1]-seed_pt[1], myPointStruct_f->pow_y[j])
                                * pow(point1[2]-seed_pt[2], myPointStruct_f->pow_z[j]);
        
        multiply_vector_matrix_columnwise_vectorised(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < m; j++)
        {
            myPointStruct_c->restr_mat[i*m +j] = temp1[j];
        }
    }
    free(temp);
    free(temp1);
    free(A);
    free(A_inv);
}

void create_prolongation_matrix(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c)
{
    // allocate memory for restriction matrix
    
    short m = myPointStruct_c->num_cloud_points;
    short n = myPointStruct_c->num_poly_terms;
    short mpn = m+n;
    
    myPointStruct_f->prol_mat = (double*)malloc(myPointStruct_f->num_nodes *m * sizeof(double));
    
    // Initialise to zeros
    for (int i = 0; i < myPointStruct_f->num_nodes * m; i++)
    {
        myPointStruct_f->prol_mat[i] = 0;
    }
    
    double *A_inv = create_matrix_vectorised(mpn,mpn);
    double *A = create_matrix_vectorised(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m);
    int i_prol;

    double point1[3], point2[3], seed_pt[3];
    for (int i = myPointStruct_f->num_boundary_nodes;
                     i < myPointStruct_f->num_nodes; i++)
    {   
        i_prol = myPointStruct_f->prolongation_points[i];
        create_A_matrix_from_cloud_indices_vectorised(myPointStruct_c, A, myPointStruct_c->cloud_index[i_prol]);
        matrixInverse_Gauss_Jordan_vectorised(A, A_inv, m+n);
        
        point1[0] = myPointStruct_f->x[i];
        point1[1] = myPointStruct_f->y[i];
        point1[2] = myPointStruct_f->z[i];
        
        for (short j = 0; j < m; j++) {
            int k = (i_prol-1)*myPointStruct_c->num_cloud_points + j;
            point2[0] = myPointStruct_c->x[myPointStruct_c->cloud_index[k]];
            point2[1] = myPointStruct_c->y[myPointStruct_c->cloud_index[k]];
            point2[2] = myPointStruct_c->z[myPointStruct_c->cloud_index[k]];
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
        
        multiply_vector_matrix_columnwise_vectorised(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < (m); j++)
            myPointStruct_f->prol_mat[i*m +j] = temp1[j];
    }
    free(temp);
    free(temp1);
    free(A);
    free(A_inv);
}

void rcm_reordering_with_boundarynodes(PointStructure* myPointstruct) {
    int m = myPointstruct->num_nodes;
    int n = myPointstruct->num_cloud_points;

    int *queue1 = (int*)malloc(m * sizeof(int));  // RCM ordered indices
    int *queue2 = (int*)malloc(m * sizeof(int));  // Maps from original to RCM
    bool *visited = (bool*)malloc(m * sizeof(bool));
    // Temporary arrays to hold reordered data
    double *temp_x = (double*)malloc(m * sizeof(double));
    double *temp_y = (double*)malloc(m * sizeof(double));
    double *temp_z = (double*)malloc(m * sizeof(double));
    double *temp_nx = (double*)malloc(m * sizeof(double));
    double *temp_ny = (double*)malloc(m * sizeof(double));
    double *temp_nz = (double*)malloc(m * sizeof(double));
    bool *temp_boundary_tag = (bool*)malloc(m * sizeof(bool));
    bool *temp_corner_tag = (bool*)malloc(m * sizeof(bool));
    BCValue *temp_node_bc = (BCValue*)malloc(m * sizeof(BCValue));
    int *temp_cloud_index = (int*)malloc(m*n * sizeof(int));

    // Initialize visited array and queue counters
    for (int i = 0; i < m; i++) {
        visited[i] = false;
    }

    int count = 0;        // Counter for filling queue1 (RCM order)
    int count_queue = 0;   // Pointer to track current node in BFS

    queue1[0] = 0;
    visited[0] = true;
    count = 1;
    count_queue = 0;

    while (count_queue < count) {
        int current = queue1[count_queue++];

        for (int i = 1; i < n; i++) {
            int neighbor = myPointstruct->cloud_index[current*n + i];

            if (neighbor < 0 || neighbor >= m)
                continue;

            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue1[count++] = neighbor;
            }
        }
    }

    /* Catch disconnected components */
    for (int i = 0; i < m; i++) {
        if (!visited[i]) {
            queue1[count++] = i;
            visited[i] = true;
        }
    }

    // Now queue1 contains nodes in the RCM order. Generate queue2 for reverse mapping
    for (int i = 0; i < m; i++) {
        queue2[queue1[i]] = i;  // Mapping original node to RCM index
    }

    // Reorder the data based on RCM order (queue1) and apply the reverse mapping (queue2)
    for (int i = 0; i < m; i++) {
        temp_x[i] = myPointstruct->x[queue1[i]];
        temp_y[i] = myPointstruct->y[queue1[i]];
        temp_z[i] = myPointstruct->z[queue1[i]];
        // Reorder the cloud index array based on queue2 mapping
        // temp_cloud_index[i*n] = queue1[i];
        for (int j = 0; j < n; j++) {
            temp_cloud_index[i*n +j] = queue2[myPointstruct->cloud_index[queue1[i]*n +j]];
        }
        temp_nx[i] = myPointstruct->x_normal[queue1[i]];
        temp_ny[i] = myPointstruct->y_normal[queue1[i]];
        temp_nz[i] = myPointstruct->z_normal[queue1[i]];
        temp_boundary_tag[i] = myPointstruct->boundary_tag[queue1[i]];
        temp_corner_tag[i] = myPointstruct->corner_tag[queue1[i]];
        temp_node_bc[i] = myPointstruct->node_bc[queue1[i]];
    }

    for (int i = 0; i < m; i++) {
        myPointstruct->x[i] = temp_x[i];
        myPointstruct->y[i] = temp_y[i];
        myPointstruct->z[i] = temp_z[i];
        myPointstruct->rcm_order[i] = queue2[i];
        for (int j = 0; j < n; j++) {
            myPointstruct->cloud_index[i*n +j] = temp_cloud_index[i*n+j];
        }

        myPointstruct->x_normal[i] = temp_nx[i];
        myPointstruct->y_normal[i] = temp_ny[i];
        myPointstruct->z_normal[i] = temp_nz[i];
        myPointstruct->boundary_tag[i] = temp_boundary_tag[i];
        myPointstruct->corner_tag[i] = temp_corner_tag[i];
        myPointstruct->node_bc[i] = temp_node_bc[i];
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
    free(temp_corner_tag);
    free(temp_cloud_index);
    free(temp_node_bc);
}


void create_prolongation_and_restriction_matrices(PointStructure* myPointStruct, short num_levels){  
    printf("Identifying prolongation and restriction points\n");
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

void calculate_point_spacings(PointStructure* myPointStruct){
    double dx, dy, dz;
    myPointStruct->d_avg = 0;
    myPointStruct->d_min = 1e10;
    myPointStruct->d_max = 0;
    double temp;
    int n = myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            dx = myPointStruct->x[i] - myPointStruct->x[myPointStruct->cloud_index[i*n +1]];
            dy = myPointStruct->y[i] - myPointStruct->y[myPointStruct->cloud_index[i*n +1]];
            if (parameters.dimension == 3)
                dz = myPointStruct->z[i] - myPointStruct->z[myPointStruct->cloud_index[i*n +1]];
            else
                dz = 0;
            temp = sqrt(dx*dx + dy*dy + dz*dz);
            myPointStruct->d_avg += temp;
            if (temp < myPointStruct->d_min)
                myPointStruct->d_min = temp;
            if (temp > myPointStruct->d_max)
                myPointStruct->d_max = temp;
        }
    }
    myPointStruct->d_avg = myPointStruct->d_avg/myPointStruct->num_nodes;
}

double calculate_dt(PointStructure* myPointStruct){
    double dt = 1e10;
    double d = myPointStruct->d_min;
    double termd =  parameters.dimension* parameters.nu/(d*d);
    double termc =    1/d;
    dt = 1/(termd + termc);
    return (fmin(dt*parameters.courant_number/parameters.dimension, parameters.dt));
}


///////////////////////////////////////////////////////////////////////////////

void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels)
{   
    for (short ii = 0; ii<num_levels ; ii = ii +1){
        myPointStruct[ii].poly_degree = 3; 
    	myPointStruct[0].poly_degree = parameters.poly_degree;
        printf("Reading mesh data for level %d\n", ii+1);
        read_PointStructure(&myPointStruct[ii]);
        printf("Identifying cloud indices\n");
        find_cloud_index(&myPointStruct[ii]);
        printf("RCM reordering of nodes\n");
        rcm_reordering_with_boundarynodes(&myPointStruct[ii]);
        printf("Correcting normal directions\n");
        correct_normal_directions(&myPointStruct[ii]);
        printf("Calculating nodal distances\n");
        calculate_point_spacings(&myPointStruct[ii]);
    }
    parameters.dt = calculate_dt(&myPointStruct[0]);
    create_prolongation_and_restriction_matrices(myPointStruct, num_levels);
}