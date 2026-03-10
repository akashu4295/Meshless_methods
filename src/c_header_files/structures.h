// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////////
// Structures
//////////////////////////////////////////////////////////////////////////////
#define MAX_BC_NAME 64
#define MAX_BC_TYPES 32

typedef enum {
    BC_INTERIOR = 0,
    BC_WALL = 1,
    BC_VELOCITY_INLET = 2,
    BC_PRESSURE_OUTLET = 3,
    BC_VELOCITY_OUTLET = 4,
    BC_ISOTHERMAL_WALL = 5,     // Wall with fixed T
    BC_ADIABATIC_WALL = 6,      // Wall with ∂T/∂n = 0
    BC_SYMMETRY = 7,            // Zero normal velocity and gradients
    BC_SUPERSONIC_INLET = 8,    // All variables specified
    BC_SUPERSONIC_OUTLET = 9,   // All variables extrapolated
    BC_SUBSONIC_INLET = 10,      // Specify T, p_total, direction
    BC_SUBSONIC_OUTLET = 11,     // Specify p_static, extrapolate others
    BC_DEFAULT = 20    // 
} BCType;

typedef struct {
    BCType type;
    double u, v, w;
    double v_n, v_t;
    double p, T, rho;
    double p_total, T_total;
} BCValue;

typedef struct {
    short physical_id;
    char  name[MAX_BC_NAME];
    BCValue bc;
} BoundaryMapEntry;

// Structure to represent the common parameters for all grids
struct parameters
    {
    short dimension; //dimension of the problem, default is 3
    short poly_degree; //degree of the polynomial basis functions
    short phs_degree; //degree of the PHS basis functions
    short cloud_size_multiplier; //multiplier for the number of cloud points
    short test;
    double courant_number; // Courant number for the time step
    double steady_state_tolerance; // Tolerance for steady state convergence
    double poisson_solver_tolerance; // Tolerance for Poisson solver
    short num_vcycles; // Number of V-cycles
    short num_relax; // Number of relaxation steps
    int num_time_steps; // Number of time steps
    short num_levels; // Number of levels in the multigrid
    int write_interval; // Interval for writing the data
    float omega; // relaxation parameter
    double dt; // time step
    short iter_momentum; // number of iterations for momentum equation
    short iter_timple; // number of iterations for timple time stepping
    double rho; // density
    double mu; // dynamic viscosity
    double Re; // Reynolds number
    double facRe;//Reynolds number factor for defect correction
    double facdt;//dual time step factor
    double nu; // kinematic viscosity
    bool fractional_step; // fractional step flag
    short poisson_solver_type; // Solver type Jacobi/Gauss Seidel/Bicgstab etc...
    bool restart; // flag to indicate whether to restart from a previous solution
    char restart_filename[250]; // filename to restart from
    bool compressible_flow;
    double gamma;          // Ratio of specific heats (1.4 for air)
    double R_gas;          // Gas constant (287 J/(kg·K) for air)
    double Pr;             // Prandtl number (0.71 for air)
    double cv;             // Specific heat at constant volume
    double cp;             // Specific heat at constant pressure (cp = gamma*cv/(gamma-1))
    double T_ref;          // Reference temperature (e.g., 300 K)
    double rho_ref;        // Reference density (e.g., 1.225 kg/m³)
    double p_ref;          // Reference pressure (e.g., 101325 Pa)
    double mu_ref;         // Reference viscosity (e.g., 1.81e-5 Pa·s)
    double T_sutherland;   // Sutherland's constant (110.4 K for air)
    int viscosity_model;   // 0=constant, 1=Sutherland, 2=power-law
    double Mach;           // Reference Mach number
    int energy_equation;   // 0=isothermal, 1=solve energy equation
    }; // child created only once and used globally

extern struct parameters parameters;

// Structure to represent the mesh data
typedef struct PointStructure {   
    char mesh_filename[250]; // name of the mesh file
    int num_nodes; // number of nodes
    int num_corners; // number of corner nodes
    int num_boundary_nodes; // number of boundary nodes
    int num_elem; // number of elements
    double d_avg; // average distance between nodes
    short num_poly_terms; //number of polynomial terms  //////NUM_POLY_TERMS
    short num_cloud_points; //number of cloud points in the domain
    short poly_degree; //degree of the polynomial basis functions
    BoundaryMapEntry boundary_map[MAX_BC_TYPES];
    short num_boundary_types;
    bool flag_outlets; // flag to indicate presence of outlets
    // BCType* node_bc_type; // size = num_nodes
    BCValue* node_bc;   // size = num_nodes
    double* x;  // x coordinates of the nodes
    double* y;  // y coordinates of the nodes
    double* z;  // z coordinates of the nodes
    int* rcm_order; // to get original order of data after rcm_reordering
    int* point_index; // index of the point in the original dataset
    double* x_normal; // x component of the normal vector
    double* y_normal; // y component of the normal vector
    double* z_normal; // z component of the normal vector
    bool* boundary_tag; // Boolean tag for boundary nodes
    bool* corner_tag; // Boolean tag for corner nodes
    short* pow_x; //power of x in the polynomial basis functions
    short* pow_y; //power of y in the polynomial basis functions
    short* pow_z; //power of z in the polynomial basis functions
    int* cloud_index; // indices of the INTERPOLATION cloud points 
    int* prolongation_points; // indices of the PROLONGATION cloud points
    int* restriction_points; // indices of the RESTRICTION cloud points
    double* Dx; // x-derivative matrix
    double* Dy; // y-derivative matrix
    double* Dz; // z-derivative matrix
    double* lap; // laplacian matrix
    double* lap_Poison; // laplacian matrix for Poisson solver
    double* restr_mat;// restriction matrix
    double* prol_mat;// prolongation matrix
}PointStructure;

typedef struct FieldVariables {
    double* u; // x-component of the velocity field
    double* v; // y-component of the velocity field
    double* w; // z-component of the velocity field
    double* p; // pressure
    double* u_new; // x-component of the intermediate velocity field
    double* v_new; // y-component of the intermediate velocity field
    double* w_new; // z-component of the intermediate velocity field
    double* u_old; // x-component of the old velocity field
    double* v_old; // y-component of the old velocity field
    double* w_old; // z-component of the old velocity field
    double* p_old; // pressure field
    double* pprime; // pressure field
    double* res; // residual field
    double* source; // source field
    double* dpdn; // dpdn field
    double* dpdx; // dpdx field
    double* dpdy; // dpdy field
    double* dpdz; // dpdz field
    double* dudx; // du/dx field
    double* dudy; // du/dy field
    double* dudz; // du/dz field
    double* dvdx; // dv/dx field
    double* dvdy; // dv/dy field
    double* dvdz; // dv/dz field
    double* dwdx; // dw/dx field
    double* dwdy; // dw/dy field
    double* dwdz; // dw/dz field
    double* lapu; // laplacian of u
    double* lapv; // laplacian of v
    double* lapw; // laplacian of w
    // Compressible flow variables
    double* rho; // density field
    double* rho_old; // density old field
    double* rho_new; // density new field
    double* T; // temperature field
    double* T_old; // Temperature old field
    double* T_new; // Temperature new field
    double* e; // Energy
    double* e_old; // Energy old
    double *drhodx, *drhody, *drhodz;
    double *dTdx, *dTdy, *dTdz;
    double *dedx, *dedy, *dedz;
    double *mu;
    double *kappa;
    double *tau_xx, *tau_yy, *tau_zz;
    double *tau_xy, *tau_xz, *tau_yz;
    double *div_tau_x, *div_tau_y, *div_tau_z;
    double *Q_visc;        // Viscous dissipation (heating)
    double *Q_source;      // External heat source term
} FieldVariables;

// All the functions, Assumes both init.c file and structures_vectors.c are linked to this header
void AllocateMemoryPointStructure(PointStructure* myPointStruct, int nodes);
void AllocateMemoryFieldVariables(FieldVariables** field, PointStructure* myPointStruct, int num_levels);
void free_PointStructure(PointStructure* myPointStruct, int num_levels);
void free_field(FieldVariables* field, int num_levels);
void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);
void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);


#endif
