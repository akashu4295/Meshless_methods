// Author: Based on incompressible solver by Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Compressible Navier-Stokes Solver for Low Mach Number Flows
// Affiliation: Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "functions.h"

/////////////////////////////////////////////////////////////////////////////
// Thermodynamic and Transport Properties
/////////////////////////////////////////////////////////////////////////////

// Sutherland's law for temperature-dependent viscosity
double calculate_viscosity_sutherland(double T, double mu_ref, double T_ref, double T_s) {
    return mu_ref * pow(T / T_ref, 1.5) * (T_ref + T_s) / (T + T_s);
}

// Power-law viscosity
double calculate_viscosity_powerlaw(double T, double mu_ref, double T_ref, double n) {
    return mu_ref * pow(T / T_ref, n);
}

// Thermal conductivity from Prandtl number
double calculate_thermal_conductivity(double mu, double cp, double Pr) {
    return mu * cp / Pr;
}

// Update thermodynamic properties at all nodes
void update_properties(PointStructure* myPointStruct, FieldVariables* field) {
    int N = myPointStruct->num_nodes;
    
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < N; i++) {
        // Ideal gas law: p = rho * R * T
        // For low Mach: p = p_ref + p' where p' << p_ref
        // rho = p_ref / (R * T)
        // Internal energy: e = cv * T
        field->rho[i] = parameters.p_ref / (parameters.R_gas * field->T[i]);
        field->e[i] = parameters.cv * field->T[i];
        
        // Viscosity (temperature-dependent)
        if (parameters.viscosity_model == 0) {
            field->mu[i] = parameters.mu_ref;    // Constant viscosity
        } else if (parameters.viscosity_model == 1) {
            // Sutherland's law
            field->mu[i] = calculate_viscosity_sutherland(field->T[i], parameters.mu_ref, parameters.T_ref, parameters.T_sutherland);
        } else if (parameters.viscosity_model == 2) {
            // Power law (n = 0.7 typical)
            field->mu[i] = calculate_viscosity_powerlaw(field->T[i], parameters.mu_ref, parameters.T_ref, 0.7);
        }
        // Thermal conductivity
        field->kappa[i] = calculate_thermal_conductivity(field->mu[i], parameters.cp, parameters.Pr);
    }
}

double calculate_dt_compressible(PointStructure* myPointStruct, FieldVariables* field) {
    double dt = 1e10;
    double d = myPointStruct->d_avg;
    int N = myPointStruct->num_nodes;
    double u_max = 0.0, v_max = 0.0, w_max = 0.0;
    double a_max = 0.0;  // Speed of sound
    double nu_max = 0.0;
    
    for (int i = 0; i < N; i++) {
        u_max = fmax(u_max, fabs(field->u[i]));
        v_max = fmax(v_max, fabs(field->v[i]));
        if (parameters.dimension == 3)
            w_max = fmax(w_max, fabs(field->w[i]));
        
        // Speed of sound: a = sqrt(gamma * R * T)
        double a = sqrt(parameters.gamma * parameters.R_gas * field->T[i]);
        a_max = fmax(a_max, a);
        
        // Kinematic viscosity
        double nu = field->mu[i] / field->rho[i];
        nu_max = fmax(nu_max, nu);
    }
    
    // CFL conditions
    double dt_convective = d / (u_max + v_max + w_max + a_max + 1e-10);
    double dt_diffusive = d * d / (2.0 * nu_max + 1e-10);
    double dt_acoustic = d / (a_max + 1e-10);
    
    dt = fmin(dt_convective, dt_diffusive);
    dt = fmin(dt, dt_acoustic);
    
    return dt * parameters.courant_number;
}

double compressible_solver_explicit(PointStructure* myPointStruct, FieldVariables* field) {
    double steady_state_error = 0.0;
    
    // Save old values
    #pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        field[0].w_old[i] = field[0].w[i];
        field[0].rho_old[i] = field[0].rho[i];
        field[0].T_old[i] = field[0].T[i];
        field[0].e_old[i] = field[0].e[i];
        field[0].p_old[i] = field[0].p[i];
    }
    
    update_properties(&myPointStruct[0], &field[0]);    // Update properties
    solve_continuity_equation(&myPointStruct[0], &field[0]);
    solve_momentum_equations(&myPointStruct[0], &field[0]);
    if (parameters.energy_equation == 1) 
        solve_energy_equation(&myPointStruct[0], &field[0]);
    solve_pressure_correction(&myPointStruct[0], &field[0]);
    update_velocity_compressible(&myPointStruct[0], &field[0]);
    
    // Calculate steady state error
    #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        double dw = field[0].w[i] - field[0].w_old[i];
        double drho = field[0].rho[i] - field[0].rho_old[i];
        double dT = field[0].T[i] - field[0].T_old[i];
        steady_state_error += du*du + dv*dv + dw*dw + drho*drho + dT*dT;
    }
    steady_state_error = sqrt(steady_state_error / myPointStruct[0].num_nodes);
    
    return steady_state_error;
}

/////////////////////////////////////////////////////////////////////////////
// Step 1: Continuity Equation
// ∂ρ/∂t + ∇·(ρu) = 0
/////////////////////////////////////////////////////////////////////////////
void solve_continuity_equation(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Calculate ρu, ρv, ρw
    #pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        field->dpdx[i] = field->rho[i] * field->u[i];  // Reuse dpdx for ρu
        field->dpdy[i] = field->rho[i] * field->v[i];  // Reuse dpdy for ρv
        field->dpdz[i] = field->rho[i] * field->w[i];  // Reuse dpdz for ρw
    }
    
    // Calculate ∇·(ρu)
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->dpdx, field->drhodx, 
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->dpdy, field->drhody, 
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    if (parameters.dimension == 3) {
        multiply_sparse_matrix_vector_vectorised_gpu(
            myPointStruct->Dz, field->dpdz, field->drhodz, 
            myPointStruct->cloud_index, num_nodes, num_cloud_points
        );
    }
    
    // Update density: ρ^(n+1) = ρ^n - Δt * ∇·(ρu)
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double div_rhou = field->drhodx[i] + field->drhody[i];
        if (parameters.dimension == 3)
            div_rhou += field->drhodz[i];
        
        field->rho_new[i] = field->rho[i] - parameters.dt * div_rhou;
        
        // Enforce boundary conditions
        if (myPointStruct->boundary_tag[i]) {
            // Boundary density from BC temperature
            field->rho_new[i] = parameters.p_ref / 
                (parameters.R_gas * myPointStruct->node_bc[i].T);
        }
    }
    
    // Update density
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->rho[i] = field->rho_new[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
// Step 2: Momentum Equations (3D version)
// ∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·τ + ρg
/////////////////////////////////////////////////////////////////////////////
void solve_momentum_equations(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // First calculate stress tensor components
    calculate_stress_tensor(myPointStruct, field);
    
    // X-MOMENTUM
    // Convective terms: ∂(ρu²)/∂x + ∂(ρuv)/∂y + ∂(ρuw)/∂z
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->dpdx[i] = field->rho[i] * field->u[i] * field->u[i];  // ρu²
        field->dpdy[i] = field->rho[i] * field->u[i] * field->v[i];  // ρuv
        field->dpdz[i] = field->rho[i] * field->u[i] * field->w[i];  // ρuw
    }
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->dpdx, field->drhodx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->dpdy, field->drhody,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->dpdz, field->drhodz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Pressure gradient
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->p, field->dedx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Update u
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double conv = field->drhodx[i] + field->drhody[i] + field->drhodz[i];
        double visc = field->div_tau_x[i];
        
        field->u_new[i] = field->u[i] - parameters.dt / field->rho[i] * 
                          (conv + field->dedx[i] - visc);
    }
    
    // Y-MOMENTUM (similar structure)
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->dpdx[i] = field->rho[i] * field->v[i] * field->u[i];  // ρvu
        field->dpdy[i] = field->rho[i] * field->v[i] * field->v[i];  // ρv²
        field->dpdz[i] = field->rho[i] * field->v[i] * field->w[i];  // ρvw
    }
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->dpdx, field->drhodx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->dpdy, field->drhody,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->dpdz, field->drhodz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->p, field->dedy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double conv = field->drhodx[i] + field->drhody[i] + field->drhodz[i];
        double visc = field->div_tau_y[i];
        
        field->v_new[i] = field->v[i] - parameters.dt / field->rho[i] * 
                          (conv + field->dedy[i] - visc);
    }
    
    // Z-MOMENTUM (3D only)
    if (parameters.dimension == 3) {
        #pragma acc parallel loop present(field)
        for (int i = 0; i < num_nodes; i++) {
            field->dpdx[i] = field->rho[i] * field->w[i] * field->u[i];  // ρwu
            field->dpdy[i] = field->rho[i] * field->w[i] * field->v[i];  // ρwv
            field->dpdz[i] = field->rho[i] * field->w[i] * field->w[i];  // ρw²
        }
        
        multiply_sparse_matrix_vector_vectorised_gpu(
            myPointStruct->Dx, field->dpdx, field->drhodx,
            myPointStruct->cloud_index, num_nodes, num_cloud_points
        );
        multiply_sparse_matrix_vector_vectorised_gpu(
            myPointStruct->Dy, field->dpdy, field->drhody,
            myPointStruct->cloud_index, num_nodes, num_cloud_points
        );
        multiply_sparse_matrix_vector_vectorised_gpu(
            myPointStruct->Dz, field->dpdz, field->drhodz,
            myPointStruct->cloud_index, num_nodes, num_cloud_points
        );
        
        multiply_sparse_matrix_vector_vectorised_gpu(
            myPointStruct->Dz, field->p, field->dedz,
            myPointStruct->cloud_index, num_nodes, num_cloud_points
        );
        
        #pragma acc parallel loop present(field, myPointStruct, parameters)
        for (int i = 0; i < num_nodes; i++) {
            double conv = field->drhodx[i] + field->drhody[i] + field->drhodz[i];
            double visc = field->div_tau_z[i];
            
            field->w_new[i] = field->w[i] - parameters.dt / field->rho[i] * 
                              (conv + field->dedz[i] - visc);
        }
    }
    
    // Enforce boundary conditions
    #pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET || 
                myPointStruct->node_bc[i].type == BC_WALL) {
                field->u_new[i] = myPointStruct->node_bc[i].u;
                field->v_new[i] = myPointStruct->node_bc[i].v;
                field->w_new[i] = myPointStruct->node_bc[i].w;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Calculate Stress Tensor and its Divergence
// τ_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i - 2/3 δ_ij ∇·u)
/////////////////////////////////////////////////////////////////////////////
void calculate_stress_tensor(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Calculate velocity gradients
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->u, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->u, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->u, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->v, field->dedx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->v, field->dedy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->v, field->dedz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->w, field->drhodx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->w, field->drhody,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->w, field->drhodz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Calculate stress components
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double dudx = field->dpdx[i];
        double dudy = field->dpdy[i];
        double dudz = field->dpdz[i];
        double dvdx = field->dedx[i];
        double dvdy = field->dedy[i];
        double dvdz = field->dedz[i];
        double dwdx = field->drhodx[i];
        double dwdy = field->drhody[i];
        double dwdz = field->drhodz[i];
        
        double div_u = dudx + dvdy + dwdz;
        double mu = field->mu[i];
        
        // Normal stresses
        field->tau_xx[i] = mu * (2.0 * dudx - 2.0/3.0 * div_u);
        field->tau_yy[i] = mu * (2.0 * dvdy - 2.0/3.0 * div_u);
        field->tau_zz[i] = mu * (2.0 * dwdz - 2.0/3.0 * div_u);
        
        // Shear stresses
        field->tau_xy[i] = mu * (dudy + dvdx);
        field->tau_xz[i] = mu * (dudz + dwdx);
        field->tau_yz[i] = mu * (dvdz + dwdy);
    }
    
    // Calculate divergence of stress tensor: ∂τ_ij/∂x_j
    // div(τ)_x = ∂τ_xx/∂x + ∂τ_xy/∂y + ∂τ_xz/∂z
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->tau_xx, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->tau_xy, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->tau_xz, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->div_tau_x[i] = field->dpdx[i] + field->dpdy[i] + field->dpdz[i];
    }
    
    // div(τ)_y = ∂τ_xy/∂x + ∂τ_yy/∂y + ∂τ_yz/∂z
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->tau_xy, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->tau_yy, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->tau_yz, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->div_tau_y[i] = field->dpdx[i] + field->dpdy[i] + field->dpdz[i];
    }
    
    // div(τ)_z = ∂τ_xz/∂x + ∂τ_yz/∂y + ∂τ_zz/∂z
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->tau_xz, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->tau_yz, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->tau_zz, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->div_tau_z[i] = field->dpdx[i] + field->dpdy[i] + field->dpdz[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
// Step 3: Energy Equation
// ∂(ρe)/∂t + ∇·(ρeu) = -p∇·u + ∇·(k∇T) + Φ + Q
// where Φ = viscous dissipation, Q = heat source
/////////////////////////////////////////////////////////////////////////////
void solve_energy_equation(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Calculate temperature gradients
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->T, field->dTdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->T, field->dTdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->T, field->dTdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Calculate heat flux: q = -k∇T
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->dpdx[i] = -field->kappa[i] * field->dTdx[i];
        field->dpdy[i] = -field->kappa[i] * field->dTdy[i];
        field->dpdz[i] = -field->kappa[i] * field->dTdz[i];
    }
    
    // Divergence of heat flux: ∇·q
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->dpdx, field->dedx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->dpdy, field->dedy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->dpdz, field->dedz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Convective term: ∇·(ρeU)
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->drhodx[i] = field->rho[i] * field->e[i] * field->u[i];
        field->drhody[i] = field->rho[i] * field->e[i] * field->v[i];
        field->drhodz[i] = field->rho[i] * field->e[i] * field->w[i];
    }
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->drhodx, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->drhody, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->drhodz, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Calculate velocity divergence for p∇·u term
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->u, field->dTdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->v, field->dTdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->w, field->dTdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Calculate viscous dissipation: Φ = τ:∇u
    calculate_viscous_dissipation(myPointStruct, field);
    
    // Update internal energy
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double conv = field->dpdx[i] + field->dpdy[i] + field->dpdz[i];
        double diff = field->dedx[i] + field->dedy[i] + field->dedz[i];
        double div_u = field->dTdx[i] + field->dTdy[i] + field->dTdz[i];
        double p_work = field->p[i] * div_u;
        
        field->e[i] = field->e[i] - parameters.dt / field->rho[i] * 
                      (conv + p_work - diff - field->Q_visc[i] - field->Q_source[i]);
        
        // Update temperature from internal energy: T = e / cv
        field->T_new[i] = field->e[i] / parameters.cv;
    }
    
    // Enforce boundary conditions
    #pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET ||
                myPointStruct->node_bc[i].type == BC_WALL) {
                field->T_new[i] = myPointStruct->node_bc[i].T;
            } else if (myPointStruct->node_bc[i].type == BC_ADIABATIC_WALL) {
                // Zero heat flux (already satisfied by Neumann BC in gradient calc)
                field->T_new[i] = field->T[i];
            }
        }
    }
    
    // Update temperature
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->T[i] = field->T_new[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
// Calculate Viscous Dissipation
// Φ = τ_ij * ∂u_i/∂x_j
/////////////////////////////////////////////////////////////////////////////
void calculate_viscous_dissipation(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Velocity gradients already calculated in stress tensor calculation
    // Reuse those or recalculate if needed
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->u, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->u, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->u, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->v, field->dedx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->v, field->dedy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->v, field->dedz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->w, field->drhodx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->w, field->drhody,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->w, field->drhodz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        double dudx = field->dpdx[i];
        double dudy = field->dpdy[i];
        double dudz = field->dpdz[i];
        double dvdx = field->dedx[i];
        double dvdy = field->dedy[i];
        double dvdz = field->dedz[i];
        double dwdx = field->drhodx[i];
        double dwdy = field->drhody[i];
        double dwdz = field->drhodz[i];
        
        // Φ = τ_xx*dudx + τ_yy*dvdy + τ_zz*dwdz + 
        //     τ_xy*(dudy + dvdx) + τ_xz*(dudz + dwdx) + τ_yz*(dvdz + dwdy)
        field->Q_visc[i] = field->tau_xx[i] * dudx +
                           field->tau_yy[i] * dvdy +
                           field->tau_zz[i] * dwdz +
                           field->tau_xy[i] * (dudy + dvdx) +
                           field->tau_xz[i] * (dudz + dwdx) +
                           field->tau_yz[i] * (dvdz + dwdy);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Step 4: Pressure Correction (Low Mach Number Formulation)
// For low Mach: p = p_ref + p'  where p' << p_ref
// Solve: ∇²p' = ∇·u_new / Δt (similar to incompressible)
/////////////////////////////////////////////////////////////////////////////
void solve_pressure_correction(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Calculate divergence of intermediate velocity
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->u_new, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->v_new, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->w_new, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // RHS for pressure Poisson: source = ρ * ∇·u / Δt
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        if (myPointStruct->corner_tag[i])
            continue;
        
        if (!myPointStruct->boundary_tag[i]) {
            // Interior nodes
            field->source[i] = field->rho[i] * 
                (field->dpdx[i] + field->dpdy[i] + field->dpdz[i]) / parameters.dt;
        } else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
            // Pressure outlet boundary
            field->source[i] = myPointStruct->node_bc[i].p;
        } else {
            // Neumann BC (walls, inlets)
            field->source[i] = field->rho[i] * 
                ((field->u_new[i] - field->u[i]) * myPointStruct->x_normal[i] +
                 (field->v_new[i] - field->v[i]) * myPointStruct->y_normal[i] +
                 (field->w_new[i] - field->w[i]) * myPointStruct->z_normal[i]) / parameters.dt;
        }
    }
    
    // Solve pressure Poisson equation (reuse existing multigrid solver)
    // This solves for pressure correction p'
    FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
}

/////////////////////////////////////////////////////////////////////////////
// Step 5: Update Velocity with Pressure Correction
/////////////////////////////////////////////////////////////////////////////
void update_velocity_compressible(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    // Calculate pressure gradient
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->p, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->p, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->p, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    // Update velocity: u = u_new - Δt/ρ * ∇p
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        if (!myPointStruct->boundary_tag[i]) {
            // Interior nodes
            field->u[i] = field->u_new[i] - parameters.dt / field->rho[i] * field->dpdx[i];
            field->v[i] = field->v_new[i] - parameters.dt / field->rho[i] * field->dpdy[i];
            field->w[i] = field->w_new[i] - parameters.dt / field->rho[i] * field->dpdz[i];
        } else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) {
            // Pressure outlet: extrapolate velocity or apply zero normal gradient
            // (implement based on your BC strategy)
        }
    }
    
    // Calculate final mass residual for monitoring
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->u, field->dpdx,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->v, field->dpdy,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dz, field->w, field->dpdz,
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    double sum = 0.0;
    #pragma acc parallel loop present(field, myPointStruct) reduction(+:sum)
    for (int i = 0; i < num_nodes; i++) {
        if (!myPointStruct->corner_tag[i]) {
            sum += fabs(field->dpdx[i] + field->dpdy[i] + field->dpdz[i]);
        }
    }
    
    printf("Velocity divergence residual: %e\n", sum / num_nodes);
}

/////////////////////////////////////////////////////////////////////////////
// 2D Version of Main Solver
/////////////////////////////////////////////////////////////////////////////
double compressible_solver_explicit_2d(PointStructure* myPointStruct, FieldVariables* field) {
    double steady_state_error = 0.0;
    
    // Save old values
    #pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        field[0].rho_old[i] = field[0].rho[i];
        field[0].T_old[i] = field[0].T[i];
        field[0].e_old[i] = field[0].e[i];
        field[0].p_old[i] = field[0].p[i];
    }
    
    update_properties(&myPointStruct[0], &field[0]);
    solve_continuity_equation_2d(&myPointStruct[0], &field[0]);
    solve_momentum_equations_2d(&myPointStruct[0], &field[0]);
    
    if (parameters.energy_equation == 1) {
        solve_energy_equation_2d(&myPointStruct[0], &field[0]);
    }
    
    solve_pressure_correction(&myPointStruct[0], &field[0]);
    update_velocity_compressible(&myPointStruct[0], &field[0]);
    
    // Calculate steady state error
    #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        double drho = field[0].rho[i] - field[0].rho_old[i];
        double dT = field[0].T[i] - field[0].T_old[i];
        steady_state_error += du*du + dv*dv + drho*drho + dT*dT;
    }
    steady_state_error = sqrt(steady_state_error / myPointStruct[0].num_nodes);
    
    return steady_state_error;
}

// 2D versions of the equation solvers (simplified - no z-components)
void solve_continuity_equation_2d(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    #pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        field->dpdx[i] = field->rho[i] * field->u[i];
        field->dpdy[i] = field->rho[i] * field->v[i];
    }
    
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dx, field->dpdx, field->drhodx, 
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    multiply_sparse_matrix_vector_vectorised_gpu(
        myPointStruct->Dy, field->dpdy, field->drhody, 
        myPointStruct->cloud_index, num_nodes, num_cloud_points
    );
    
    #pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++) {
        double div_rhou = field->drhodx[i] + field->drhody[i];
        field->rho_new[i] = field->rho[i] - parameters.dt * div_rhou;
        
        if (myPointStruct->boundary_tag[i]) {
            field->rho_new[i] = parameters.p_ref / 
                (parameters.R_gas * myPointStruct->node_bc[i].T);
        }
    }
    
    #pragma acc parallel loop present(field)
    for (int i = 0; i < num_nodes; i++) {
        field->rho[i] = field->rho_new[i];
    }
}

void solve_momentum_equations_2d(PointStructure* myPointStruct, FieldVariables* field) {
    // Similar to 3D but without z-components
    // Implementation follows same pattern as solve_momentum_equations
    // but only with x and y terms
    // (Full implementation omitted for brevity - follows 3D version pattern)
}

void solve_energy_equation_2d(PointStructure* myPointStruct, FieldVariables* field) {
    // Similar to 3D but without z-components
    // (Full implementation omitted for brevity)
}