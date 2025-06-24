# IN-WORK
1. A multigrid (multi-level) accelerated Navier Stokes solver to achieve steady state using TIMPLE a meshless framework
2. A compressible flow solver
   
# Current Status: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Heat conduction problem added with multigrid implementation and for Dirchlet boundary conditions
5. SOR solver used for Pressure Poisson and Heat conduction equations
6. Heat conduction equations solved for Neumann boundary condition
7. Navier Stokes Solver implemented with Fractional Step
8. Navier Stokes Solver implemented with Time Implicit solver
9. Multigrided Poisson solver

To run the code use a terminal and compile with gcc compiler:
command: gcc <code_file_name.c> -lm

**Make necessary changes in the grid_filenames.csv and flow_parameters.csv
grid_filenames.csv has the mesh filenames in the order from finest grid to coarse grid
flow_parameters.csv has the details like polynomial degree, phs degree, cloud_size_multiplier and others**

It looks like the `flow_parameters.csv` file has inconsistent formatting—possibly due to extra commas, missing quotes, or unstructured comments that are not in proper CSV format.

I’ll first open the file as plain text and parse it manually to extract the key-value pairs. Then I’ll read `grid_filenames.csv` similarly.


Here's what we can infer:

---

### ✅ `flow_parameters.csv`:

The file contains key simulation parameters (e.g., `Re`, `time_step`, `poly_deg`) in a simple `key,value` format. After the parameters, there are explanations written as comments using slashes (`///////////`), which are not valid CSV but give useful documentation.

---

### ✅ `grid_filenames.csv`:

The file starts with the number of mesh levels (`num_levels,1`) and then lists several Gmsh mesh files (`*.msh`), one per line, indicating supported geometries for the solver.

---

I’ll now generate a clean and informative `README.md` that includes:

* Project overview
* Folder structure
* Compilation instructions
* File descriptions for the uploaded CSVs
* Instructions on how to use your own Gmsh mesh

Let me prepare and export it now.


It seems like I can’t do more advanced data analysis right now. Please try again later.

However, I’ve already drafted the full content for the `README.md` file based on the structure of your code and the information in the uploaded `.csv` files. Here's the complete content that you can copy and save into a `README.md` file yourself:

---

### 📄 `README.md`

````markdown
# Meshless Method Navier-Stokes Solver (C Code)

This repository contains a C-based implementation of a **Meshless Method** for solving the Navier-Stokes equations. The solver utilizes a radial basis function (RBF) based meshless approach with optional GPU or multicore CPU acceleration using OpenACC.

---

## 🔧 Compilation

To compile the code for different hardware architectures:

**Single-core CPU:**
```bash
gcc mg_NS_solver.c -lm
````

**Multi-core CPU:**

```bash
nvc -acc -ta=multicore mg_NS_solver.c
```

**GPU (NVIDIA):**

```bash
nvc -acc -gpu=managed mg_NS_solver.c
```

---

## ▶️ Running the Code

After compiling:

```bash
./a.out
```

---

## 📁 Folder Structure

* `header_files/`: Contains all the custom C header files used in the solver.
* `init/`: Holds sample initialization scripts.
* `flow_parameters.csv`: Contains simulation parameters. Can be edited by the user.
* `grid_filenames.csv`: Lists available mesh files for use in the simulation.

---

## 📄 `flow_parameters.csv`

This file defines the core physical and numerical parameters for the simulation. Each line is of the format `parameter_name,value`. Here are some key parameters:

* `domain_dimensions`: Dimensionality of the domain (2 or 3).
* `poly_deg`: Degree of the appended polynomial (recommended: \[2–15]).
* `phs_deg`: Degree of the polyharmonic spline (odd integers: 3, 5, 7, 9...).
* `cloud_size_multiplier`: Controls the local stencil size; typically between 1.5 and 2.5.
* `courant_number`: CFL condition; controls the time step size.
* `steady_tolerance`: Tolerance for steady-state convergence.
* `poisson_solver_tolerance`: Tolerance for pressure Poisson solver.
* `num_vcycles`: Number of V-cycles in multigrid.
* `num_relax`: Number of relaxation steps.
* `num_time_steps`: Total simulation time steps.
* `write_interval`: How often to write output.
* `Re`: Reynolds number.
* `restart_filename`: If `restart` is 1, this file will be loaded.

You can customize these settings to suit your problem configuration.

---

## 📄 `grid_filenames.csv`

* First line: `num_levels` – Number of mesh refinement levels or hierarchies.
* Following lines: Paths to `.msh` mesh files generated using Gmsh.

### 🛠 Using Custom Mesh Files

To use your own mesh:

1. Create a mesh in Gmsh.
2. Export it in **ASCII format (version 2)**.
3. Add its path to `grid_filenames.csv`.
4. Ensure the mesh is consistent with expected domain dimensions.

---

## 🧠 Notes

* Ensure all headers in `header_files/` are correctly referenced in your code.
* Modify `flow_parameters.csv` as needed before compilation.
* For accurate results, tune parameters like `cloud_size_multiplier`, `phs_deg`, and `poisson_solver_tolerance`.

---

## 🧪 Sample Cases

Several mesh cases are included under the paths:

* `mesh/2d/TC/`
* `mesh/2d/SQ/`
* `mesh/3d/CUBE/`
* `mesh/3d/SP_in_SP/`

---

## ✍️ Author


# Meshless Method Navier-Stokes Solver (C Code)
This repository contains a C-based implementation of a **Meshless Method** for solving the Navier-Stokes equations. The solver utilizes a radial basis function (RBF) based meshless approach with optional GPU or multicore CPU acceleration using OpenACC.

## 🔧 Compilation
To compile the code for different hardware architectures:
**Single-core CPU:**
   gcc mg_NS_solver.c -lm
**Multi-core CPU:**
   nvc -acc -ta=multicore mg_NS_solver.c
**GPU (NVIDIA):**
   nvc -acc -gpu=managed mg_NS_solver.c

## ▶️ Running the Code
After compiling:
   ./a.out

## 📁 Folder Structure
* `header_files/`: Contains all the custom C header files used in the solver.
* `init/`: Holds sample initialization scripts.
* `flow_parameters.csv`: Contains simulation parameters. Can be edited by the user.
* `grid_filenames.csv`: Lists available mesh files for use in the simulation.


## 📄 `flow_parameters.csv`
This file defines the core physical and numerical parameters for the simulation. Each line is of the format `parameter_name,value`. Here are some key parameters:
* `domain_dimensions`: Dimensionality of the domain (2 or 3).
* `poly_deg`: Degree of the appended polynomial (recommended: \[2–15]).
* `phs_deg`: Degree of the polyharmonic spline (odd integers: 3, 5, 7, 9...).
* `cloud_size_multiplier`: Controls the local stencil size; typically between 1.5 and 2.5.
* `courant_number`: CFL condition; controls the time step size.
* `steady_tolerance`: Tolerance for steady-state convergence.
* `poisson_solver_tolerance`: Tolerance for pressure Poisson solver.
* `num_vcycles`: Number of V-cycles in multigrid.
* `num_relax`: Number of relaxation steps.
* `num_time_steps`: Total simulation time steps.
* `write_interval`: How often to write output.
* `Re`: Reynolds number.
* `restart_filename`: If `restart` is 1, this file will be loaded.

You can customize these settings to suit your problem configuration.


## 📄 `grid_filenames.csv`
* First line: `num_levels` – Number of mesh refinement levels or hierarchies.
* Following lines: Paths to `.msh` mesh files generated using Gmsh from refined to coarse mesh.

### 🛠 Using Custom Mesh Files
To use your own mesh:
1. Create a mesh in Gmsh.
2. Export it in **ASCII format (version 2)**.
3. Add its path to `grid_filenames.csv`.
4. Ensure the mesh is consistent with expected domain dimensions.


## 🧠 Notes
* Modify `flow_parameters.csv` as needed before compilation..


## 🧪 Sample Cases
Several mesh cases are included under the paths:
* `mesh/2d/TC/`
* `mesh/2d/SQ/`
* `mesh/3d/CUBE/`
* `mesh/3d/SP_in_SP/`

---

## ✍️ Author

