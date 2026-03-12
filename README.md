# Meshless Method Navier-Stokes Solver (C Code)

This repository contains a C-based implementation of a **Meshless Method** for solving the Navier-Stokes equations. The solver utilizes a radial basis function (RBF) based meshless approach with optional multicore CPU acceleration using OpenACC.

---

### Quick Notes on status of the code

### IN-WORK/ TODO
1. A multigrid (multi-level) accelerated Navier Stokes solver to achieve steady state using TIMPLE in meshless framework
2. A compressible flow solver
   
### Current Status: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Heat conduction problem added with multigrid implementation and for Dirchlet boundary conditions
5. SOR solver used for Pressure Poisson and Heat conduction equations
6. Heat conduction equations solved for Neumann boundary condition
7. Navier Stokes Solver implemented with Fractional Step
8. Navier Stokes Solver implemented with Time Implicit solver
9. Multigrided Poisson solver
10. Various boundary conditions implemented.
11. GPU parallelised with OpenACC

**Make necessary changes in the grid_filenames.csv and flow_parameters.csv
grid_filenames.csv has the mesh filenames in the order from finest grid to coarse grid
flow_parameters.csv has the details like polynomial degree, phs degree, cloud_size_multiplier and others**

---

### Compilation

To compile the code for different hardware architectures:

**Single-core CPU:**
```bash
gcc @sources.txt -lm
```

**Multi-core CPU:**

```bash
nvc -acc -ta=multicore @sources.txt
```

**GPU (NVIDIA):**

```bash
nvc -acc -gpu=managed @sources.txt
```

### Running the Code

After compiling:

```bash
./a.out
```

## Installation of GUI wrapper with python
Since we are in the developing stage, we are using a conda environment. If you don't have conda installed in the system, please do the following:

### Linux
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

rm ./Miniconda3-latest-Linux-x86_64.sh
```

### Windows (Powershell commands)
```bash
curl -o Miniconda3-latest-Windows-x86_64.exe `
     https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

Start-Process -Wait .\Miniconda3-latest-Windows-x86_64.exe `
  -ArgumentList "/InstallationType=JustMe",
                "/AddToPath=0",
                "/RegisterPython=0",
                "/S",
                "/D=$env:USERPROFILE\Miniconda3"

$env:USERPROFILE\Miniconda3\Scripts\conda.exe init powershell

Remove-Item Miniconda3-latest-Windows-x86_64.exe 

conda install -c conda-forge m2w64-gcc
```

After installing conda environment run this (for both Linux and Windows system):
```bash
conda env create -f environment.yml
conda activate memphys_gui
```
To start the GUI:
```bash
python gui.py
```
---

### Folder Structure

* `header_files/`: Contains all the custom C header files used in the solver.
* `init/`: Holds sample initialization scripts.
* `flow_parameters.csv`: Contains simulation parameters. Can be edited by the user.
* `grid_filenames.csv`: Lists available mesh files for use in the simulation.

---

### `flow_parameters.csv`

This file defines the core physical and numerical parameters for the simulation. Each line is of the format `parameter_name,value`. Here are some key parameters:

* `domain_dimensions`: Dimensionality of the domain (2 or 3).
* `poly_deg`: Degree of the appended polynomial (recommended: \[2–15]).
* `phs_deg`: Degree of the polyharmonic spline (odd integers: 3, 5, 7, 9...).
* `cloud_size_multiplier`: Controls the local stencil size; typically between 1.5 and 2.5 for 2d cases and between 3 to 4 for 3d cases.
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

### `grid_filenames.csv`

* First line: `num_levels` – Number of mesh refinement levels or hierarchies.
* Following lines: Paths to `.msh` mesh files generated using Gmsh from fine grid to coarse grid.

#### Using Custom Mesh Files

To use your own mesh:

1. Create a mesh in Gmsh.
2. Export it in **ASCII format (version 2)**.
3. Add its path to `grid_filenames.csv`.
4. Ensure the mesh is consistent with expected domain dimensions.

---

### Notes

* Modify `flow_parameters.csv` as needed before compilation.

---

### Sample Cases

Several 2d and 3d mesh cases are included under the paths:

* `mesh/2d/TC/`
* `mesh/2d/SQ/`
* `mesh/3d/CUBE/`
* `mesh/3d/SP_in_SP/`


# MeMPhyS GUI Guide

## [TODO]MeMPhyS GUI Restructuring 

### High Priority
- [o] Test full solver workflow (compile → run → plot)
- [o] Verify mesh file validation
- [ ] Test all callbacks thoroughly
- [ ] Remove debug print statements from main.py

### Medium Priority
- [ ] Enable custom fonts (set `ENABLE_CUSTOM_FONTS = True`)
- [ ] Implement config save/load functionality
- [ ] Add unit tests for modules
- [ ] Add more validation rules

### Low Priority
- [ ] Advanced plot settings dialog
- [ ] Keyboard shortcuts
- [ ] Drag-and-drop mesh files
- [ ] Recent files menu

## Known Limitations

1. **Font Changes**: Require application restart (DearPyGUI limitation)
2. **Screenshot Feature**: External plotter window only
3. **Config Save/Load**: Not yet implemented (placeholders exist)

## New Directory Structure

```
MeMPhyS/
├── main.py                      # New entry point (replaces old gui.py)
├── src/
│   ├── config/                  # Configuration and constants
│   │   ├── __init__.py
│   │   ├── constants.py         # All constants and parameters
│   │   └── themes.py            # DearPyGUI themes
│   │
│   ├── core/                    # Core functionality
│   │   ├── __init__.py
│   │   ├── state.py             # Application state management
│   │   └── logger.py            # Logging system
│   │
│   ├── utils/                   # Utility functions
│   │   ├── __init__.py
│   │   ├── fonts.py             # Font management
│   │   ├── file_io.py           # File operations
│   │   └── platform_utils.py   # OS-specific utilities
│   │
│   ├── solver/                  # Solver management
│   │   ├── __init__.py
│   │   ├── runner.py            # Compilation and execution
│   │   └── monitoring.py        # Convergence monitoring
│   │
│   ├── callbacks/               # Callback functions
│   │   ├── __init__.py
│   │   ├── solver_callbacks.py  # Solver operations
│   │   ├── mesh_callbacks.py    # Mesh management
│   │   ├── plot_callbacks.py    # Visualization
│   │   └── menu_callbacks.py    # Menu actions
│   │
│   ├── ui/                      # User interface
│   │   ├── __init__.py
│   │   ├── main_window.py       # Main window assembly
│   │   ├── menu_bar.py          # Menu bar
│   │   ├── parameters_panel.py  # Left panel
│   │   ├── visualization_panel.py # Right panel
│   │   └── dialogs.py           # Modal dialogs
│   │
│   └── plotting/
│       └── plotter.py           # External plotting (existing)
│
├── header_files/                # C header files (unchanged)
├── logs/                        # Log files
└── requirements.txt             # Python dependencies
```

## Quick Start

### Running the Application

```bash
python gui.py
```

## Module Descriptions

### 1. **src/config/** - Configuration

#### constants.py
Contains all application constants:
- Solver parameters (BASE_PARAMETERS, IMPLICIT_PARAMETERS, MULTIGRID_PARAMETERS)
- GUI dimensions and settings
- File paths and names
- Colors and styling
- Validation rules (PARAMETER_CONSTRAINTS)

#### themes.py
All DearPyGUI themes:
- Button themes (primary, secondary, success, error)
- Input themes
- Plot themes
- Modal themes
- `initialize_all_themes()` - Create all themes
- `apply_global_theme()` - Apply app-wide styling

### 2. **src/core/** - Core Functionality

#### state.py
Application state management via `AppState` class:
- Process management (plotter, solver)
- File handle management
- Font registry
- UI state (multigrid enabled, mesh levels, etc.)
- `app_state` - Global instance

**Key Methods:**
```python
app_state.solver_running = True
app_state.set_mesh_file(1, "mesh.msh")
app_state.cleanup()  # Call on exit
```

#### logger.py
Structured logging system via `Logger` class:
- Multiple outputs (file, GUI, console)
- Log levels (INFO, WARNING, ERROR, SUCCESS, DEBUG)
- Colored console output
- Auto-scrolling GUI log
- `logger` - Global instance

**Key Methods:**
```python
logger.info("Information message")
logger.error("Error message")
logger.success("Success message")
logger.log_exception(e, "Context")
logger.separator()
```

### 3. **src/utils/** - Utilities

#### fonts.py
Font management:
- `find_system_font()` - Cross-platform font discovery
- `initialize_fonts()` - Load all fonts at startup
- `change_font()` - Change current font
- Platform-specific font paths

#### file_io.py
File operations:
- `write_parameters_csv()` - Write parameters to CSV
- `write_grid_csv()` - Write grid configuration
- `validate_mesh_file()` - Validate mesh files
- `read_csv_file()` / `write_csv_file()` - CSV utilities

#### platform_utils.py
OS-specific utilities:
- `open_folder()` - Open in file explorer
- `open_url()` - Open in browser
- `get_platform()` - Get OS name
- `build_compile_command()` - Build compilation command

### 4. **src/solver/** - Solver Management

#### runner.py
Solver compilation and execution via `SolverRunner` class:
- `compile_and_run()` - Main entry point
- `stop_solver()` - Terminate running solver
- `check_dependencies()` - Verify compiler, files exist
- Asynchronous execution with live output
- `solver_runner` - Global instance

**Usage:**
```python
from src.solver import solver_runner

success = solver_runner.compile_and_run(
    init_file="init_cavity.c",
    button_tag="run_button",
    on_complete=my_callback
)
```

#### monitoring.py
Convergence monitoring via `ConvergenceMonitor` class:
- Real-time CSV file monitoring
- Automatic plot updates
- Smart axis scaling
- Convergence status tracking
- `convergence_monitor` - Global instance

**Usage:**
```python
from src.solver import convergence_monitor

convergence_monitor.start()
status = convergence_monitor.get_convergence_status()
convergence_monitor.stop()
```

### 5. **src/callbacks/** - Callbacks

All callback functions organized by function:

#### solver_callbacks.py
- `run_solver_callback()` - Run solver
- `validate_numeric_input()` - Parameter validation
- `show_implicit_callback()` - Show/hide implicit params
- `stop_solver_callback()` - Stop solver

#### mesh_callbacks.py
- `show_multigrid_callback()` - Enable/disable multigrid
- `update_mesh_inputs_callback()` - Update mesh input visibility
- `select_mesh_file_callback()` - File selection
- `validate_all_mesh_files_callback()` - Validate all meshes
- `auto_fill_mesh_levels_callback()` - Auto-generate mesh names

#### plot_callbacks.py
- `update_plot_callback()` - Launch plotter
- `reset_convergence_plot_callback()` - Reset plot
- `export_convergence_data_callback()` - Export data to CSV

#### menu_callbacks.py
- `open_logs_callback()` - Open logs folder
- `show_about_callback()` - Show About dialog
- `show_preferences_callback()` - Show Preferences
- `exit_application_callback()` - Clean exit

### 6. **src/ui/** - User Interface

#### main_window.py
- `create_main_window()` - Assemble main window

#### parameters_panel.py
- `create_parameters_panel()` - Left panel with parameters

#### visualization_panel.py
- `create_visualization_panel()` - Right panel with plots

#### dialogs.py
- `create_about_dialog()` - About dialog
- `create_preferences_dialog()` - Preferences dialog
- `show_error_dialog()` - Error message
- `show_confirmation_dialog()` - Yes/No confirmation


## Configuration

### Enable Custom Fonts
In `main.py`, line ~47:
```python
ENABLE_CUSTOM_FONTS = True  # Change to True
```

### Adjust Convergence Update Interval
In `src/config/constants.py`:
```python
CONVERGENCE_UPDATE_INTERVAL = 2.0  # seconds
```

### Change Log Level
In your code:
```python
logger.set_enable_console(True)   # Console output
logger.set_enable_file(True)      # File output
logger.set_enable_gui(True)       # GUI output
```

## Code Examples

### Using the Logger
```python
from src.core import logger

logger.info("Information message")
logger.success("Success message")
logger.error("Error message")
logger.warning("Warning message")
logger.separator()
```

### Using App State
```python
from src.core import app_state

app_state.solver_running = True
app_state.set_mesh_file(1, "mesh.msh")
app_state.cleanup()  # On exit
```

### Using Solver
```python
from src.solver import solver_runner

solver_runner.compile_and_run(
    init_file="init_cavity.c",
    button_tag="run_button"
)
```

### Some default parameters (non-dimensionalisation)
Density, rho = 1
Viscosity, nu = 1/Re

## Troubleshooting

### Convergence plot not updating
- Check if `Convergence.csv` exists
- Verify convergence monitor started (check logs)
- Run solver to generate data

### Fonts not loading
- Check available fonts for your OS
- See `src/utils/fonts.py` for font paths
- Try different font from Preferences

### Solver won't compile
- Verify gcc is installed: `gcc --version`
- Check init file path is correct
- Check mesh files exist
- Review logs for compilation errors

## Getting Help

1. Check logs: `logs/log_YYYY-MM-DD.txt`
2. Check GitHub issues
3. Enable debug logging


### Author

Dr. Akash Unnikrishnan developed this code as part of his PhD work with the guidance of Prof. Surya Pratap Vanka from University of Illinois at Urbana Champaign and Prof. Vinod Narayanan at Indian Institute of Technology Gandhinagar. Special mention has to go for Dr. Shantanu Shahane, who developed the first version of memphys which is available here: github.com/shahaneshantanu/memphys


---

**Version**: 2.2 Restructured  
**Date**: December 2025  
**Maintainer**: Akash Unnikrishnan

---
