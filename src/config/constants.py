"""
Configuration constants for MeMPhyS GUI
Contains all default parameters, paths, and application settings
"""

import os
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt

# ============================================================
# Directory Configuration
# ============================================================

LOG_DIR = "./logs"
HEADER_DIR = "src/c_header_files"
OUTPUT_DIR = "./output"

# Create directories if they don't exist
for directory in [LOG_DIR, OUTPUT_DIR]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# ============================================================
# Application Options (User Configurable)
# ============================================================

DEFAULT_OPTIONS = {
    # Logging options
    "enable_gui_logging": True,
    "enable_file_logging": True,
    "enable_console_logging": True,
    
    # File management
    "auto_move_outputs": True,
    "create_dated_folders": True,
    "append_number_if_exists": True,
    
    # Monitoring
    "auto_start_convergence_monitor": True,
    "convergence_plot_auto_scale": True,
    
    # Validation
    "validate_parameters_on_run": True,
    "validate_mesh_files_on_run": True,
    
    # UI options
    "show_tooltips": True,
    "confirm_on_exit": False,
    "dark_theme": True,  # Dark theme by default
    
    # Solver options
    "auto_open_plot_on_complete": False,
    "play_sound_on_complete": False,
}

# ============================================================
# Gmsh Configuration
# ============================================================

GMSH_EXECUTABLE = "gmsh"  # Command to launch Gmsh
GMSH_FILE_EXTENSION = ".geo"
DEFAULT_GMSH_PATH = ""  # Will be set when user creates/opens a geometry

# ============================================================
# Boundary Condition Configuration
# ============================================================

BC_CSV_FILE = "bc.csv"

# Available boundary condition types
BC_TYPES = [
    "velocity_inlet",
    "pressure_outlet",
    "wall",
    "symmetry",
    "periodic",
    "outflow",
    "interior",
]

# Variables that can be set for each BC type
BC_VARIABLES = {
    "velocity_inlet": ["u", "v", "w"],
    "pressure_outlet": ["p"],
    "wall": ["u", "v", "w"],
    "symmetry": [],
    "periodic": [],
    "outflow": [],
    "interior": [],
}

# ============================================================
# Logging Configuration
# ============================================================

TODAY_STR = datetime.today().strftime("%Y-%m-%d")
LOG_FILE_PATH = os.path.join(LOG_DIR, f"log_{TODAY_STR}.txt")

# ============================================================
# Solver Parameters
# ============================================================

BASE_PARAMETERS = {
    "domain_dimensions": 2,
    "poly_deg": 3,
    "phs_deg": 3,
    "cloud_size_multiplier": 2,
    "test_derivative": 0,
    "courant_number": 0.3,
    "steady_tolerance": 1e-8,
    "poisson_solver_tolerance": 1e-8,
    "sor_parameter": 0.5,
    "time_step": 0.1,
    "num_time_steps": 100000,
    "write_interval": 50,
    "Re": 10,
}

MULTIGRID_PARAMETERS = {
    "num_vcycles": 10,
    "num_relax": 50,
}

IMPLICIT_PARAMETERS = {
    "iter_momentum": 5,
    "iter_timple": 1,
}

# Additional fixed parameters
FIXED_PARAMETERS = {
    "neumann_flag_boundary": "1",
    "facRe": "1",
    "facdt": "1",
}

RESTART_PARAMETERS = {
    "restart": "0",
    "restart_filename": "Solution.csv",
}
# ============================================================
# Solver Method Options
# ============================================================

SOLVER_METHODS = ["Fractional Step", "Time Implicit"]
Poisson_SOLVER_METHODS = ["Jacobi", "Gauss-Seidel", "BiCGStab"]
DEFAULT_SOLVER_METHOD = "Fractional Step"

# ============================================================
# File Paths and Names
# ============================================================

PARAMS_CSV = "flow_parameters.csv"
GRID_CSV = "grid_filenames.csv"
CONVERGENCE_CSV = "Convergence.csv"
SOLUTION_CSV = "Solution.csv"
DEFAULT_VTK = "Solution.vtk"

# Solver executable names
SOLVER_EXECUTABLE_WINDOWS = "solver.exe"
SOLVER_EXECUTABLE_UNIX = "solver"
SOLVER_SOURCE = "mg_NS_solver.c"

# ============================================================
# Multigrid Configuration
# ============================================================

MAX_MESH_LEVELS = 10
MIN_MESH_LEVELS = 1
DEFAULT_MESH_LEVELS = 1

# ============================================================
# GUI Configuration
# ============================================================

# Window dimensions
MAIN_WINDOW_WIDTH = 1280
MAIN_WINDOW_HEIGHT = 800
LEFT_PANEL_WIDTH = 350

# Plot dimensions
CONVERGENCE_PLOT_HEIGHT = 300

# Input widths
PARAMETER_INPUT_WIDTH = 100
MESH_PATH_INPUT_WIDTH = 200
VTK_PATH_INPUT_WIDTH = 150
COMBO_WIDTH = 150
SAVE_PATH_INPUT_WIDTH = 200

# Boundary Conditions Panel dimension
BC_WINDOW_WIDTH = 800
BC_WINDOW_HEIGHT = 400

# Button dimensions
RUN_BUTTON_WIDTH = 250
RUN_BUTTON_HEIGHT = 45
STANDARD_BUTTON_WIDTH = 80

# Log window
LOG_WINDOW_HEIGHT = 250

# Convergence plot update interval (seconds)
CONVERGENCE_UPDATE_INTERVAL = 2.0

# Auto-scroll threshold for log window
LOG_SCROLL_THRESHOLD = 5

# ============================================================
# Font Configuration
# ============================================================

FONT_SIZES = [12, 14, 16, 18, 20, 24]
DEFAULT_FONT_SIZE = 16

# Preferred system fonts by platform
FONT_PREFERENCES = {
    "Windows": ["Arial.ttf", "Calibri.ttf", "SegoeUI.ttf"],
    "Darwin": ["Helvetica.ttc", "Arial.ttf", "Menlo.ttc"],  # macOS
    "Linux": ["DejaVuSans.ttf", "LiberationSans-Regular.ttf", "DejaVuSans-Bold.ttf"]
}

# Font search paths by platform
FONT_SEARCH_PATHS = {
    "Windows": [Path(os.environ.get('WINDIR', 'C:/Windows')) / "Fonts"],
    "Darwin": [
        Path("/System/Library/Fonts"),
        Path("/Library/Fonts"),
        Path.home() / "Library/Fonts"
    ],
    "Linux": [
        Path("/usr/share/fonts/truetype"),
        Path("/usr/local/share/fonts"),
        Path.home() / ".fonts"
    ]
}

# ============================================================
# Visualization Configuration
# ============================================================

# Available colormaps from matplotlib
COLORMAPS = sorted([m for m in plt.colormaps()])
DEFAULT_COLORMAP = "viridis"

# Available variables for plotting
PLOT_VARIABLES = [
    "u",
    "v", 
    "w",
    "velocity magnitude",
    "p"
]
DEFAULT_PLOT_VARIABLE = "velocity magnitude"

# Image save formats
IMAGE_FORMATS = [".png", ".jpg", ".tif", ".pdf"]
DEFAULT_IMAGE_FORMAT = ".png"
DEFAULT_SAVE_PATH = "contour.png"

# Screenshot resolution multiplier
SCREENSHOT_SCALE = 3

# ============================================================
# File Dialog Configuration
# ============================================================

# File extensions and their colors (RGBA)
FILE_EXTENSIONS = {
    ".msh": (255, 255, 255, 255),  # Mesh files
    ".c": (255, 255, 255, 255),     # C source files
    ".vtk": (255, 255, 255, 255),   # VTK files
}

FILE_DIALOG_WIDTH = 600
FILE_DIALOG_HEIGHT = 400

# ============================================================
# Color Scheme (RGB values 0-255)
# ============================================================

COLORS_DARK = {
    "header":    (200, 220, 255),
    "subheader": (160, 185, 220),
    "success":   (145, 210, 145),
    "warning":   (255, 220, 100),
    "error":     (255, 100, 100),
    "info":      (220, 230, 255),
    "label":     (220, 230, 255),
}

COLORS_LIGHT = {
    "header":    (40,  80,  180),
    "subheader": (60,  100, 160),
    "success":   (30,  130, 30),
    "warning":   (160, 110, 0),
    "error":     (180, 30,  30),
    "info":      (50,  80,  160),
    "label":     (50,  80,  160),
}

COLORS = COLORS_DARK
# ============================================================
# Application Metadata
# ============================================================

APP_NAME = "MeMPhyS"
APP_VERSION = "2.5.1"
APP_FULL_NAME = "MeMPhyS"
APP_SUBTITLE = "Meshless Multi-Physics Solver"

# Team information
DEVELOPERS = {
    "gui": ["Akash Unnikrishnan"],
    "backend": [
        "Akash Unnikrishnan",
        "Prof. Surya Pratap Vanka",
        "Prof. Vinod Narayanan"
    ],
    "affiliation": [
        "Indian Institute of Technology Gandhinagar",
        "Faculty of Physics, University of Warsaw",
        "University of Illinois Urbana-Champaign"
    ]
}

# Help URL
HELP_URL = "https://github.com/akashu4295/Meshless_methods"

# ============================================================
# Validation Rules
# ============================================================

# Parameter constraints (name: (min, max, type))
PARAMETER_CONSTRAINTS = {
    "domain_dimensions": (2, 3, int),
    "poly_deg": (1, 10, int),
    "phs_deg": (1, 10, int),
    "cloud_size_multiplier": (1, 10, int),
    "test_derivative": (0, 3, int),
    "courant_number": (0.0, 1.0, float),
    "steady_tolerance": (1e-12, 1e-4, float),
    "poisson_solver_tolerance": (1e-12, 1e-4, float),
    "sor_parameter": (1.0, 2.0, float),
    "time_step": (1e-6, 10.0, float),
    "num_time_steps": (1, 1000000, int),
    "write_interval": (1, 10000, int),
    "Re": (0.1, 10000.0, float),
    "num_vcycles": (1, 100, int),
    "num_relax": (1, 1000, int),
    "iter_momentum": (1, 100, int),
    "iter_timple": (1, 100, int),
}

# ============================================================
# Compiler Configuration
# ============================================================

COMPILER = "gcc"
COMPILER_FLAGS = ["-lm"]

# ============================================================
# Thread Configuration
# ============================================================

SOLVER_THREAD_DAEMON = True
CONVERGENCE_THREAD_DAEMON = True

# ============================================================
# Subprocess Configuration
# ============================================================

# Buffer size for subprocess output
SUBPROCESS_BUFFER_SIZE = 1