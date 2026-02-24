"""
Configuration package for MeMPhyS GUI

This package contains all configuration constants, parameters, and theme definitions.
Import from this package to access application settings.

Example:
    from src.config.constants import BASE_PARAMETERS, LOG_DIR
    from src.config.themes import initialize_all_themes
"""

from .constants import (
    # Directories
    LOG_DIR,
    HEADER_DIR,
    OUTPUT_DIR,
    
    # Logging
    LOG_FILE_PATH,
    TODAY_STR,
    
    # Options
    DEFAULT_OPTIONS,
    
    # Parameters
    BASE_PARAMETERS,
    MULTIGRID_PARAMETERS,
    IMPLICIT_PARAMETERS,
    FIXED_PARAMETERS,
    
    # Solver
    SOLVER_METHODS,
    DEFAULT_SOLVER_METHOD,
    Poisson_SOLVER_METHODS,
    # File names
    PARAMS_CSV,
    GRID_CSV,
    CONVERGENCE_CSV,
    SOLUTION_CSV,
    DEFAULT_VTK,
    
    # GUI dimensions
    MAIN_WINDOW_WIDTH,
    MAIN_WINDOW_HEIGHT,
    LEFT_PANEL_WIDTH,
    CONVERGENCE_PLOT_HEIGHT,
    LOG_WINDOW_HEIGHT,
    
    # Boundary Condtions panel
    BC_WINDOW_HEIGHT,
    BC_WINDOW_WIDTH,

    # Input dimensions
    PARAMETER_INPUT_WIDTH,
    MESH_PATH_INPUT_WIDTH,
    VTK_PATH_INPUT_WIDTH,
    COMBO_WIDTH,
    SAVE_PATH_INPUT_WIDTH,
    
    # Button dimensions
    RUN_BUTTON_WIDTH,
    RUN_BUTTON_HEIGHT,
    STANDARD_BUTTON_WIDTH,
    
    # Font configuration
    FONT_SIZES,
    DEFAULT_FONT_SIZE,
    FONT_PREFERENCES,
    FONT_SEARCH_PATHS,
    
    # Visualization
    COLORMAPS,
    DEFAULT_COLORMAP,
    PLOT_VARIABLES,
    DEFAULT_PLOT_VARIABLE,
    DEFAULT_SAVE_PATH,
    SCREENSHOT_SCALE,
    
    # Colors
    COLORS,
    COLORS_DARK,
    COLORS_LIGHT,
    
    # Application metadata
    APP_NAME,
    APP_VERSION,
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    HELP_URL,
    
    # Multigrid
    MAX_MESH_LEVELS,
    MIN_MESH_LEVELS,
    DEFAULT_MESH_LEVELS,
    
    # Validation
    PARAMETER_CONSTRAINTS,
    
    # Gmsh and BC
    GMSH_EXECUTABLE,
    GMSH_FILE_EXTENSION,
    DEFAULT_GMSH_PATH,
    BC_CSV_FILE,
    BC_TYPES,
    BC_VARIABLES,
    
    # Compiler
    COMPILER,
    COMPILER_FLAGS,
    SOLVER_EXECUTABLE_WINDOWS,
    SOLVER_EXECUTABLE_UNIX,
    SOLVER_SOURCE,
    
    # Threading
    CONVERGENCE_UPDATE_INTERVAL,
    LOG_SCROLL_THRESHOLD,
    SUBPROCESS_BUFFER_SIZE,
    
    # File dialog
    FILE_EXTENSIONS,
    FILE_DIALOG_WIDTH,
    FILE_DIALOG_HEIGHT,
)

from .theme_registry import (
    _THEMED_TEXTS,
    themed_texts,
    update_themed_texts,
)

from .themes import (
    create_button_theme,
    create_secondary_button_theme,
    create_menu_bar_theme,
    create_input_theme,
    create_window_theme,
    create_plot_theme,
    create_combo_theme,
    create_checkbox_theme,
    create_separator_theme,
    create_modal_theme,
    create_disabled_theme,
    create_success_button_theme,
    create_error_button_theme,
    initialize_all_themes,
    apply_dark_theme,
    apply_light_theme,
    toggle_theme,
)

__all__ = [
    # Constants
    'LOG_DIR',
    'HEADER_DIR',
    'OUTPUT_DIR',
    'LOG_FILE_PATH',
    'TODAY_STR',
    'DEFAULT_OPTIONS',
    'BASE_PARAMETERS',
    'MULTIGRID_PARAMETERS',
    'IMPLICIT_PARAMETERS',
    'FIXED_PARAMETERS',
    'SOLVER_METHODS',
    'DEFAULT_SOLVER_METHOD',
    'PARAMS_CSV',
    'GRID_CSV',
    'CONVERGENCE_CSV',
    'SOLUTION_CSV',
    'DEFAULT_VTK',
    'MAIN_WINDOW_WIDTH',
    'MAIN_WINDOW_HEIGHT',
    'LEFT_PANEL_WIDTH',
    'CONVERGENCE_PLOT_HEIGHT',
    'LOG_WINDOW_HEIGHT',
    'PARAMETER_INPUT_WIDTH',
    'MESH_PATH_INPUT_WIDTH',
    'VTK_PATH_INPUT_WIDTH',
    'COMBO_WIDTH',
    'SAVE_PATH_INPUT_WIDTH',
    'RUN_BUTTON_WIDTH',
    'RUN_BUTTON_HEIGHT',
    'STANDARD_BUTTON_WIDTH',
    'FONT_SIZES',
    'DEFAULT_FONT_SIZE',
    'FONT_PREFERENCES',
    'FONT_SEARCH_PATHS',
    'COLORMAPS',
    'DEFAULT_COLORMAP',
    'PLOT_VARIABLES',
    'DEFAULT_PLOT_VARIABLE',
    'DEFAULT_SAVE_PATH',
    'SCREENSHOT_SCALE',
    'COLORS',
    'COLORS_DARK',
    'COLORS_LIGHT',
    'APP_NAME',
    'APP_VERSION',
    'APP_FULL_NAME',
    'APP_SUBTITLE',
    'DEVELOPERS',
    'HELP_URL',
    'MAX_MESH_LEVELS',
    'MIN_MESH_LEVELS',
    'DEFAULT_MESH_LEVELS',
    'PARAMETER_CONSTRAINTS',
    'COMPILER',
    'COMPILER_FLAGS',
    'SOLVER_EXECUTABLE_WINDOWS',
    'SOLVER_EXECUTABLE_UNIX',
    'SOLVER_SOURCE',
    'CONVERGENCE_UPDATE_INTERVAL',
    'LOG_SCROLL_THRESHOLD',
    'SUBPROCESS_BUFFER_SIZE',
    'FILE_EXTENSIONS',
    'FILE_DIALOG_WIDTH',
    'FILE_DIALOG_HEIGHT',
    
    # Themes
    'create_button_theme',
    'create_secondary_button_theme',
    'create_menu_bar_theme',
    'create_input_theme',
    'create_window_theme',
    'create_plot_theme',
    'create_combo_theme',
    'create_checkbox_theme',
    'create_separator_theme',
    'create_modal_theme',
    'create_disabled_theme',
    'create_success_button_theme',
    'create_error_button_theme',
    'initialize_all_themes',
    'apply_dark_theme',
    'apply_light_theme',
    'toggle_theme',

    # REGISTRY
    '_THEMED_TEXTS',
    'themed_texts',
    'update_themed_texts',

]