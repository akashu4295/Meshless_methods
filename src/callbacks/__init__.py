"""
Callbacks package for MeMPhyS GUI

This package contains organized callback functions for different
parts of the application:
- solver_callbacks: Solver execution and parameter handling
- mesh_callbacks: Mesh file management and multigrid
- plot_callbacks: Visualization and plotting
- menu_callbacks: Menu bar actions

Example:
    from src.callbacks import (
        run_solver_callback,
        show_multigrid_callback,
        update_plot_callback
    )
"""

# Solver callbacks
from .solver_callbacks import (
    run_solver_callback,
    on_solver_complete,
    validate_numeric_input,
    show_implicit_callback,
    write_files_callback,
    stop_solver_callback,
    check_solver_status_callback,
    validate_all_parameters_callback,
)

# Mesh callbacks
from .mesh_callbacks import (
    show_multigrid_callback,
    update_mesh_inputs_callback,
    open_file_dialog_callback,
    select_mesh_file_callback,
    validate_selected_mesh_file,
    validate_all_mesh_files_callback,
    clear_mesh_files_callback,
    browse_init_file_callback,
    select_init_file_callback,
    copy_mesh_path_callback,
    auto_fill_mesh_levels_callback,
)

# Plot callbacks
from .plot_callbacks import (
    update_plot_callback,
    save_plot_image_callback,
    reset_convergence_plot_callback,
    force_convergence_update_callback,
    toggle_convergence_monitor_callback,
    show_convergence_status_callback,
    export_convergence_data_callback,
    browse_vtk_file_callback,
    select_vtk_file_callback,
    change_colormap_callback,
    change_variable_callback,
    open_plot_settings_callback,
    close_all_plot_windows_callback,
    browse_output_folder_callback,
    set_vtk_from_latest_output_callback,
    open_in_paraview_callback,
    paraview_selected_callback,
)

# Menu callbacks
from .menu_callbacks import (
    open_logs_callback,
    open_help_callback,
    show_about_callback,
    show_preferences_callback,
    show_options_callback,
    apply_preferences_callback,
    exit_application_callback,
    open_config_callback,
    save_config_callback,
    load_config_file_callback,
    save_config_file_callback,
    show_system_info_callback,
    show_application_state_callback,
    check_for_updates_callback,
    report_bug_callback,
    toggle_theme_callback,
)

# Geometry callbacks
from .geometry_callbacks import (
    launch_gmsh_callback,
    browse_geometry_file_callback,
    select_geometry_file_callback,
    set_mesh_from_geometry_callback,
)

# Boundary Conditions callbacks
from .bc_window_callbacks import (
    show_bc_window_callback,
)

__all__ = [
    # Solver callbacks
    'run_solver_callback',
    'on_solver_complete',
    'validate_numeric_input',
    'show_implicit_callback',
    'write_files_callback',
    'stop_solver_callback',
    'check_solver_status_callback',
    'validate_all_parameters_callback',
    
    # Mesh callbacks
    'show_multigrid_callback',
    'update_mesh_inputs_callback',
    'open_file_dialog_callback',
    'select_mesh_file_callback',
    'validate_selected_mesh_file',
    'validate_all_mesh_files_callback',
    'clear_mesh_files_callback',
    'browse_init_file_callback',
    'select_init_file_callback',
    'copy_mesh_path_callback',
    'auto_fill_mesh_levels_callback',
    
    # Plot callbacks
    'update_plot_callback',
    'save_plot_image_callback',
    'reset_convergence_plot_callback',
    'force_convergence_update_callback',
    'toggle_convergence_monitor_callback',
    'show_convergence_status_callback',
    'export_convergence_data_callback',
    'browse_vtk_file_callback',
    'select_vtk_file_callback',
    'change_colormap_callback',
    'change_variable_callback',
    'open_plot_settings_callback',
    'close_all_plot_windows_callback',
    'browse_output_folder_callback',
    'set_vtk_from_latest_output_callback',
    'paraview_selected_callback',
    'open_in_paraview_callback',
    
    # Menu callbacks
    'open_logs_callback',
    'open_help_callback',
    'show_about_callback',
    'show_preferences_callback',
    'show_options_callback',
    'apply_preferences_callback',
    'exit_application_callback',
    'open_config_callback',
    'save_config_callback',
    'load_config_file_callback',
    'save_config_file_callback',
    'show_system_info_callback',
    'show_application_state_callback',
    'check_for_updates_callback',
    'report_bug_callback',
    'toggle_theme_callback',

    # boundary conditions callbacks
    'show_bc_window_callback',
]