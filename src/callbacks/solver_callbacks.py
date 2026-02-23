"""
Solver-related callbacks for MeMPhyS GUI

Handles callbacks for:
- Running the solver
- Parameter validation
- File writing
- Solver method changes
"""

import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils import write_parameters_csv, write_grid_csv
from src.solver import solver_runner, convergence_monitor
from src.config import BASE_PARAMETERS, IMPLICIT_PARAMETERS, MULTIGRID_PARAMETERS


def run_solver_callback(sender, app_data, user_data):
    """
    Main callback to compile and run the solver
    
    Args:
        sender: Button tag that triggered the callback
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.header("Starting New Solver Run")
    init_on = dpg.get_value("init_toggle")
    gpu_on = dpg.get_value("gpu_toggle")
    # Check if already running
    if solver_runner.is_solver_running():
        logger.warning("Solver is already running. Please wait for completion.")
        return
    
    # Get initialization file path
    if not dpg.does_item_exist("init_path") and init_on:
        logger.error("Initialization file path widget not found")
        return
    
    init_file = dpg.get_value("init_path")
    
    if not init_file and init_on:
        logger.error("No initialization file specified")
        logger.info("Please select an initialization file using the Browse button")
        return
    
    # Write configuration files
    logger.info("Writing configuration files...")
    
    grid_success = write_grid_csv()
    if not grid_success:
        logger.error("Failed to write grid configuration. Cannot continue.")
        return
    
    params_success = write_parameters_csv()
    if not params_success:
        logger.error("Failed to write parameters. Cannot continue.")
        return
    
    logger.success("Configuration files written successfully")
    
    # Check dependencies
    deps_ok, missing = solver_runner.check_dependencies()
    if not deps_ok:
        logger.error("Missing required dependencies:")
        for item in missing:
            logger.error(f"  - {item}")
        return
    
    # Start compilation and execution
    if gpu_on:
        compiler = "nvc"
    else:
        compiler = "gcc"
    success = solver_runner.compile_and_run(
        init_file=init_file,
        compiler=compiler,
        button_tag=sender,
        on_complete=on_solver_complete
    )
    
    if not success:
        logger.error("Failed to start solver")

def on_solver_complete(returncode: int):
    """
    Callback when solver completes
    
    Args:
        returncode: Process return code
    """
    if returncode == 0:
        logger.success("Solver run completed successfully")
        
        # Organize output files
        from src.utils.output_manager import organize_solver_outputs
        organize_solver_outputs()
    else:
        logger.error(f"Solver run failed with return code {returncode}")


def validate_numeric_input(sender, app_data, user_data):
    """
    Validate that a parameter input contains a valid number
    
    Args:
        sender: Input widget tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    try:
        value = dpg.get_value(sender)
        
        # Try to convert to float
        float(value)
        
    except ValueError:
        # Invalid input - revert to default
        param_name = sender.replace("param_", "")
        
        # Find the parameter in one of the dictionaries
        default_value = None
        
        if param_name in BASE_PARAMETERS:
            default_value = BASE_PARAMETERS[param_name]
        elif param_name in IMPLICIT_PARAMETERS:
            default_value = IMPLICIT_PARAMETERS[param_name]
        elif param_name in MULTIGRID_PARAMETERS:
            default_value = MULTIGRID_PARAMETERS[param_name]
        
        if default_value is not None:
            dpg.set_value(sender, str(default_value))
            logger.warning(f"Invalid input for {param_name}. Reverted to default: {default_value}")
        else:
            logger.error(f"Invalid input for {param_name} and no default found")


def show_implicit_callback(sender, app_data, user_data):
    """
    Show/hide implicit solver parameters based on solver method selection
    
    Args:
        sender: Combo box tag
        app_data: Selected value
        user_data: User data (unused)
    """
    solver_method = dpg.get_value("solver_method")
    show_implicit = (solver_method == "Time Implicit")
    
    # Show/hide implicit parameters
    for pname in IMPLICIT_PARAMETERS.keys():
        tag = f"param_{pname}"
        if dpg.does_item_exist(tag):
            dpg.configure_item(tag, show=show_implicit)
    
    # Update app state
    app_state.current_solver_method = solver_method
    
    logger.debug(f"Solver method changed to: {solver_method}")


def write_files_callback(sender, app_data, user_data):
    """
    Write parameter and grid files without running solver
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Writing configuration files...")
    
    params_success = write_parameters_csv()
    grid_success = write_grid_csv()
    
    if params_success and grid_success:
        logger.success("All files written successfully")
    elif params_success:
        logger.warning("Parameters written, but grid file failed")
    elif grid_success:
        logger.warning("Grid file written, but parameters failed")
    else:
        logger.error("Failed to write configuration files")


def stop_solver_callback(sender, app_data, user_data):
    """
    Stop the currently running solver
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if solver_runner.is_solver_running():
        logger.warning("Stopping solver...")
        success = solver_runner.stop_solver()
        
        if success:
            logger.info("Solver stopped successfully")
        else:
            logger.error("Failed to stop solver")
    else:
        logger.info("No solver is currently running")


def check_solver_status_callback(sender, app_data, user_data):
    """
    Check and display solver status
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if solver_runner.is_solver_running():
        logger.info("Solver Status: RUNNING")
        
        # Get convergence status
        status = convergence_monitor.get_convergence_status()
        
        if status["has_data"]:
            logger.info(f"Iterations: {status['iterations']}")
            if status['current_error']:
                logger.info(f"Current Error: {status['current_error']:.6e}")
            if status['min_error']:
                logger.info(f"Minimum Error: {status['min_error']:.6e}")
            logger.info(f"Converging: {'Yes' if status['is_converging'] else 'No'}")
    else:
        logger.info("Solver Status: NOT RUNNING")


def validate_all_parameters_callback(sender, app_data, user_data):
    """
    Validate all parameter inputs
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Validating all parameters...")
    
    from src.config import PARAMETER_CONSTRAINTS
    
    invalid_params = []
    
    # Check all parameters with constraints
    for param_name, (min_val, max_val, param_type) in PARAMETER_CONSTRAINTS.items():
        tag = f"param_{param_name}"
        
        if not dpg.does_item_exist(tag):
            continue
        
        try:
            value_str = dpg.get_value(tag)
            value = param_type(value_str)
            
            if not (min_val <= value <= max_val):
                invalid_params.append(
                    f"{param_name}: {value} not in range [{min_val}, {max_val}]"
                )
                logger.warning(f"Parameter {param_name} out of range: {value}")
        
        except (ValueError, TypeError):
            invalid_params.append(f"{param_name}: invalid value '{value_str}'")
            logger.error(f"Invalid value for {param_name}: {value_str}")
    
    if invalid_params:
        logger.warning(f"Found {len(invalid_params)} invalid parameters")
        for msg in invalid_params:
            logger.warning(f"  - {msg}")
    else:
        logger.success("All parameters are valid")

def show_restart_callback(sender, app_data, user_data):
    dpg.configure_item("restart_group", show=app_data)

def show_init_callback(sender, app_data, user_data):
    dpg.configure_item("init_group", show=app_data)
