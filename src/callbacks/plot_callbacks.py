"""
Visualization and plotting callbacks for MeMPhyS GUI

Handles callbacks for:
- Contour plotting
- Plot updates
- Image saving
- Convergence plot control
"""

import os
import sys
import subprocess
import dearpygui.dearpygui as dpg
import sys
import numpy as np

from src.core import logger, app_state
from src.solver import convergence_monitor
from src.config import DEFAULT_VTK, DEFAULT_SAVE_PATH
from src.utils import validate_file_path

def update_plot_callback(sender, app_data, user_data):
    """
    Update contour plot using a subprocess (no external plotter.py file)
    """
    import os
    import sys
    import subprocess

    vtk_path = dpg.get_value("contour_vtk_path") or DEFAULT_VTK
    var_choice = dpg.get_value("contour_var") or "velocity magnitude"
    cmap = dpg.get_value("contour_cmap") or "viridis"

    dimension = (
        dpg.get_value("param_domain_dimensions")
        if dpg.does_item_exist("param_domain_dimensions")
        else "2"
    )

    if not validate_file_path(vtk_path, must_exist=True):
        from src.utils.output_manager import get_latest_output_folder

        latest_folder = get_latest_output_folder()
        if latest_folder:
            alt_path = os.path.join(latest_folder, os.path.basename(vtk_path))
            if os.path.exists(alt_path):
                vtk_path = alt_path
                logger.info(f"Using VTK file from output folder: {vtk_path}")
            else:
                logger.error(f"VTK file not found: {os.path.basename(vtk_path)}")
                return
        else:
            logger.error(f"VTK file not found: {vtk_path}")
            return

    logger.info(f"Plotting {var_choice} from {vtk_path}")
    
    # Kill previous process
    if app_state.plotter_process is not None:
        try:
            if app_state.plotter_process.poll() is None:
                app_state.plotter_process.terminate()
                app_state.plotter_process.wait(timeout=2)
                logger.info("Previous plotter closed")
        except Exception as e:
            logger.warning(f"Error closing previous plotter: {e}")

    plot_code = f"""
import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter
import sys

vtk_path = r\"{vtk_path}\"
var_choice = r\"{var_choice}\"
cmap = r\"{cmap}\"
dimension = r\"{dimension}\"

mesh = pv.read(vtk_path)

if var_choice in ("u", "v", "w"):
    idx = {{"u": 0, "v": 1, "w": 2}}[var_choice]
    plot_data = mesh["velocity"][:, idx]
elif var_choice == "velocity magnitude":
    vectors = mesh["velocity"]
    plot_data = np.linalg.norm(vectors, axis=1)
elif var_choice == "p":
    plot_data = mesh["pressure"]
else:
    raise RuntimeError(f"Variable '{{var_choice}}' not found")

plotter = BackgroundPlotter()
plotter.add_mesh(
    mesh,
    scalars=plot_data,
    cmap=cmap,
    scalar_bar_args={{"title": var_choice}},
)

if dimension == "2":
    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.reset_camera()

plotter.show()
plotter.app.exec_()
"""

    try:
        proc = subprocess.Popen(
            [sys.executable, "-c", plot_code],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        app_state.plotter_process = proc
        logger.success("Plotter window opened")

    except Exception as e:
        logger.log_exception(e, "Error launching plotter")

def open_in_paraview_callback(sender, app_data, user_data):
    """
    Open current VTK file in ParaView.
    If ParaView is not found, show browse popup.
    """
    import os
    import subprocess
    import shutil

    vtk_path = dpg.get_value("contour_vtk_path") or DEFAULT_VTK

    # Validate file
    if not validate_file_path(vtk_path, must_exist=True):
        logger.error(f"VTK file not found: {vtk_path}")
        return

    # If user previously selected ParaView
    if app_state.paraview_path and os.path.exists(app_state.paraview_path):
        paraview_exe = app_state.paraview_path
    else:
        # Try to find in PATH
        paraview_exe = shutil.which("paraview")

    if paraview_exe is None:
        dpg.configure_item("paraview_popup", show=True)
        return

    try:
        subprocess.Popen(
            [paraview_exe, vtk_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        logger.success("Opened file in ParaView")

    except Exception as e:
        logger.log_exception(e, "Failed to launch ParaView")

def paraview_selected_callback(app_data):
    """
    Store selected ParaView executable and close popup.
    """
    selected_path = app_data["file_path_name"]

    if not selected_path:
        return

    app_state.paraview_path = selected_path
    logger.success(f"ParaView set to: {selected_path}")

    dpg.configure_item("paraview_popup", show=False)

def reset_convergence_plot_callback(sender, app_data, user_data):
    """
    Reset the convergence plot
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    convergence_monitor.reset()
    logger.info("Convergence plot reset")


def force_convergence_update_callback(sender, app_data, user_data):
    """
    Force an immediate update of the convergence plot
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    convergence_monitor.force_update()
    logger.info("Convergence plot updated")


def toggle_convergence_monitor_callback(sender, app_data, user_data):
    """
    Toggle convergence monitoring on/off
    
    Args:
        sender: Checkbox tag
        app_data: Checkbox state
        user_data: User data (unused)
    """
    enabled = app_data
    
    if enabled:
        if not convergence_monitor.is_running():
            convergence_monitor.start()
            logger.info("Convergence monitoring enabled")
    else:
        if convergence_monitor.is_running():
            convergence_monitor.stop()
            logger.info("Convergence monitoring disabled")


def show_convergence_status_callback(sender, app_data, user_data):
    """
    Display detailed convergence status
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    status = convergence_monitor.get_convergence_status()
    
    logger.separator()
    logger.info("Convergence Status:")
    
    if status["has_data"]:
        logger.info(f"  Iterations: {status['iterations']}")
        
        if status['current_error'] is not None:
            logger.info(f"  Current Error: {status['current_error']:.6e}")
        
        if status['min_error'] is not None:
            logger.info(f"  Minimum Error: {status['min_error']:.6e}")
        
        converging_str = "Yes" if status.get('is_converging', False) else "No"
        logger.info(f"  Converging: {converging_str}")
    else:
        logger.info("  No convergence data available")
    
    logger.separator()


def export_convergence_data_callback(sender, app_data, user_data):
    """
    Export convergence data to a file
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    data = convergence_monitor.get_latest_data()
    
    if data is None:
        logger.warning("No convergence data available to export")
        return
    
    x, y = data
    
    # Create export filename with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"convergence_export_{timestamp}.csv"
    
    try:
        import pandas as pd
        df = pd.DataFrame({"Iteration": x, "Error": y})
        df.to_csv(filename, index=False)
        
        logger.success(f"Convergence data exported to {filename}")
        logger.info(f"Exported {len(x)} data points")
    
    except Exception as e:
        logger.log_exception(e, "Error exporting convergence data")


def browse_output_folder_callback(sender, app_data, user_data):
    """
    Open the output folder in file explorer
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils import open_folder
    from src.utils.output_manager import get_latest_output_folder
    from src.config import OUTPUT_DIR
    
    # Try to open latest output folder
    latest_folder = get_latest_output_folder()
    
    if latest_folder and os.path.exists(latest_folder):
        success = open_folder(latest_folder)
        if success:
            logger.info(f"Opened latest output folder: {latest_folder}")
    else:
        # Fall back to main output directory
        if os.path.exists(OUTPUT_DIR):
            success = open_folder(OUTPUT_DIR)
            if success:
                logger.info(f"Opened output directory: {OUTPUT_DIR}")
        else:
            logger.warning("No output folder found yet. Run the solver first.")


def set_vtk_from_latest_output_callback(sender, app_data, user_data):
    """
    Automatically set VTK path from latest output folder
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils.output_manager import get_latest_output_folder
    from src.config import DEFAULT_VTK
    
    latest_folder = get_latest_output_folder()
    
    if not latest_folder:
        logger.warning("No output folder found. Run the solver first.")
        return
    
    # Look for VTK files in latest output
    vtk_files = []
    for file in os.listdir(latest_folder):
        if file.endswith('.vtk'):
            vtk_files.append(file)
    
    if not vtk_files:
        logger.warning(f"No VTK files found in {latest_folder}")
        return
    
    # Use the first VTK file found (or Solution.vtk if it exists)
    if DEFAULT_VTK in vtk_files:
        vtk_file = DEFAULT_VTK
    else:
        vtk_file = vtk_files[0]
    
    vtk_path = os.path.join(latest_folder, vtk_file)
    
    if dpg.does_item_exist("contour_vtk_path"):
        dpg.set_value("contour_vtk_path", vtk_path)
        logger.success(f"VTK path set to: {vtk_path}")
    else:
        logger.error("Could not set VTK path widget")


def select_vtk_file_callback(sender, app_data, user_data):
    """
    Callback when VTK file is selected
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: Target input field tag
    """
    if 'file_path_name' not in app_data:
        return
    
    file_path = app_data['file_path_name']
    
    if dpg.does_item_exist("contour_vtk_path"):
        dpg.set_value("contour_vtk_path", file_path)
        app_state.current_vtk_file = file_path
        logger.info(f"VTK file selected: {file_path}")


def change_colormap_callback(sender, app_data, user_data):
    """
    Callback when colormap is changed
    
    Args:
        sender: Combo box tag
        app_data: Selected colormap
        user_data: User data (unused)
    """
    app_state.current_colormap = app_data
    logger.debug(f"Colormap changed to: {app_data}")


def change_variable_callback(sender, app_data, user_data):
    """
    Callback when plot variable is changed
    
    Args:
        sender: Combo box tag
        app_data: Selected variable
        user_data: User data (unused)
    """
    app_state.current_variable = app_data
    logger.debug(f"Plot variable changed to: {app_data}")


def open_plot_settings_callback(sender, app_data, user_data):
    """
    Open a dialog with advanced plot settings
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # TODO: Implement advanced plot settings dialog
    logger.info("Advanced plot settings not yet implemented")


def close_all_plot_windows_callback(sender, app_data, user_data):
    """
    Close all external plot windows
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    closed = False
    
    if app_state.plotter_process is not None:
        try:
            if app_state.plotter_process.poll() is None:
                app_state.plotter_process.terminate()
                app_state.plotter_process.wait(timeout=2)
                logger.info("Plotter window closed")
                closed = True
        except Exception as e:
            logger.warning(f"Error closing plotter: {e}")
        finally:
            app_state.plotter_process = None
    
    if not closed:
        logger.info("No plot windows to close")


def browse_vtk_file_callback(sender, app_data, user_data):
    """
    Open file dialog for VTK file selection
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if dpg.does_item_exist("file_dialog_vtk"):
        dpg.configure_item("file_dialog_vtk", show=True)


def browse_output_folder_callback(sender, app_data, user_data):
    """
    Open the output folder in file explorer
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils import open_folder
    from src.utils.output_manager import get_latest_output_folder
    from src.config import OUTPUT_DIR
    
    # Try to open latest output folder
    latest_folder = get_latest_output_folder()
    
    if latest_folder and os.path.exists(latest_folder):
        success = open_folder(latest_folder)
        if success:
            logger.info(f"Opened latest output folder: {latest_folder}")
    else:
        # Fall back to main output directory
        if os.path.exists(OUTPUT_DIR):
            success = open_folder(OUTPUT_DIR)
            if success:
                logger.info(f"Opened output directory: {OUTPUT_DIR}")
        else:
            logger.warning("No output folder found yet. Run the solver first.")


def set_vtk_from_latest_output_callback(sender, app_data, user_data):
    """
    Automatically set VTK path from latest output folder
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils.output_manager import get_latest_output_folder
    from src.config import DEFAULT_VTK
    
    latest_folder = get_latest_output_folder()
    
    if not latest_folder:
        logger.warning("No output folder found. Run the solver first.")
        return
    
    # Look for VTK files in latest output
    vtk_files = []
    for file in os.listdir(latest_folder):
        if file.endswith('.vtk'):
            vtk_files.append(file)
    
    if not vtk_files:
        logger.warning(f"No VTK files found in {latest_folder}")
        return
    
    # Use the first VTK file found (or Solution.vtk if it exists)
    if DEFAULT_VTK in vtk_files:
        vtk_file = DEFAULT_VTK
    else:
        vtk_file = vtk_files[0]
    
    vtk_path = os.path.join(latest_folder, vtk_file)
    
    if dpg.does_item_exist("contour_vtk_path"):
        dpg.set_value("contour_vtk_path", vtk_path)
        logger.success(f"VTK path set to: {vtk_path}")
    else:
        logger.error("Could not set VTK path widget")