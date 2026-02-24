"""
Main window construction for MeMPhyS GUI

Assembles all UI components into the main application window.
"""

import dearpygui.dearpygui as dpg

from src.config import APP_FULL_NAME
from .menu_bar import create_menu_bar
from .parameters_panel import create_parameters_panel
from .visualization_panel import create_visualization_panel


def create_main_window(themes: dict) -> int:
    """
    Create the main application window
    Args: themes: Dictionary of theme IDs
    Returns: Main window tag/ID
    """
    with dpg.window(label=APP_FULL_NAME,
        tag="MainWindow", no_close=False,  # Allow closing
        no_collapse=True) as main_window:
        # Menu Bar
        create_menu_bar(themes)
        # Main Content Area (horizontal split)
        with dpg.group(horizontal=True):
            # Left Panel: Parameters and Controls
            create_parameters_panel(themes)
            # Right Panel: Visualization, BC, and Logs
            with dpg.child_window(width=-1, border=True):
                # Visualization panel (plots and controls)
                create_visualization_panel(themes)
                dpg.add_spacer(height=10)                    
    return main_window