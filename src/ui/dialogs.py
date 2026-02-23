"""
Dialog windows for MeMPhyS GUI

Creates modal dialogs for About and Preferences.
Note: These are created on-demand by callbacks, not at startup.
"""

import dearpygui.dearpygui as dpg

from src.config import (
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    FONT_SIZES,
    FONT_PREFERENCES,
    COLORS,
    BC_TYPES, 
    BC_VARIABLES,
)
from src.utils import get_platform
from src.core import app_state
from src.callbacks import apply_preferences_callback

from src.utils.gmsh_bc_manager import (
    get_physical_names,
    get_boundary_condition,
    set_boundary_condition,
)


def create_about_dialog():
    """
    Create About dialog window
    This function is called by the menu callback when needed.
    """
    # Check if already exists
    if dpg.does_item_exist("about_window"):
        dpg.configure_item("about_window", show=True)
        dpg.focus_item("about_window")
        return
    
    with dpg.window(label="About", tag="about_window",
        modal=True, width=620, height=400,
        no_resize=True, pos=(330, 200)):
        
        dpg.add_text(APP_FULL_NAME, color=COLORS["header"])
        dpg.add_text(APP_SUBTITLE, color=COLORS["header"])
        dpg.add_separator()
        
        dpg.add_text("GUI Development:", color=COLORS["subheader"])
        for dev in DEVELOPERS["gui"]:
            dpg.add_text(f"\t{dev}")
        dpg.add_spacer(height=6)
        
        dpg.add_text("Backend / Simulation:", color=COLORS["subheader"])
        for dev in DEVELOPERS["backend"]:
            dpg.add_text(f"\t{dev}")
        dpg.add_spacer(height=6)
        
        dpg.add_text("Affiliation:", color=COLORS["subheader"])
        for affil in DEVELOPERS["affiliation"]:
            dpg.add_text(f"\t{affil}")
        dpg.add_spacer(height=12)
        
        dpg.add_button(label="Close", width=80,
            callback=lambda: dpg.configure_item("about_window", show=False))


def create_preferences_dialog():
    """
    Create Preferences dialog window    
    This function is called by the menu callback when needed.
    """
    if dpg.does_item_exist("preferences_window"):
        dpg.configure_item("preferences_window", show=True)
        dpg.focus_item("preferences_window")
        return
    
    current_font_name = app_state.current_font_name
    current_font_size = app_state.current_font_size

    platform = get_platform()
    available_fonts = FONT_PREFERENCES.get(platform, [])

    if not current_font_name or current_font_name not in available_fonts:
        current_font_name = available_fonts[0] if available_fonts else "DejaVuSans.ttf"
    
    with dpg.window(label="Preferences", tag="preferences_window",
        modal=True, width=450, height=250,
        no_resize=True, pos=(415, 275)):
        dpg.add_text("Preferences", color=COLORS["header"])
        dpg.add_separator()
        
        # Font family selector
        dpg.add_combo(label="Font family", items=available_fonts,
            default_value=current_font_name,
            tag="pref_font_family", width=200)
        
        # Font size selector
        dpg.add_combo(label="Font size", items=FONT_SIZES,
            default_value=current_font_size, tag="pref_font_size",
            width=200)
        
        dpg.add_spacer(height=12)
        
        with dpg.group(horizontal=True):
            dpg.add_button(label="Apply", width=80,
                callback=apply_preferences_callback)
            dpg.add_button(label="Close", width=80,
                callback=lambda: dpg.configure_item("preferences_window", show=False))

def show_error_dialog(title: str, message: str):
    """
    Show an error dialog    
    Args:
        title: Dialog title
        message: Error message to display
    """
    dialog_tag = f"error_dialog_{id(message)}"
    
    with dpg.window(label=title, tag=dialog_tag,
        modal=True, width=400, height=200,
        no_resize=True, pos=(440, 300)):
        dpg.add_text("Error:", color=COLORS["error"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        dpg.add_button(label="OK", width=80,
            callback=lambda: dpg.delete_item(dialog_tag))

def show_info_dialog(title: str, message: str):
    """
    Show an information dialog    
    Args:
        title: Dialog title
        message: Information message to display
    """
    dialog_tag = f"info_dialog_{id(message)}"
    
    with dpg.window(label=title, tag=dialog_tag,
        modal=True, width=400, height=200,
        no_resize=True,pos=(440, 300)):
        dpg.add_text("Information:", color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        dpg.add_button(label="OK",width=80,
            callback=lambda: dpg.delete_item(dialog_tag))

def show_confirmation_dialog(title: str, message: str, on_confirm, on_cancel=None):
    """
    Show a confirmation dialog with Yes/No buttons    
    Args:
        title: Dialog title
        message: Confirmation message
        on_confirm: Callback function to call on Yes
        on_cancel: Optional callback function to call on No
    """
    dialog_tag = f"confirm_dialog_{id(message)}"
    
    def confirm_callback():
        dpg.delete_item(dialog_tag)
        if on_confirm:
            on_confirm()
    
    def cancel_callback():
        dpg.delete_item(dialog_tag)
        if on_cancel:
            on_cancel()
    
    with dpg.window(label=title, tag=dialog_tag,
        modal=True, width=400, height=200,
        no_resize=True, pos=(440, 300)):
        dpg.add_text("Confirm:", color=COLORS["warning"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        with dpg.group(horizontal=True):
            dpg.add_button(label="Yes", width=80,
                callback=confirm_callback)
            dpg.add_button(label="No", width=80,
                callback=cancel_callback)


def create_options_dialog():
    """
    Create Options dialog window
    
    This function is called by the menu callback when needed.
    """
    from src.config import DEFAULT_OPTIONS
    
    # Check if already exists
    if dpg.does_item_exist("options_window"):
        dpg.configure_item("options_window", show=True)
        dpg.focus_item("options_window")
        return
    
    with dpg.window(
        label="Options",
        tag="options_window",
        modal=True,
        width=600,
        height=600,
        no_resize=True,
        pos=(340, 100)
    ):
        dpg.add_text("Application Options", color=COLORS["header"])
        dpg.add_separator()
        dpg.add_text("Configure application behavior and features", color=COLORS["info"])
        dpg.add_spacer(height=10)
        
        # Scrollable area for options
        with dpg.child_window(height=450, border=True):
            # Logging Options
            dpg.add_text("Logging Options:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Enable GUI Logging",
                tag="opt_enable_gui_logging",
                default_value=app_state.get_option("enable_gui_logging", DEFAULT_OPTIONS["enable_gui_logging"])
            )
            with dpg.tooltip("opt_enable_gui_logging"):
                dpg.add_text("Show log messages in the GUI log window")
            
            dpg.add_checkbox(
                label="Enable File Logging",
                tag="opt_enable_file_logging",
                default_value=app_state.get_option("enable_file_logging", DEFAULT_OPTIONS["enable_file_logging"])
            )
            with dpg.tooltip("opt_enable_file_logging"):
                dpg.add_text("Write log messages to log files")
            
            dpg.add_checkbox(
                label="Enable Console Logging",
                tag="opt_enable_console_logging",
                default_value=app_state.get_option("enable_console_logging", DEFAULT_OPTIONS["enable_console_logging"])
            )
            with dpg.tooltip("opt_enable_console_logging"):
                dpg.add_text("Print log messages to terminal console")
            
            dpg.add_spacer(height=15)
            dpg.add_separator()
            dpg.add_spacer(height=15)
            
            # File Management Options
            dpg.add_text("File Management:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Auto-move Output Files",
                tag="opt_auto_move_outputs",
                default_value=app_state.get_option("auto_move_outputs", DEFAULT_OPTIONS["auto_move_outputs"])
            )
            with dpg.tooltip("opt_auto_move_outputs"):
                dpg.add_text("Automatically move solver outputs to organized folders after completion")
            
            dpg.add_checkbox(
                label="Create Dated Folders",
                tag="opt_create_dated_folders",
                default_value=app_state.get_option("create_dated_folders", DEFAULT_OPTIONS["create_dated_folders"])
            )
            with dpg.tooltip("opt_create_dated_folders"):
                dpg.add_text("Organize outputs into folders named by date (YYYY-MM-DD)")
            
            dpg.add_checkbox(
                label="Append Number if File Exists",
                tag="opt_append_number_if_exists",
                default_value=app_state.get_option("append_number_if_exists", DEFAULT_OPTIONS["append_number_if_exists"])
            )
            with dpg.tooltip("opt_append_number_if_exists"):
                dpg.add_text("Add _1, _2, _3 etc. to filenames to avoid overwriting existing files")
            
            dpg.add_spacer(height=15)
            dpg.add_separator()
            dpg.add_spacer(height=15)
            
            # Monitoring Options
            dpg.add_text("Monitoring:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Auto-start Convergence Monitor",
                tag="opt_auto_start_convergence_monitor",
                default_value=app_state.get_option("auto_start_convergence_monitor", DEFAULT_OPTIONS["auto_start_convergence_monitor"])
            )
            with dpg.tooltip("opt_auto_start_convergence_monitor"):
                dpg.add_text("Automatically start monitoring convergence data")
            
            dpg.add_checkbox(
                label="Convergence Plot Auto-scale",
                tag="opt_convergence_plot_auto_scale",
                default_value=app_state.get_option("convergence_plot_auto_scale", DEFAULT_OPTIONS["convergence_plot_auto_scale"])
            )
            with dpg.tooltip("opt_convergence_plot_auto_scale"):
                dpg.add_text("Automatically adjust plot axes based on data")
            
            dpg.add_spacer(height=15)
            dpg.add_separator()
            dpg.add_spacer(height=15)
            
            # Validation Options
            dpg.add_text("Validation:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Validate Parameters on Run",
                tag="opt_validate_parameters_on_run",
                default_value=app_state.get_option("validate_parameters_on_run", DEFAULT_OPTIONS["validate_parameters_on_run"])
            )
            with dpg.tooltip("opt_validate_parameters_on_run"):
                dpg.add_text("Check all parameters are valid before running solver")
            
            dpg.add_checkbox(
                label="Validate Mesh Files on Run",
                tag="opt_validate_mesh_files_on_run",
                default_value=app_state.get_option("validate_mesh_files_on_run", DEFAULT_OPTIONS["validate_mesh_files_on_run"])
            )
            with dpg.tooltip("opt_validate_mesh_files_on_run"):
                dpg.add_text("Verify mesh files exist and are readable before running")
            
            dpg.add_spacer(height=15)
            dpg.add_separator()
            dpg.add_spacer(height=15)
            
            # UI Options
            dpg.add_text("User Interface:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Show Tooltips",
                tag="opt_show_tooltips",
                default_value=app_state.get_option("show_tooltips", DEFAULT_OPTIONS["show_tooltips"])
            )
            with dpg.tooltip("opt_show_tooltips"):
                dpg.add_text("Display helpful tooltips when hovering over UI elements")
            
            dpg.add_checkbox(
                label="Confirm on Exit",
                tag="opt_confirm_on_exit",
                default_value=app_state.get_option("confirm_on_exit", DEFAULT_OPTIONS["confirm_on_exit"])
            )
            with dpg.tooltip("opt_confirm_on_exit"):
                dpg.add_text("Ask for confirmation before closing the application")
            
            dpg.add_spacer(height=15)
            dpg.add_separator()
            dpg.add_spacer(height=15)
            
            # Solver Options
            dpg.add_text("Solver:", color=COLORS["subheader"])
            dpg.add_spacer(height=5)
            
            dpg.add_checkbox(
                label="Auto-open Plot on Complete",
                tag="opt_auto_open_plot_on_complete",
                default_value=app_state.get_option("auto_open_plot_on_complete", DEFAULT_OPTIONS["auto_open_plot_on_complete"])
            )
            with dpg.tooltip("opt_auto_open_plot_on_complete"):
                dpg.add_text("Automatically open visualization when solver completes")
            
            dpg.add_checkbox(
                label="Play Sound on Complete",
                tag="opt_play_sound_on_complete",
                default_value=app_state.get_option("play_sound_on_complete", DEFAULT_OPTIONS["play_sound_on_complete"])
            )
            with dpg.tooltip("opt_play_sound_on_complete"):
                dpg.add_text("Play notification sound when solver finishes")
        
        dpg.add_spacer(height=12)
        
        # Buttons
        with dpg.group(horizontal=True):
            dpg.add_button(
                label="Apply",
                width=100,
                callback=lambda: apply_options_callback()
            )
            
            dpg.add_button(
                label="Reset to Defaults",
                width=150,
                callback=lambda: reset_options_to_defaults()
            )
            
            dpg.add_button(
                label="Close",
                width=100,
                callback=lambda: dpg.configure_item("options_window", show=False)
            )


def apply_options_callback():
    """Apply options from the dialog"""
    from src.core import logger, app_state
    from src.config import DEFAULT_OPTIONS
    from src.utils.config_manager import save_app_options
    
    # Read all checkbox values and save to app_state
    for option_key in DEFAULT_OPTIONS.keys():
        tag = f"opt_{option_key}"
        if dpg.does_item_exist(tag):
            value = dpg.get_value(tag)
            app_state.set_option(option_key, value)
    
    # Save options to file
    save_app_options()
    
    # Apply logging options immediately
    from src.core import logger
    logger.set_enable_gui(app_state.get_option("enable_gui_logging", True))
    logger.set_enable_file(app_state.get_option("enable_file_logging", True))
    logger.set_enable_console(app_state.get_option("enable_console_logging", True))
    
    logger.success("Options applied and saved successfully")


def reset_options_to_defaults():
    """Reset all options to default values"""
    from src.core import logger, app_state
    from src.config import DEFAULT_OPTIONS
    
    # Reset app_state options
    app_state.load_default_options()
    
    # Update UI checkboxes
    for option_key, default_value in DEFAULT_OPTIONS.items():
        tag = f"opt_{option_key}"
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, default_value)
    
    logger.info("Options reset to defaults")