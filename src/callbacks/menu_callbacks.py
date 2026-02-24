"""
Menu bar callbacks for MeMPhyS GUI

Handles callbacks for:
- File menu (open, save, exit)
- Edit menu (preferences)
- Help menu (help, about)
"""

import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils import open_folder, open_url, change_font, get_platform
from src.config import (
    LOG_DIR,
    HELP_URL,
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    FONT_SIZES,
    FONT_PREFERENCES,
)
from src.config import COLORS, COLORS_DARK, COLORS_LIGHT, _THEMED_TEXTS, themed_texts, update_themed_texts, toggle_theme


def open_logs_callback(sender, app_data, user_data):
    """
    Open the logs directory
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    success = open_folder(LOG_DIR)
    
    if not success:
        logger.error(f"Failed to open logs directory: {LOG_DIR}")


def open_help_callback(sender, app_data, user_data):
    """
    Open the help documentation in web browser
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    success = open_url(HELP_URL)
    
    if success:
        logger.info("Opened help documentation in browser")
    else:
        logger.error("Failed to open help URL")


def show_about_callback(sender, app_data, user_data):
    """
    Show the About dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # If window already exists, just show it
    if dpg.does_item_exist("about_window"):
        dpg.configure_item("about_window", show=True)
        dpg.focus_item("about_window")
        return
    
    # Create About window
    with dpg.window(
        label="About",
        tag="about_window",
        modal=True,
        width=620,
        height=400,
        no_resize=True,
        pos=(330, 200)
    ):
        dpg.add_text(APP_FULL_NAME, color=(200, 220, 255))
        dpg.add_text(APP_SUBTITLE, color=(200, 220, 255))
        dpg.add_separator()
        
        dpg.add_text("GUI Development:", color=(255, 220, 160))
        for dev in DEVELOPERS["gui"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        dpg.add_text("Backend / Simulation:", color=(255, 220, 160))
        for dev in DEVELOPERS["backend"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        dpg.add_text("Affiliation:", color=(255, 220, 160))
        for affil in DEVELOPERS["affiliation"]:
            dpg.add_text(f"\t{affil}")
        
        dpg.add_spacer(height=12)
        
        close_btn = dpg.add_button(
            label="Close",
            width=80,
            callback=lambda: dpg.configure_item("about_window", show=False)
        )


def show_preferences_callback(sender, app_data, user_data):
    """
    Show the Preferences dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # If window already exists, just show it
    if dpg.does_item_exist("preferences_window"):
        dpg.configure_item("preferences_window", show=True)
        dpg.focus_item("preferences_window")
        return
    
    # Get current font settings
    current_font_name = app_state.current_font_name
    current_font_size = app_state.current_font_size
    
    # Get available fonts for current platform
    platform = get_platform()
    available_fonts = FONT_PREFERENCES.get(platform, [])
    
    # Create Preferences window
    with dpg.window(
        label="Preferences",
        tag="preferences_window",
        modal=True,
        width=450,
        height=250,
        no_resize=True,
        pos=(415, 275)
    ):
        dpg.add_text("Preferences", color=(200, 220, 255))
        dpg.add_separator()
        
        dpg.add_combo(
            label="Font family",
            items=available_fonts,
            default_value=current_font_name,
            tag="pref_font_family",
            width=200
        )
        
        dpg.add_combo(
            label="Font size",
            items=FONT_SIZES,
            default_value=current_font_size,
            tag="pref_font_size",
            width=200
        )
        
        dpg.add_spacer(height=12)
        
        dpg.add_button(
            label="Apply",
            width=80,
            callback=apply_preferences_callback
        )
        
        dpg.add_same_line()
        
        dpg.add_button(
            label="Close",
            width=80,
            callback=lambda: dpg.configure_item("preferences_window", show=False)
        )


def apply_preferences_callback(sender, app_data, user_data):
    """
    Apply preferences changes
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if not dpg.does_item_exist("pref_font_family") or not dpg.does_item_exist("pref_font_size"):
        logger.error("Preference widgets not found")
        return
    
    font_name = dpg.get_value("pref_font_family")
    font_size = dpg.get_value("pref_font_size")
    
    # Validate font name is not empty
    if not font_name or font_name.strip() == "":
        logger.error("Invalid font name. Please select a valid font.")
        return
    
    logger.info(f"Attempting to change font to {font_name} @ {font_size}px")
    
    # Update app state
    app_state.current_font_name = font_name
    app_state.current_font_size = font_size
    
    # Save preferences
    from src.utils.config_manager import save_user_preferences
    save_user_preferences()
    
    # Try to change font
    success = change_font(font_name, font_size)
    
    if success:
        logger.success(f"Font preference saved: {font_name} @ {font_size}px")
        logger.info("Changes will apply fully on next restart.")
    else:
        logger.warning("Could not apply font immediately")
        logger.info("Font preference saved. Will be applied on next restart.")


def show_options_callback(sender, app_data, user_data):
    """
    Show the Options dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.ui.dialogs import create_options_dialog
    create_options_dialog()


# def toggle_theme_callback(sender, app_data, user_data):
#     """
#     Toggle between dark and light themes
    
#     Args:
#         sender: Menu item tag
#         app_data: Application data (unused)
#         user_data: User data (unused)
#     """
#     from src.config.themes import toggle_theme
    
#     # Get current theme state
#     current_dark = app_state.get_option("dark_theme", True)
    
#     # Toggle
#     new_dark = toggle_theme(current_dark)
    
#     # Save new state
#     app_state.set_option("dark_theme", new_dark)
    
#     theme_name = "Dark" if new_dark else "Light"
#     logger.success(f"Switched to {theme_name} theme")
#     logger.info("Theme preference will be saved on exit")


def toggle_theme_callback():
    # global _is_dark, COLORS
    # _is_dark = not _is_dark
    current_dark = app_state.get_option("dark_theme", True)
    
    new_palette = COLORS_DARK if current_dark else COLORS_LIGHT

    # Update the live COLORS dict in-place so any future widgets use it
    COLORS.update(new_palette)
    update_themed_texts(new_palette)
    # Re-color every tracked widget instantly
    for tag, color_key in _THEMED_TEXTS.items():
        if dpg.does_item_exist(tag):
            dpg.configure_item(tag, color=new_palette[color_key])

    # Swap the global DPG theme as before
    new_dark = toggle_theme(current_dark)
    
    # Save new state
    app_state.set_option("dark_theme", new_dark)
    
    theme_name = "Dark" if new_dark else "Light"
    logger.success(f"Switched to {theme_name} theme")
    logger.info("Theme preference will be saved on exit")


def exit_application_callback(sender, app_data, user_data):
    """
    Exit the application with cleanup
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Exiting application...")
    
    # Auto-save configuration
    from src.utils.config_manager import auto_save_on_exit
    auto_save_on_exit()
    
    # Perform cleanup
    app_state.cleanup()
    
    # Stop DearPyGUI
    dpg.stop_dearpygui()


def open_config_callback(sender, app_data, user_data):
    """
    Open a configuration file (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Configuration loading not yet implemented")
    
    if dpg.does_item_exist("file_dialog_config"):
        dpg.configure_item("file_dialog_config", show=True)


def save_config_callback(sender, app_data, user_data):
    """
    Save current configuration (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Configuration saving not yet implemented")
    
    if dpg.does_item_exist("file_dialog_config_save"):
        dpg.configure_item("file_dialog_config_save", show=True)


def show_system_info_callback(sender, app_data, user_data):
    """
    Show system information dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils import get_system_info_summary
    
    info = get_system_info_summary()
    
    logger.separator()
    logger.info("System Information:")
    for line in info.split('\n'):
        if line.strip():
            logger.info(line)
    logger.separator()


def show_application_state_callback(sender, app_data, user_data):
    """
    Show current application state (for debugging)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    state_summary = app_state.get_state_summary()
    
    logger.separator()
    logger.info("Application State:")
    for key, value in state_summary.items():
        logger.info(f"  {key}: {value}")
    logger.separator()


def check_for_updates_callback(sender, app_data, user_data):
    """
    Check for application updates (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Update checking not yet implemented")
    logger.info(f"Current version: {APP_FULL_NAME}")


def report_bug_callback(sender, app_data, user_data):
    """
    Open bug report page
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    bug_url = f"{HELP_URL}/issues"
    success = open_url(bug_url)
    
    if success:
        logger.info("Opened bug report page in browser")
    else:
        logger.error("Failed to open bug report URL")
        logger.info(f"Please visit: {bug_url}")


def save_config_file_callback(sender, app_data, user_data):
    """
    Callback when user selects save location
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: User data (unused)
    """
    if 'file_path_name' not in app_data:
        logger.error("No file path in save dialog")
        return
    
    file_path = app_data['file_path_name']
    
    # Ensure .json extension
    if not file_path.endswith('.json'):
        file_path += '.json'
    
    from src.utils.config_manager import save_configuration
    success = save_configuration(file_path)
    
    if success:
        logger.success(f"Configuration saved to {file_path}")
    else:
        logger.error(f"Failed to save configuration to {file_path}")


def load_config_file_callback(sender, app_data, user_data):
    """
    Callback when user selects config file to load
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: User data (unused)
    """
    if 'file_path_name' not in app_data:
        logger.error("No file path in load dialog")
        return
    
    file_path = app_data['file_path_name']
    
    from src.utils.config_manager import load_configuration
    success = load_configuration(file_path)
    
    if success:
        logger.success(f"Configuration loaded from {file_path}")
        logger.info("Configuration applied. Some changes (like fonts) require restart.")
        
        # Trigger UI updates for loaded settings
        from src.callbacks import show_multigrid_callback, show_implicit_callback
        
        # Update multigrid UI
        if dpg.does_item_exist("multigrid_toggle"):
            show_multigrid_callback(None, None, None)
        
        # Update implicit parameters UI
        if dpg.does_item_exist("solver_method"):
            show_implicit_callback(None, None, None)
        
        logger.info("UI updated with loaded configuration")
    else:
        logger.error(f"Failed to load configuration from {file_path}")

