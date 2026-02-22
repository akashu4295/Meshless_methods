"""
MeMPhyS GUI

Meshless Multi-Physics Solver
Version 2.5.1

This is the main entry point for the MeMPhyS GUI application.
"""

import sys
import dearpygui.dearpygui as dpg

# Import configuration
from src.config import (
    MAIN_WINDOW_WIDTH,
    MAIN_WINDOW_HEIGHT,
    APP_FULL_NAME,
    LOG_FILE_PATH,
    DEFAULT_FONT_SIZE,
)
from src.config.themes import initialize_all_themes, apply_dark_theme
from src.core import app_state, logger  # Import core modules
from src.utils import initialize_fonts, set_default_font    # Import utilities
from src.solver import convergence_monitor  # Import solver modules
from src.ui import create_main_window   # Import UI

def initialize_application():
    """
    Initialize the application
    
    Returns:
        Dictionary of themes
    """
    print("=== Starting initialization ===")
    
    # Disable GUI logging during initialization (no GUI exists yet)
    logger.set_enable_gui(False)
    
    print("Step 1: Opening log file...")
    logger.info(f"Starting {APP_FULL_NAME}")
    logger.separator()
    
    # Open log file
    log_handle = app_state.open_log_file(LOG_FILE_PATH, mode='a')
    if log_handle:
        logger.set_file_handle(log_handle)
        logger.info("Log file opened successfully")
        print("Log file opened")
    else:
        logger.warning("Could not open log file")
        print("Warning: Could not open log file")
    
    # Create DearPyGUI context
    print("Step 2: Creating DearPyGUI context...")
    dpg.create_context()
    logger.info("DearPyGUI context created")
    print("DearPyGUI context created")
    
    # Load default options
    print("Step 2b: Loading default options...")
    app_state.load_default_options()
    
    # Try to load saved options
    from src.utils.config_manager import load_app_options, apply_app_options
    saved_options = load_app_options()
    if saved_options:
        apply_app_options(saved_options)
    
    logger.info("Application options loaded")
    print("Options loaded")
    
    # Load user preferences
    print("Step 2c: Loading user preferences...")
    from src.utils.config_manager import load_user_preferences, apply_user_preferences
    saved_prefs = load_user_preferences()
    if saved_prefs:
        apply_user_preferences(saved_prefs)
    
    logger.info("User preferences loaded")
    print("User preferences loaded")
    
    # Initialize fonts (optional - can be disabled if causing issues)
    ENABLE_CUSTOM_FONTS = False  # Set to True once font issues are resolved
    
    if ENABLE_CUSTOM_FONTS:
        try:
            logger.info("Initializing fonts...")
            font_registry = initialize_fonts()
            logger.success(f"Loaded {len(font_registry)} font variants")
            
            # Set default font
            # Try to use the first available font
            if font_registry:
                first_font = list(font_registry.keys())[0]
                font_name, size = first_font
                
                # Try to set the default size
                try:
                    if set_default_font(font_name, DEFAULT_FONT_SIZE):
                        logger.success(f"Set default font: {font_name} @ {DEFAULT_FONT_SIZE}px")
                    else:
                        # Fall back to first available
                        set_default_font(font_name, size)
                        logger.info(f"Using font: {font_name} @ {size}px")
                except Exception as font_error:
                    logger.warning(f"Could not set default font: {font_error}")
                    logger.info("Continuing without custom font")
        
        except Exception as e:
            logger.log_exception(e, "Error initializing fonts")
            logger.warning("Continuing with default fonts")
    else:
        logger.info("Custom fonts disabled, using DearPyGUI defaults")
    
    # Initialize themes
    print("Step 3: Initializing themes...")
    themes = initialize_all_themes()
    print("Themes initialized")
    
    print("Step 4: Applying dark theme...")
    apply_dark_theme()
    logger.success("Themes initialized")
    print("Dark theme applied")
    
    # Mark application as initialized
    app_state.mark_initialized()
    print("=== Initialization complete ===")
    
    return themes


def create_gui(themes: dict):
    """
    Create the GUI
    
    Args:
        themes: Dictionary of theme IDs
    """
    print("Step A: Starting GUI creation...")
    logger.info("Creating GUI components...")
    
    print("Step B: Creating main window...")
    # Create main window
    main_window = create_main_window(themes)
    print("Step C: Main window created")
    logger.success("Main window created")
    
    print("Step D: Setting primary window...")
    # Set as primary window (fills entire viewport)
    dpg.set_primary_window("MainWindow", True)
    print("Step E: Primary window set")
    
    # Re-enable GUI logging now that log window exists
    logger.set_enable_gui(True)
    logger.info("GUI logging enabled")
    print("Step F: GUI creation complete")
    return main_window


def setup_viewport():
    """Setup the application viewport"""
    dpg.create_viewport(
        title=APP_FULL_NAME,
        width=MAIN_WINDOW_WIDTH,
        height=MAIN_WINDOW_HEIGHT,
        resizable=True
    )
    logger.info(f"Viewport created: {MAIN_WINDOW_WIDTH}x{MAIN_WINDOW_HEIGHT}")


def start_background_services():
    """Start background services (convergence monitor, etc.)"""
    logger.info("Starting background services...")
    
    # Start convergence monitor
    try:
        convergence_monitor.start()
        logger.success("Convergence monitor started")
    except Exception as e:
        logger.warning(f"Could not start convergence monitor: {e}")
        logger.info("Convergence monitoring will be unavailable")


def cleanup_and_exit():
    """Cleanup resources and exit"""
    print("Starting cleanup...")
    logger.info("Shutting down application...")
    logger.separator()
    
    try:
        # Stop convergence monitor
        print("Stopping convergence monitor...")
        convergence_monitor.cleanup()
    except Exception as e:
        print(f"Warning: Error during convergence monitor cleanup: {e}")
    
    try:
        # Cleanup application state
        print("Cleaning up application state...")
        app_state.cleanup()
    except Exception as e:
        print(f"Warning: Error during app state cleanup: {e}")
    
    try:
        # Destroy DearPyGUI context
        print("Destroying DearPyGUI context...")
        dpg.destroy_context()
    except Exception as e:
        print(f"Warning: Error destroying DearPyGUI context: {e}")
    
    print("Cleanup complete")
    logger.info("Application closed successfully")


def main():
    """Main application entry point"""
    try:
        print("=== MeMPhyS GUI Starting ===")
        
        # Initialize application
        print("\n--- Initialization Phase ---")
        themes = initialize_application()
        
        # Create GUI
        print("\n--- GUI Creation Phase ---")
        create_gui(themes)
        
        # Setup viewport
        print("\n--- Viewport Setup Phase ---")
        setup_viewport()
        
        # Setup DearPyGUI
        print("\n--- DearPyGUI Setup Phase ---")
        dpg.setup_dearpygui()
        print("DearPyGUI setup complete")
        
        # Show viewport
        print("\n--- Showing Viewport ---")
        dpg.show_viewport()
        print("Viewport shown")
        
        logger.success("Application initialized successfully")
        logger.separator()
        logger.info("Ready for input")
        print("\n=== Application Ready ===\n")
        
        # Main render loop (this blocks until window closes)
        print("Starting main render loop...")
        
        # Start convergence monitor in a frame callback (after GUI is running)
        def delayed_start():
            try:
                convergence_monitor.start()
                logger.success("Convergence monitor started")
            except Exception as e:
                logger.warning(f"Could not start convergence monitor: {e}")
            
            # Restore last session
            from src.utils.config_manager import load_session_state, restore_session_state
            session = load_session_state()
            if session:
                restore_session_state(session)
        
        dpg.set_frame_callback(2, delayed_start)
        
        dpg.start_dearpygui()
        print("Main loop ended")
        
    except Exception as e:
        logger.log_exception(e, "Critical error during startup")
        print(f"Critical error: {e}")
        sys.exit(1)
    
    finally:
        # Cleanup
        cleanup_and_exit()


if __name__ == "__main__":
    main()