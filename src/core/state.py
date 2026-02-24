"""
Application state management for MeMPhyS GUI

This module manages the global application state including:
- Active processes (solver, plotter)
- File handles (log files)
- Font registry
- UI state (current selections, running status)
"""

import os
from typing import Optional, Dict, Any
from pathlib import Path


class AppState:
    """
    Centralized application state manager
    
    This class encapsulates all mutable state in the application,
    providing a single source of truth for:
    - Running processes
    - File handles
    - Font configurations
    - UI state
    """
    
    def __init__(self):
        # Process management
        self._plotter_process: Optional[Any] = None
        self._solver_process: Optional[Any] = None
        self._convergence_monitor_running: bool = False
        
        # File handles
        self._log_file_handle: Optional[Any] = None
        self._log_file_path: Optional[str] = None
        
        # Font management
        self._font_registry: Dict[tuple, int] = {}  # (name, size) -> font_id
        self._current_font_name: str = "DejaVuSans.ttf"  # Default font
        self._current_font_size: int = 16
        
        # UI state
        self._solver_running: bool = False
        self._multigrid_enabled: bool = False
        self._num_mesh_levels: int = 1
        self._current_solver_method: str = "Fractional Step"
        
        # Mesh file paths
        self._mesh_files: Dict[int, str] = {}  # level -> path
        self._init_file_path: str = ""
        
        # Visualization state
        self._current_vtk_file: str = ""
        self._current_colormap: str = "viridis"
        self._current_variable: str = "velocity magnitude"
        
        # Flags
        self._gui_initialized: bool = False
        self._cleanup_done: bool = False
        
        # User options
        self._options: Dict[str, bool] = {}

        # Paraview
        self.paraview_path = None
    
    # ==================== Options Management ====================
    
    def set_option(self, key: str, value: bool):
        """Set an option value"""
        self._options[key] = value
    
    def get_option(self, key: str, default: bool = False) -> bool:
        """Get an option value"""
        return self._options.get(key, default)
    
    def get_all_options(self) -> Dict[str, bool]:
        """Get all options"""
        return self._options.copy()
    
    def load_default_options(self):
        """Load default options"""
        from src.config import DEFAULT_OPTIONS
        self._options = DEFAULT_OPTIONS.copy()
    
    # ==================== Process Management ====================
    
    @property
    def plotter_process(self):
        """Get the current plotter process"""
        return self._plotter_process
    
    @plotter_process.setter
    def plotter_process(self, process):
        """Set the plotter process"""
        self._plotter_process = process
    
    @property
    def solver_process(self):
        """Get the current solver process"""
        return self._solver_process
    
    @solver_process.setter
    def solver_process(self, process):
        """Set the solver process"""
        self._solver_process = process
    
    @property
    def convergence_monitor_running(self) -> bool:
        """Check if convergence monitor thread is running"""
        return self._convergence_monitor_running
    
    @convergence_monitor_running.setter
    def convergence_monitor_running(self, value: bool):
        """Set convergence monitor running state"""
        self._convergence_monitor_running = value
    
    def terminate_plotter_process(self):
        """Safely terminate the plotter process if running"""
        if self._plotter_process is not None:
            try:
                if self._plotter_process.poll() is None:
                    self._plotter_process.terminate()
                    self._plotter_process.wait(timeout=2)
                    print("Plotter process terminated successfully")
            except Exception as e:
                print(f"Error terminating plotter process: {e}")
            finally:
                self._plotter_process = None
    
    def terminate_solver_process(self):
        """Safely terminate the solver process if running"""
        if self._solver_process is not None:
            try:
                if self._solver_process.poll() is None:
                    self._solver_process.terminate()
                    self._solver_process.wait(timeout=2)
                    print("Solver process terminated successfully")
            except Exception as e:
                print(f"Error terminating solver process: {e}")
            finally:
                self._solver_process = None
    
    # ==================== File Handle Management ====================
    
    @property
    def log_file_handle(self):
        """Get the log file handle"""
        return self._log_file_handle
    
    def open_log_file(self, file_path: str, mode: str = 'a'):
        """
        Open a log file for writing
        
        Args:
            file_path: Path to the log file
            mode: File opening mode ('a' for append, 'w' for write)
        
        Returns:
            File handle or None if failed
        """
        try:
            # Close existing file if open
            self.close_log_file()
            
            # Open new file
            self._log_file_handle = open(file_path, mode)
            self._log_file_path = file_path
            return self._log_file_handle
        except Exception as e:
            print(f"Error opening log file: {e}")
            return None
    
    def close_log_file(self):
        """Safely close the log file"""
        if self._log_file_handle is not None:
            try:
                if not self._log_file_handle.closed:
                    self._log_file_handle.write("\n=== Session ended ===\n")
                    self._log_file_handle.close()
                    print(f"Log file closed: {self._log_file_path}")
            except Exception as e:
                print(f"Error closing log file: {e}")
            finally:
                self._log_file_handle = None
                self._log_file_path = None
    
    def write_to_log_file(self, message: str):
        """
        Write a message to the log file
        
        Args:
            message: Message to write
        """
        if self._log_file_handle is not None and not self._log_file_handle.closed:
            try:
                self._log_file_handle.write(message + "\n")
                self._log_file_handle.flush()
            except Exception as e:
                print(f"Error writing to log file: {e}")
    
    # ==================== Font Management ====================
    
    def register_font(self, name: str, size: int, font_id: int):
        """
        Register a loaded font
        
        Args:
            name: Font name
            size: Font size
            font_id: DearPyGUI font ID
        """
        self._font_registry[(name, size)] = font_id
    
    def get_font_id(self, name: str, size: int) -> Optional[int]:
        """
        Get font ID for a given name and size
        
        Args:
            name: Font name
            size: Font size
        
        Returns:
            Font ID or None if not found
        """
        return self._font_registry.get((name, size))
    
    @property
    def current_font_name(self) -> str:
        """Get current font name"""
        return self._current_font_name
    
    @current_font_name.setter
    def current_font_name(self, name: str):
        """Set current font name"""
        self._current_font_name = name
    
    @property
    def current_font_size(self) -> int:
        """Get current font size"""
        return self._current_font_size
    
    @current_font_size.setter
    def current_font_size(self, size: int):
        """Set current font size"""
        self._current_font_size = size
    
    # ==================== UI State ====================
    
    @property
    def solver_running(self) -> bool:
        """Check if solver is currently running"""
        return self._solver_running
    
    @solver_running.setter
    def solver_running(self, value: bool):
        """Set solver running state"""
        self._solver_running = value
    
    @property
    def multigrid_enabled(self) -> bool:
        """Check if multigrid is enabled"""
        return self._multigrid_enabled
    
    @multigrid_enabled.setter
    def multigrid_enabled(self, value: bool):
        """Set multigrid enabled state"""
        self._multigrid_enabled = value
    
    @property
    def num_mesh_levels(self) -> int:
        """Get number of mesh levels"""
        return self._num_mesh_levels
    
    @num_mesh_levels.setter
    def num_mesh_levels(self, value: int):
        """Set number of mesh levels"""
        self._num_mesh_levels = max(1, min(10, value))  # Clamp between 1 and 10
    
    @property
    def current_solver_method(self) -> str:
        """Get current solver method"""
        return self._current_solver_method
    
    @current_solver_method.setter
    def current_solver_method(self, method: str):
        """Set current solver method"""
        self._current_solver_method = method
    
    # ==================== Mesh Files ====================
    
    def set_mesh_file(self, level: int, path: str):
        """
        Set mesh file path for a given level
        
        Args:
            level: Mesh level (1-indexed)
            path: File path
        """
        self._mesh_files[level] = path
    
    def get_mesh_file(self, level: int) -> Optional[str]:
        """
        Get mesh file path for a given level
        
        Args:
            level: Mesh level (1-indexed)
        
        Returns:
            File path or None if not set
        """
        return self._mesh_files.get(level)
    
    def get_all_mesh_files(self) -> Dict[int, str]:
        """Get all mesh files"""
        return self._mesh_files.copy()
    
    def clear_mesh_files(self):
        """Clear all mesh file paths"""
        self._mesh_files.clear()
    
    @property
    def init_file_path(self) -> str:
        """Get initialization file path"""
        return self._init_file_path
    
    @init_file_path.setter
    def init_file_path(self, path: str):
        """Set initialization file path"""
        self._init_file_path = path
    
    # ==================== Visualization State ====================
    
    @property
    def current_vtk_file(self) -> str:
        """Get current VTK file path"""
        return self._current_vtk_file
    
    @current_vtk_file.setter
    def current_vtk_file(self, path: str):
        """Set current VTK file path"""
        self._current_vtk_file = path
    
    @property
    def current_colormap(self) -> str:
        """Get current colormap"""
        return self._current_colormap
    
    @current_colormap.setter
    def current_colormap(self, cmap: str):
        """Set current colormap"""
        self._current_colormap = cmap
    
    @property
    def current_variable(self) -> str:
        """Get current visualization variable"""
        return self._current_variable
    
    @current_variable.setter
    def current_variable(self, var: str):
        """Set current visualization variable"""
        self._current_variable = var
    
    # ==================== Cleanup ====================
    
    def cleanup(self):
        """
        Perform cleanup of all resources
        Should be called before application exit
        """
        if self._cleanup_done:
            return
        
        print("Performing application cleanup...")
        
        # Stop convergence monitor
        self._convergence_monitor_running = False
        
        # Terminate processes
        self.terminate_plotter_process()
        self.terminate_solver_process()
        
        # Close log file
        self.close_log_file()
        
        # Mark cleanup as done
        self._cleanup_done = True
        print("Cleanup completed")
    
    # ==================== Utility Methods ====================
    
    def is_initialized(self) -> bool:
        """Check if GUI has been initialized"""
        return self._gui_initialized
    
    def mark_initialized(self):
        """Mark GUI as initialized"""
        self._gui_initialized = True
    
    def get_state_summary(self) -> Dict[str, Any]:
        """
        Get a summary of current application state
        
        Returns:
            Dictionary with state information
        """
        return {
            "solver_running": self._solver_running,
            "multigrid_enabled": self._multigrid_enabled,
            "num_mesh_levels": self._num_mesh_levels,
            "solver_method": self._current_solver_method,
            "init_file": self._init_file_path,
            "mesh_files_count": len(self._mesh_files),
            "log_file_open": self._log_file_handle is not None and not self._log_file_handle.closed,
            "plotter_active": self._plotter_process is not None and self._plotter_process.poll() is None,
            "solver_active": self._solver_process is not None and self._solver_process.poll() is None,
        }
    
    def __repr__(self) -> str:
        """String representation of state"""
        return f"AppState({self.get_state_summary()})"


# Global application state instance
# Other modules should import and use this instance
app_state = AppState()