"""
Logging utilities for MeMPhyS GUI

Handles both file logging and GUI log window updates with proper
thread safety and error handling.
"""

import os
from datetime import datetime
from typing import Optional
import dearpygui.dearpygui as dpg
from src.config import LOG_SCROLL_THRESHOLD
import queue

class Logger:
    """
    Centralized logging system for both file and GUI output
    
    Features:
    - Thread-safe logging
    - Automatic file and GUI logging
    - Log levels (INFO, WARNING, ERROR, SUCCESS)
    - Auto-scrolling log window
    - Timestamping
    """
    
    # Log level constants
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    SUCCESS = "SUCCESS"
    DEBUG = "DEBUG"
    
    # ANSI color codes for console output
    COLORS = {
        INFO: "\033[0m",      # Default
        WARNING: "\033[93m",  # Yellow
        ERROR: "\033[91m",    # Red
        SUCCESS: "\033[92m",  # Green
        DEBUG: "\033[94m",    # Blue
    }
    RESET = "\033[0m"
    
    def __init__(self, log_window_tag: str = "log_window", 
                 log_child_tag: str = "log_child",
                 file_handle=None):
        """
        Initialize the logger
        
        Args:
            log_window_tag: DearPyGUI tag for the log text widget
            log_child_tag: DearPyGUI tag for the scrollable child window
            file_handle: File handle for log file (optional)
        """
        self.log_window_tag = log_window_tag
        self.log_child_tag = log_child_tag
        self.file_handle = file_handle
        self._enable_timestamps = True
        self._enable_console = True
        self._enable_file = True
        self._enable_gui = True
        self._gui_queue = queue.Queue()
    
    def set_file_handle(self, file_handle):
        """Update the file handle for logging"""
        self.file_handle = file_handle
    
    def set_enable_timestamps(self, enable: bool):
        """Enable or disable timestamps in logs"""
        self._enable_timestamps = enable
    
    def set_enable_console(self, enable: bool):
        """Enable or disable console output"""
        self._enable_console = enable
    
    def set_enable_file(self, enable: bool):
        """Enable or disable file logging"""
        self._enable_file = enable
    
    def set_enable_gui(self, enable: bool):
        """Enable or disable GUI logging"""
        self._enable_gui = enable
    
    def _format_message(self, message: str, level: str = INFO, 
                       include_timestamp: bool = True) -> str:
        """
        Format a log message with optional timestamp and level
        
        Args:
            message: The message to log
            level: Log level (INFO, WARNING, ERROR, SUCCESS, DEBUG)
            include_timestamp: Whether to include timestamp
        
        Returns:
            Formatted message string
        """
        parts = []
        
        if include_timestamp and self._enable_timestamps:
            timestamp = datetime.now().strftime("%H:%M:%S")
            parts.append(f"[{timestamp}]")
        
        if level != self.INFO:
            parts.append(f"[{level}]")
        
        parts.append(message)
        
        return " ".join(parts)
    
    def _write_to_file(self, message: str):
        """Write message to log file"""
        if not self._enable_file:
            return
            
        if self.file_handle is not None:
            try:
                if not self.file_handle.closed:
                    self.file_handle.write(message + "\n")
                    self.file_handle.flush()
            except Exception as e:
                print(f"Error writing to log file: {e}")
    
    def _write_to_gui(self, message: str):
        """Queue message for GUI - safe to call from any thread"""
        if not self._enable_gui:
            return
        self._gui_queue.put(message)

    def flush_gui_queue(self):
        try:
            if not dpg.is_dearpygui_running():
                return
            if not dpg.does_item_exist(self.log_child_tag):
                return

            messages = []
            while not self._gui_queue.empty():
                try:
                    messages.append(self._gui_queue.get_nowait())
                except queue.Empty:
                    break

            if not messages:
                return

            for msg in messages:
                dpg.add_text(msg, parent=self.log_child_tag)

            dpg.set_y_scroll(self.log_child_tag, -1.0)

        except Exception as e:
            print(f"flush error: {e}")
    
    def _write_to_console(self, message: str, level: str = INFO):
        """Write message to console with color"""
        if not self._enable_console:
            return
            
        color = self.COLORS.get(level, self.COLORS[self.INFO])
        print(f"{color}{message}{self.RESET}")
    
    def _auto_scroll(self):
        """Auto-scroll log window if user is near bottom"""
        try:
            # Check if DearPyGUI is running
            if not dpg.is_dearpygui_running():
                return
                
            if not dpg.does_item_exist(self.log_child_tag):
                return
            
            max_scroll = dpg.get_y_scroll_max(self.log_child_tag)
            current_scroll = dpg.get_y_scroll(self.log_child_tag)
            
            # If within threshold of bottom, scroll to bottom
            if abs(max_scroll - current_scroll) < LOG_SCROLL_THRESHOLD:
                dpg.set_y_scroll(self.log_child_tag, -1)
        except Exception:
            # Silently fail if scrolling doesn't work
            pass
    
    def log(self, message: str, level: str = INFO, 
            to_file: bool = True, to_gui: bool = True, 
            to_console: bool = True):
        """
        Log a message to all enabled outputs
        
        Args:
            message: The message to log
            level: Log level (INFO, WARNING, ERROR, SUCCESS, DEBUG)
            to_file: Write to file
            to_gui: Write to GUI
            to_console: Write to console
        """
        formatted = self._format_message(message, level)
        
        if to_file and self._enable_file:
            self._write_to_file(formatted)
        
        if to_gui and self._enable_gui:
            self._write_to_gui(formatted)
        
        if to_console and self._enable_console:
            self._write_to_console(formatted, level)
    
    def info(self, message: str):
        """Log an info message"""
        self.log(message, self.INFO)
    
    def warning(self, message: str):
        """Log a warning message"""
        self.log(message, self.WARNING)
    
    def error(self, message: str):
        """Log an error message"""
        self.log(message, self.ERROR)
    
    def success(self, message: str):
        """Log a success message"""
        self.log(message, self.SUCCESS)
    
    def debug(self, message: str):
        """Log a debug message"""
        self.log(message, self.DEBUG)
    
    def separator(self, char: str = "-", length: int = 60):
        """Log a separator line"""
        self.log(char * length, self.INFO)
    
    def header(self, message: str, char: str = "=", length: int = 60):
        """
        Log a header with separator lines
        
        Args:
            message: Header text
            char: Character for separator
            length: Length of separator
        """
        self.separator(char, length)
        self.log(message, self.INFO)
        self.separator(char, length)

    def clear(self):
        """Clear the GUI log window"""
        try:
            if dpg.does_item_exist(self.log_child_tag):
                dpg.delete_item(self.log_child_tag, children_only=True)
                dpg.set_y_scroll(self.log_child_tag, 0.0)
        except Exception as e:
            print(f"Error clearing log window: {e}")
    
    def log_exception(self, exception: Exception, context: str = ""):
        """
        Log an exception with context
        
        Args:
            exception: The exception to log
            context: Additional context about where the exception occurred
        """
        if context:
            self.error(f"{context}: {type(exception).__name__}: {str(exception)}")
        else:
            self.error(f"{type(exception).__name__}: {str(exception)}")
    
    def log_dict(self, data: dict, title: str = "Data"):
        """
        Log a dictionary in a formatted way
        
        Args:
            data: Dictionary to log
            title: Title for the data
        """
        self.info(f"{title}:")
        for key, value in data.items():
            self.info(f"  {key}: {value}")
    
    def log_list(self, items: list, title: str = "Items"):
        """
        Log a list in a formatted way
        
        Args:
            items: List to log
            title: Title for the list
        """
        self.info(f"{title}:")
        for i, item in enumerate(items, 1):
            self.info(f"  {i}. {item}")
    
    def log_file_operation(self, operation: str, filepath: str, 
                          success: bool = True):
        """
        Log a file operation
        
        Args:
            operation: Type of operation (read, write, delete, etc.)
            filepath: Path to the file
            success: Whether operation was successful
        """
        filename = os.path.basename(filepath)
        if success:
            self.success(f"{operation} {filename}")
        else:
            self.error(f"Failed to {operation.lower()} {filename}")
    
    def log_process(self, process_name: str, status: str):
        """
        Log a process status
        
        Args:
            process_name: Name of the process
            status: Status (starting, running, completed, failed, etc.)
        """
        status_lower = status.lower()
        if "fail" in status_lower or "error" in status_lower:
            self.error(f"{process_name}: {status}")
        elif "success" in status_lower or "complet" in status_lower:
            self.success(f"{process_name}: {status}")
        elif "warn" in status_lower:
            self.warning(f"{process_name}: {status}")
        else:
            self.info(f"{process_name}: {status}")
    
    def log_parameter(self, param_name: str, value, unit: str = ""):
        """
        Log a parameter value
        
        Args:
            param_name: Parameter name
            value: Parameter value
            unit: Optional unit
        """
        unit_str = f" {unit}" if unit else ""
        self.info(f"{param_name}: {value}{unit_str}")


def clear_logs():
    """Clear the log window (backward compatible)"""
    logger.clear()

# Create a global logger instance
# Other modules should import this instance
logger = Logger(log_window_tag="log_window", log_child_tag="log_child")