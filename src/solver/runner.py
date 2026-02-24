"""
Solver compilation and execution module for MeMPhyS GUI

Handles compiling C source files and running the solver with
proper error handling, logging, and process management.
"""

import os
import glob
import subprocess
import threading
from pathlib import Path
from typing import List, Optional, Tuple, Callable
import dearpygui.dearpygui as dpg
from src.solver.monitoring import convergence_monitor

from src.config import (
    HEADER_DIR,
    SOLVER_SOURCE,
    SOLVER_EXECUTABLE_WINDOWS,
    SOLVER_EXECUTABLE_UNIX,
    COMPILER,
    COMPILER_FLAGS,
    SUBPROCESS_BUFFER_SIZE,
)
from src.core import logger, app_state
from src.utils import (
    get_platform,
    is_windows,
    get_executable_name,
    build_compile_command,
    validate_file_path,
)


class SolverRunner:
    """
    Manages compilation and execution of the C solver
    
    Features:
    - Automatic source file collection
    - Cross-platform compilation
    - Asynchronous execution with live output
    - Process management
    - Error handling and reporting
    """
    
    def __init__(self):
        self.compile_process: Optional[subprocess.Popen] = None
        self.run_process: Optional[subprocess.Popen] = None
        self.is_running: bool = False
        self.button_tag: Optional[str] = None
    
    def collect_source_files(self, init_file: str) -> List[str]:
        """
        Collect all C source files needed for compilation
        
        Args:
            init_file: Path to initialization .c file
        
        Returns:
            List of source file paths
        
        Raises:
            FileNotFoundError: If required files are missing
        """
        source_files = []
        init_on = dpg.get_value("init_toggle")
        # Add header files
        if os.path.exists(HEADER_DIR):
            header_c_files = glob.glob(os.path.join(HEADER_DIR, "*.c"))
            source_files.extend(header_c_files)
            logger.debug(f"Found {len(header_c_files)} header files")
        else:
            logger.warning(f"Header directory not found: {HEADER_DIR}")
        
        # Add initialization file
        if init_on:
            if init_file:
                if validate_file_path(init_file, must_exist=True):
                    source_files.append(init_file)
                    logger.debug(f"Added init file: {init_file}")
                else:
                    raise FileNotFoundError(f"Initialization file not found: {init_file}")
            else:
                raise FileNotFoundError("No initialization file specified")
        
        # Add main solver file
        if validate_file_path(SOLVER_SOURCE, must_exist=True):
            source_files.append(SOLVER_SOURCE)
            logger.debug(f"Added solver source: {SOLVER_SOURCE}")
        else:
            raise FileNotFoundError(f"Solver source not found: {SOLVER_SOURCE}")
        
        return source_files
    
    def get_output_executable(self) -> str:
        """
        Get the output executable name for the current platform
        
        Returns:
            Executable name (solver.exe on Windows, solver on Unix)
        """
        if is_windows():
            return SOLVER_EXECUTABLE_WINDOWS
        else:
            return SOLVER_EXECUTABLE_UNIX
    
    def build_compile_command(self, source_files: List[str], compiler: str = "gcc") -> List[str]:
        """
        Build the compilation command
        
        Args:
            source_files: List of source files to compile
        
        Returns:
            Command as list of arguments
        """
        output_exe = self.get_output_executable()
        if compiler == "nvc":
            # Example for NVIDIA HPC SDK compiler
            cmd = [compiler] + ["-acc", "-Minfo=accel", "-O3"] + source_files + COMPILER_FLAGS + ["-o", output_exe]
        else:
            cmd = [compiler] + source_files + COMPILER_FLAGS + ["-o", output_exe]
        
        logger.debug(f"Compile command: {' '.join(cmd)}")
        return cmd
    
    def compile_solver(self, init_file: str, compiler: str = "gcc") -> Tuple[bool, str, str]:
        """
        Compile the solver
        
        Args:
            init_file: Path to initialization file
        
        Returns:
            Tuple of (success, stdout, stderr)
        """
        try:
            logger.info("Starting compilation...")
            
            # Collect source files
            init_on = dpg.get_value("init_toggle")
            
            source_files = self.collect_source_files(init_file)
            logger.info(f"Compiling {len(source_files)} source files")
            
            # Build compile command
            compile_cmd = self.build_compile_command(source_files, compiler=compiler)
            
            # Run compilation
            self.compile_process = subprocess.Popen(
                compile_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            stdout, stderr = self.compile_process.communicate()
            returncode = self.compile_process.returncode
            
            if returncode == 0:
                logger.success("Compilation successful")
                return True, stdout, stderr
            else:
                logger.error(f"Compilation failed with code {returncode}")
                if stderr:
                    logger.error(f"Compilation errors:\n{stderr}")
                return False, stdout, stderr
        
        except FileNotFoundError as e:
            error_msg = str(e)
            logger.error(error_msg)
            return False, "", error_msg
        
        except Exception as e:
            logger.log_exception(e, "Compilation error")
            return False, "", str(e)
        
        finally:
            self.compile_process = None
    
    def run_solver_async(self, on_complete: Optional[Callable] = None):
        """
        Run the solver asynchronously in a separate thread
        
        Args:
            on_complete: Optional callback function called when solver completes
        """
        def solver_thread():
            try:
                executable = self.get_output_executable()
                
                # Check if executable exists
                if not os.path.exists(executable):
                    logger.error(f"Executable not found: {executable}")
                    self._on_solver_complete(-1, on_complete)
                    return
                
                # Run executable
                run_cmd = [f"./{executable}"] if not is_windows() else [executable]
                
                logger.info("Starting solver execution...")
                app_state.solver_running = True
                
                self.run_process = subprocess.Popen(
                    run_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=SUBPROCESS_BUFFER_SIZE,
                    universal_newlines=True
                )
                convergence_monitor.attach(self.run_process)
                convergence_monitor.start()

                app_state.solver_process = self.run_process
                
                # Wait for completion
                self.run_process.wait()
                returncode = self.run_process.returncode
                
                # Handle completion
                self._on_solver_complete(returncode, on_complete)
            
            except Exception as e:
                logger.log_exception(e, "Solver execution error")
                self._on_solver_complete(-1, on_complete)
            
            finally:
                app_state.solver_running = False
                self.run_process = None
                self.is_running = False
        
        # Start thread
        thread = threading.Thread(target=solver_thread, daemon=True)
        thread.start()
        self.is_running = True
    
    def _on_solver_complete(self, returncode: int, callback: Optional[Callable] = None):
        """
        Handle solver completion
        
        Args:
            returncode: Process return code
            callback: Optional callback function
        """
        if not dpg.is_dearpygui_running():
            return
        
        if returncode == 0:
            logger.success("Solver completed successfully")
            logger.info("Solution file saved as Solution.csv")
            
            # Update button state
            if self.button_tag:
                dpg.configure_item(self.button_tag, label="Done! Run Again")
                dpg.enable_item(self.button_tag)
        else:
            logger.error(f"Solver exited with code {returncode}")
            
            # Update button state
            if self.button_tag:
                dpg.configure_item(self.button_tag, label="Failed! Try Again")
                dpg.enable_item(self.button_tag)
        
        # Call user callback if provided
        if callback:
            try:
                callback(returncode)
            except Exception as e:
                logger.log_exception(e, "Error in solver completion callback")
    
    def compile_and_run(self, init_file: str, compiler: str = "gcc", button_tag: Optional[str] = None,
                       on_complete: Optional[Callable] = None) -> bool:
        """
        Compile and run the solver (main entry point)
        
        Args:
            init_file: Path to initialization file
            button_tag: Optional DearPyGUI tag for run button (for status updates)
            on_complete: Optional callback when solver completes
        
        Returns:
            True if compilation succeeded and solver started, False otherwise
        """
        if self.is_running:
            logger.warning("Solver is already running")
            return False
        
        self.button_tag = button_tag
        
        # Disable button
        if button_tag:
            dpg.disable_item(button_tag)
            dpg.configure_item(button_tag, label="Compiling...")
        
        # Compile
        success, stdout, stderr = self.compile_solver(init_file, compiler=compiler)
        
        if not success:
            logger.error("Cannot run solver: compilation failed")
            
            # Re-enable button
            if button_tag:
                dpg.configure_item(button_tag, label="Compilation Failed! Try Again")
                dpg.enable_item(button_tag)
            
            return False
        
        # Update button
        if button_tag:
            dpg.configure_item(button_tag, label="Running Solver...")
        
        # Run solver asynchronously
        self.run_solver_async(on_complete)
        
        return True
    
    def stop_solver(self) -> bool:
        """
        Stop the currently running solver
        
        Returns:
            True if stopped successfully, False otherwise
        """
        if not self.is_running:
            logger.warning("No solver is currently running")
            return False
        
        try:
            if self.run_process and self.run_process.poll() is None:
                logger.warning("Terminating solver process...")
                self.run_process.terminate()
                
                # Wait up to 5 seconds for graceful termination
                try:
                    self.run_process.wait(timeout=5)
                    logger.info("Solver terminated gracefully")
                except subprocess.TimeoutExpired:
                    # Force kill if it doesn't terminate
                    logger.warning("Force killing solver process...")
                    self.run_process.kill()
                    self.run_process.wait()
                    logger.info("Solver force killed")
                
                return True
            else:
                logger.info("Solver process already terminated")
                return True
        
        except Exception as e:
            logger.log_exception(e, "Error stopping solver")
            return False
        
        finally:
            self.is_running = False
            app_state.solver_running = False
            self.run_process = None
    
    def is_solver_running(self) -> bool:
        """
        Check if solver is currently running
        
        Returns:
            True if running, False otherwise
        """
        return self.is_running
    
    def get_executable_path(self) -> str:
        """
        Get the full path to the solver executable
        
        Returns:
            Absolute path to executable
        """
        return os.path.abspath(self.get_output_executable())
    
    def check_dependencies(self) -> Tuple[bool, List[str]]:
        """
        Check if all dependencies are available
        
        Returns:
            Tuple of (all_ok, list_of_missing_items)
        """
        missing = []
        
        # Check compiler
        from src.utils import check_command_exists
        if not check_command_exists(COMPILER):
            missing.append(f"Compiler: {COMPILER}")
        
        # Check solver source
        if not os.path.exists(SOLVER_SOURCE):
            missing.append(f"Solver source: {SOLVER_SOURCE}")
        
        # Check header directory
        if not os.path.exists(HEADER_DIR):
            missing.append(f"Header directory: {HEADER_DIR}")
        
        all_ok = len(missing) == 0
        
        if all_ok:
            logger.success("All solver dependencies found")
        else:
            logger.error(f"Missing dependencies: {', '.join(missing)}")
        
        return all_ok, missing
    
    def cleanup(self):
        """Clean up resources"""
        if self.is_running:
            self.stop_solver()


# Global solver runner instance
solver_runner = SolverRunner()