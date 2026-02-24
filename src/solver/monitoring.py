"""
Solver monitoring module for MeMPhyS GUI

Monitors convergence data and updates plots in real-time
during solver execution.
"""
import re
import threading
import subprocess
import numpy as np
import dearpygui.dearpygui as dpg
from typing import Optional

from src.config import CONVERGENCE_CSV, CONVERGENCE_UPDATE_INTERVAL
from src.core import logger, app_state
from src.utils import read_csv_file


class ConvergenceMonitor:
    """
    Reads convergence data from solver stdout instead of CSV polling.
    """

    _PATTERN = re.compile(
        r"Time step:\s*(\d+),\s*Steady state error:\s*([0-9eE+\-\.]+)"
    )

    def __init__(self,
                 csv_file: str = None,  # kept for compatibility
                 series_tag: str = "conv_series",
                 x_axis_tag: str = "x_axis_conv",
                 y_axis_tag: str = "y_axis_conv",
                 update_interval: float = 0.1):  # unused but kept for compatibility

        # Keep signature identical for drop-in compatibility
        self.csv_file = csv_file
        self.series_tag = series_tag
        self.x_axis_tag = x_axis_tag
        self.y_axis_tag = y_axis_tag
        self.update_interval = update_interval

        self.is_monitoring = False
        self.monitor_thread: Optional[threading.Thread] = None
        self._process: Optional[subprocess.Popen] = None

        self._x: list[float] = []
        self._y: list[float] = []

    def attach(self, process: subprocess.Popen):
        """
        Attach a running solver process.
        Call immediately after launching subprocess.
        """
        self._process = process

    def start(self):
        """Start monitoring stdout in background thread."""
        if self.is_monitoring:
            return

        if self._process is None:
            raise RuntimeError("SolverOutputMonitor: No process attached.")

        self.is_monitoring = True
        app_state.convergence_monitor_running = True

        self.monitor_thread = threading.Thread(
            target=self._read_loop,
            daemon=True,
            name="SolverOutputMonitor"
        )
        self.monitor_thread.start()

    def stop(self):
        """Stop monitoring."""
        if not self.is_monitoring:
            return

        self.is_monitoring = False
        app_state.convergence_monitor_running = False

        if self.monitor_thread:
            self.monitor_thread.join(timeout=1.0)

    def cleanup(self):
        """Cleanup resources."""
        self.stop()

    def is_running(self) -> bool:
        return self.is_monitoring

    def reset(self):
        """Reset convergence plot to initial state."""
        self._x.clear()
        self._y.clear()

        if dpg.does_item_exist(self.series_tag):
            dpg.set_value(self.series_tag, [[], []])

        if dpg.does_item_exist(self.x_axis_tag):
            dpg.set_axis_limits(self.x_axis_tag, 0, 100)

        if dpg.does_item_exist(self.y_axis_tag):
            dpg.set_axis_limits(self.y_axis_tag, 1e-10, 1)

    def force_update(self):
        """Force GUI update with current buffered data."""
        self._update_plot()

    def get_latest_data(self) -> Optional[tuple]:
        if len(self._x) == 0:
            return None
        return (np.array(self._x), np.array(self._y))

    def get_convergence_status(self) -> dict:
        if len(self._y) == 0:
            return {
                "has_data": False,
                "iterations": 0,
                "current_error": None,
                "min_error": None,
            }

        y = np.array(self._y)
        x = np.array(self._x)

        return {
            "has_data": True,
            "iterations": int(x[-1]) if len(x) > 0 else 0,
            "current_error": float(y[-1]),
            "min_error": float(np.min(y[y > 0])) if np.any(y > 0) else None,
            "is_converging": self._is_converging(y),
        }


    def _read_loop(self):
        """Continuously read solver stdout."""
        for raw_line in self._process.stdout:
            if not self.is_monitoring:
                break

            line = raw_line.strip()

            # Forward to GUI log if desired
            if dpg.is_dearpygui_running():
                logger.info(line)

            match = self._PATTERN.search(line)
            if match:
                timestep = float(match.group(1))
                error = float(match.group(2))

                self._x.append(timestep)
                self._y.append(error)

                self._update_plot()

    def _update_plot(self):
        if not dpg.is_dearpygui_running():
            return

        if len(self._x) == 0:
            return

        if dpg.does_item_exist(self.series_tag):
            dpg.set_value(self.series_tag, [self._x, self._y])

        self._update_axis_limits()

    def _update_axis_limits(self):
        if not self._x or not self._y:
            return

        # X-axis
        if dpg.does_item_exist(self.x_axis_tag):
            dpg.set_axis_limits(self.x_axis_tag, 0, self._x[-1] + 10)

        # Y-axis (log scale)
        positive_y = [v for v in self._y if v > 0]
        if not positive_y:
            return

        import math
        y_min = 10 ** (math.floor(math.log10(min(positive_y))) - 1)
        y_max = 10 ** (math.ceil(math.log10(max(positive_y))) + 0.2)

        if dpg.does_item_exist(self.y_axis_tag):
            dpg.set_axis_limits(self.y_axis_tag, y_min, y_max)

    def _is_converging(self, y: np.ndarray, window: int = 10) -> bool:
        if len(y) < window:
            return False

        recent = y[-window:]
        return recent[-1] < recent[0]

# Global instance
convergence_monitor = ConvergenceMonitor()