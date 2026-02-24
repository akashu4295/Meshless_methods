"""
Solver package for MeMPhyS GUI

This package handles solver compilation, execution, and monitoring:
- SolverRunner: Compiles and runs the C solver
- ConvergenceMonitor: Monitors convergence in real-time

"""

from .runner import SolverRunner, solver_runner
from .monitoring import convergence_monitor, ConvergenceMonitor

__all__ = [
    # Classes
    'SolverRunner',
    'ConvergenceMonitor',
    
    # Global instances
    'solver_runner',
    'convergence_monitor',
]