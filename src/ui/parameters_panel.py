"""
Parameters panel construction for MeMPhyS GUI

Creates the left panel with solver parameters, mesh file selection,
and run controls. Sections are collapsible like a VS Code file explorer.
"""

import dearpygui.dearpygui as dpg

from src.config import (
    BASE_PARAMETERS,
    IMPLICIT_PARAMETERS,
    MULTIGRID_PARAMETERS,
    PARAMETER_LABELS,
    PARAMETER_TOOLTIPS,
    SOLVER_METHODS,
    Poisson_SOLVER_METHODS,
    DEFAULT_SOLVER_METHOD,
    MAX_MESH_LEVELS,
    LEFT_PANEL_WIDTH,
    PARAMETER_INPUT_WIDTH,
    MESH_PATH_INPUT_WIDTH,
    RUN_BUTTON_WIDTH,
    RUN_BUTTON_HEIGHT,
    COLORS,
    FILE_EXTENSIONS,
    FILE_DIALOG_WIDTH,
    FILE_DIALOG_HEIGHT,
    themed_texts,
)
from src.callbacks import (
    run_solver_callback,
    validate_numeric_input,
    show_implicit_callback,
    show_multigrid_callback,
    update_mesh_inputs_callback,
    open_file_dialog_callback,
    select_mesh_file_callback,
    launch_gmsh_callback,
    set_mesh_from_geometry_callback,
    browse_geometry_file_callback,
    select_geometry_file_callback,
    show_bc_window_callback,
    show_restart_callback,
    show_init_callback,
    stop_solver_callback,
)


def create_parameters_panel(themes: dict) -> int:
    with dpg.child_window(width=LEFT_PANEL_WIDTH, border=True,
        tag="parameters_panel") as panel:
        with _explorer_section("GEOMETRY AND GRIDS", default_open=True):
            _create_geometry_section(themes)
            _section_gap()
            _create_multigrid_section()
            _section_gap()
            _create_mesh_files_section(themes)

        with _explorer_section("SOLVER CONFIGURATION", default_open=True):
            _create_solver_method_section(themes)
            _section_gap()
            _create_poisson_solver_method_section()

        with _explorer_section("FLOW PARAMETERS", default_open=True):
            _create_flow_parameters_section()
            _create_implicit_parameters_section()

        with _explorer_section("INITIAL AND RESTART CONDITIONS", default_open=False):
            _create_init_file_section()
            dpg.add_spacer(height=6)
            _create_restart_file_section()

        _section_gap()
        _create_run_button(themes)

    return panel


# ── Visual helpers ─────────────────────────────────────────────────────────────

class _explorer_section:
    """
    Context manager that wraps a dpg.collapsing_header styled to a
    sidebar section header (all-caps label, subtle separator,
    indented content).
    """

    def __init__(self, title: str, default_open: bool = True):
        self._title = title
        self._default_open = default_open
        self._header = None

    def __enter__(self):
        self._header = dpg.add_collapsing_header(
            label=self._title,
            default_open=self._default_open,
            bullet=False,
        )
        dpg.push_container_stack(self._header)
        dpg.add_spacer(height=3)
        self._indent = dpg.add_group(indent=8)
        dpg.push_container_stack(self._indent)
        return self._header

    def __exit__(self, *_):
        dpg.add_spacer(height=4)
        dpg.pop_container_stack()  # indent group
        dpg.pop_container_stack()  # collapsing header


def _section_gap():
    """Consistent vertical breathing room between sub-sections."""
    dpg.add_spacer(height=6)

# ── Sections ───────────────────────────────────────────────────────────────────

def _create_solver_method_section(themes: dict):
    with dpg.group(horizontal=True):
        themed_texts("Method","subheader")
        dpg.add_combo(
            items=SOLVER_METHODS,
            default_value=DEFAULT_SOLVER_METHOD,
            tag="solver_method",
            width=PARAMETER_INPUT_WIDTH + 80,
            callback=show_implicit_callback
        )
    with dpg.tooltip("solver_method"):
        themed_texts("Select the solver method for the Navier-Stokes equations","info")
        dpg.add_separator()
        themed_texts("Implicit: more stable for large time steps, solves a linear system each step.","info")
        themed_texts("Explicit: faster for small meshes/time steps, needs careful stability tuning.","info")


def _create_poisson_solver_method_section():
    with dpg.group(horizontal=True, tag="poisson_solver_method_group"):
        themed_texts("Poisson Solver", "subheader")
        dpg.add_combo(
            items=Poisson_SOLVER_METHODS,
            default_value=Poisson_SOLVER_METHODS[0],
            tag="poisson_solver_method",
            width=PARAMETER_INPUT_WIDTH + 25,
        )
    with dpg.tooltip("poisson_solver_method_group"):
        themed_texts("Select the solver type for the Poisson equation",
                     "info")
        
        dpg.add_separator()
        themed_texts("Multigrid: efficient for large problems.","info")
        themed_texts("Direct: faster for smaller meshes.","info")

def _create_flow_parameters_section():
    for pname, pval in BASE_PARAMETERS.items():
        _add_parameter_input(pname, pval)


def _create_implicit_parameters_section():
    for pname, pval in IMPLICIT_PARAMETERS.items():
        _add_parameter_input(pname, pval, show=False)


def _create_multigrid_section():
    dpg.add_checkbox(
        label="Enable Multigrid",
        tag="multigrid_toggle",
        callback=show_multigrid_callback
    )
    with dpg.tooltip("multigrid_toggle"):
        themed_texts("Enable multigrid solver with multiple mesh levels", "info")
        dpg.add_separator()
        themed_texts("Provide mesh files finest -> coarsest.", "info")
        themed_texts("Number of files must match the mesh levels specified below.", "info")

    with dpg.group(horizontal=False, tag="multigrid_parameters_section", show=False):
        dpg.add_spacer(height=4)
        themed_texts("Multigrid Parameters", "subheader")
        for pname, pval in MULTIGRID_PARAMETERS.items():
            _add_parameter_input(pname, pval)
        dpg.add_input_int(
            label="Mesh Levels",
            default_value=1,
            min_value=1,
            max_value=MAX_MESH_LEVELS,
            width=PARAMETER_INPUT_WIDTH,
            tag="num_mesh_levels",
            callback=update_mesh_inputs_callback
        )
        dpg.add_spacer(height=4)

def _create_geometry_section(themes: dict):
    themed_texts("Create geometry in Gmsh, use a .geo file, or browse for an existing .msh file",
        "success", wrap=LEFT_PANEL_WIDTH - 36)
    dpg.add_spacer(height=4)

    with dpg.group(horizontal=True):
        new_geo_btn = dpg.add_button(
            label="Open Gmsh GUI",
            callback=launch_gmsh_callback,
            width=130
        )
        open_geo_btn = dpg.add_button(
            label="Use .geo File",
            callback=browse_geometry_file_callback,
            width=130
        )

    with dpg.tooltip(new_geo_btn):
        themed_texts("Open the Gmsh GUI for geometry creation", "info")
        dpg.add_separator()
        themed_texts("1. Name boundaries for correct BC association", "info")
        themed_texts("2. Save mesh as .msh (ASCII v2)", "info")
        themed_texts("3. Load the .msh file below", "info")

    with dpg.tooltip(open_geo_btn):
        themed_texts("Select a .geo file to generate the mesh from", "info")
        dpg.add_separator()
        themed_texts("1. Name boundaries for correct BC association", "info")
        themed_texts("2. Save and select the .geo file here", "info")
        themed_texts("3. Mesh is generated and loaded automatically", "info")

    if "button_secondary" in themes:
        dpg.bind_item_theme(new_geo_btn, themes["button_secondary"])
        dpg.bind_item_theme(open_geo_btn, themes["button_secondary"])

    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_geometry",
        callback=select_geometry_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT
    )
    dpg.add_file_extension(".geo",          parent="file_dialog_geometry", color=(150, 255, 150, 255))
    dpg.add_file_extension(".geo_unrolled", parent="file_dialog_geometry", color=(150, 255, 150, 255))
    dpg.add_spacer(height=4)


def _create_mesh_files_section(themes: dict):
    with dpg.group(tag="multigrid_mesh_section"):
        themed_texts("Mesh files  (finest to coarsest)",
                     "success", tag="multigrid_hint_text", show=False)
        dpg.add_spacer(height=3)

        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            show_initial = (i == 0)
            with dpg.group(horizontal=True, show=show_initial, tag=f"mesh_group_{idx}"):
                dpg.add_input_text(
                    tag=f"mesh_file_{idx}",
                    hint=f"Mesh {idx}",
                    width=MESH_PATH_INPUT_WIDTH,
                    show=show_initial
                )
                dpg.add_button(
                    label="Browse",
                    tag=f"browse_{idx}",
                    callback=open_file_dialog_callback,
                    user_data=idx,
                    show=show_initial
                )

        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            dialog_tag = f"file_dialog_{idx}"
            dpg.add_file_dialog(
                directory_selector=False,
                tag=dialog_tag,
                user_data=f"mesh_file_{idx}",
                callback=select_mesh_file_callback,
                show=False,
                width=FILE_DIALOG_WIDTH,
                height=FILE_DIALOG_HEIGHT
            )
            if ".msh" in FILE_EXTENSIONS:
                dpg.add_file_extension(".msh", parent=dialog_tag, color=FILE_EXTENSIONS[".msh"])

    dpg.add_spacer(height=6)
    associate_bc_btn = dpg.add_button(
        label="Set Boundary Conditions",
        callback=show_bc_window_callback,
        width=270
    )
    with dpg.tooltip(associate_bc_btn):
        themed_texts("Open the boundary condition association window",
                     "info")
    if "button_secondary" in themes:
        dpg.bind_item_theme(associate_bc_btn, themes["button_secondary"])


def _create_init_file_section():
    dpg.add_checkbox(
        label="Initialise with user-defined conditions",
        tag="init_toggle",
        callback=show_init_callback
    )
    with dpg.tooltip("init_toggle"):
        themed_texts("Provide a .c file defining custom initial conditions",
                     "info")
        dpg.add_separator()
        themed_texts("See 'init.c' in the repository for a template.", "info")

    with dpg.group(horizontal=True, show=False, tag="init_group"):
        dpg.add_input_text(hint="Path to initialisation file", tag="init_path",
                           width=MESH_PATH_INPUT_WIDTH)
        dpg.add_button(label="Browse", tag="init_browse",
                       callback=open_file_dialog_callback, user_data="init")

    dpg.add_file_dialog(
        directory_selector=False, tag="file_dialog_init",
        user_data="init_path", callback=select_mesh_file_callback,
        show=False, width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT
    )
    if ".c" in FILE_EXTENSIONS:
        dpg.add_file_extension(".c", parent="file_dialog_init",
                               color=FILE_EXTENSIONS[".c"])


def _create_restart_file_section():
    dpg.add_checkbox(
        label="Restart from previous run",
        tag="restart_toggle",
        callback=show_restart_callback
    )
    with dpg.tooltip("restart_toggle"):
        themed_texts("Provide a .csv restart file from a previous run",
                     "info")
        dpg.add_separator()
        themed_texts("Must use the same mesh file as the original run.", "info")

    with dpg.group(horizontal=True, show=False, tag="restart_group"):
        dpg.add_input_text(hint="Path to restart file", tag="restart_path",
                           width=MESH_PATH_INPUT_WIDTH)
        dpg.add_button(label="Browse", tag="restart_browse",
                       callback=open_file_dialog_callback, user_data="restart")

    dpg.add_file_dialog(
        directory_selector=False, tag="file_dialog_restart",
        user_data="restart_path", callback=select_mesh_file_callback,
        show=False, width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT
    )
    if ".csv" in FILE_EXTENSIONS:
        dpg.add_file_extension(".csv", parent="file_dialog_restart",
                               color=FILE_EXTENSIONS[".csv"])

def _create_run_button(themes: dict):
    dpg.add_checkbox(label="Run on GPU (if supported)", tag="gpu_toggle")
    dpg.add_spacer(height=6)
    run_btn = dpg.add_button(
        label="Compile and Run Solver",
        width=270,
        height=RUN_BUTTON_HEIGHT,
        callback=run_solver_callback,
        tag="run_solver_button"
    )
    with dpg.tooltip(run_btn):
        themed_texts("Compile and execute the solver", "info")
        dpg.add_separator()
        themed_texts("1. Writes configuration files", "info")
        themed_texts("2. Compiles C sources with gcc", "info")
        themed_texts("3. Runs the solver executable", "info")

    dpg.add_spacer(height=4)
    stop_btn = dpg.add_button(
        label="Stop Solver",
        width=270,
        height=RUN_BUTTON_HEIGHT,
        callback=stop_solver_callback,
        tag="stop_solver_button",
        enabled=False  # disabled until solver starts
    )
    with dpg.tooltip(stop_btn):
        themed_texts("Terminate the running solver process", "info")

    if "button" in themes:
        dpg.bind_item_theme(run_btn, themes["button"])
    if "button_danger" in themes:  # optional red theme
        dpg.bind_item_theme(stop_btn, themes["button_danger"])

def _add_parameter_input(pname: str, pval, show: bool = True):
    """Helper to add a labelled input with tooltip for a single parameter."""
    tag = f"param_{pname}"
    label = PARAMETER_LABELS.get(pname, pname)  # fall back to raw name if missing
    
    dpg.add_input_text(
        label=label,
        tag=tag,
        default_value=str(pval),
        width=PARAMETER_INPUT_WIDTH,
        callback=validate_numeric_input,
        on_enter=True,
        show=show,
    )

    if pname in PARAMETER_TOOLTIPS:
        description, suggestion = PARAMETER_TOOLTIPS[pname]
        with dpg.tooltip(tag):
            themed_texts(description, "info")
            dpg.add_separator()
            themed_texts(suggestion, "info")