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
    COLORS,
    BC_TYPES,
    BC_VARIABLES,
    BC_WINDOW_WIDTH,
    BC_WINDOW_HEIGHT,
)

from src.utils.gmsh_bc_manager import (
    get_physical_names,
    get_boundary_condition,
    set_boundary_condition,
)


# ---------------------------------------------------------------------------
# Velocity input mode constants
# ---------------------------------------------------------------------------

#: Available velocity specification modes shown in the per-row dropdown
VELOCITY_INPUT_MODES = [
    "Components (u, v, w)",       # Cartesian x/y/z components
    "Normal + Tangential (n, t)", # 1 normal + 1 tangential (2D-style)
]

#: Variable names produced by each velocity input mode
VELOCITY_MODE_VARS = {
    "Components (u, v, w)":       ["u", "v", "w"],
    "Normal + Tangential (n, t)": ["v_n", "v_t"],
}

#: Human-readable labels shown next to each input field, keyed by variable name
VELOCITY_VAR_LABELS = {
    "u":   "u",
    "v":   "v",
    "w":   "w",
    "v_n": "n",
    "v_t": "t",
}

# Column widths for the BC table layout
_COL_NAME_W    = 100  # physical entity name column
_COL_BCTYPE_W  = 130  # BC type combo column
_COL_VELMODE_W = 165  # velocity mode combo column
_COL_VAR_INPUT_W = 72 # float input width


# ---------------------------------------------------------------------------
# Window creation
# ---------------------------------------------------------------------------

def show_bc_window_callback(sender, app_data, user_data):
    """
    Show the Boundary Conditions dialog.

    If the window already exists it is brought to focus; otherwise it is
    created fresh and the BC rows are populated immediately from the mesh.

    Args:
        sender:    Button tag that triggered this callback.
        app_data:  Unused.
        user_data: Unused.
    """
    if dpg.does_item_exist("bc_window"):
        dpg.configure_item("bc_window", show=True)
        dpg.focus_item("bc_window")
        return

    with dpg.window(
        label="Boundary Conditions",
        tag="bc_window",
        modal=True,
        width=BC_WINDOW_WIDTH,
        height=BC_WINDOW_HEIGHT,
        no_resize=True,
        pos=(415, 275),
    ):
        dpg.add_spacer(height=4)
        dpg.add_text(
            "Reading boundary conditions from mesh...",
            tag="bc_dialog_instructions",
            color=COLORS["info"],
        )
        dpg.add_spacer(height=3)
        dpg.add_separator()
        dpg.add_spacer(height=3)

        with dpg.group(horizontal=True, tag="bc_column_headers"):
            dpg.add_text("Entity",        color=COLORS["subheader"], indent=4)
            dpg.add_spacer(width=_COL_NAME_W - 46)
            dpg.add_text("BC Type",       color=COLORS["subheader"])
            dpg.add_spacer(width=_COL_BCTYPE_W - 44)
            dpg.add_text("Velocity Mode", color=COLORS["subheader"])
            dpg.add_spacer(width=_COL_VELMODE_W - 88)
            dpg.add_text("Values",        color=COLORS["subheader"])

        dpg.add_separator()
        dpg.add_spacer(height=2)

        with dpg.child_window(
            tag="bc_widgets_dialog_container",
            height=BC_WINDOW_HEIGHT - 185,
            border=False,
            horizontal_scrollbar=False,
        ):
            pass  # rows injected by refresh_bc_callback()

        dpg.add_spacer(height=3)
        dpg.add_separator()
        dpg.add_spacer(height=3)

        with dpg.group(horizontal=True):
            dpg.add_button(
                label="Load bc.csv",
                tag="load_bc_dialog_button",
                callback=lambda: load_bc_csv_dialog_callback(),
            )
            dpg.add_spacer(width=BC_WINDOW_WIDTH - 320)
            dpg.add_button(
                label="Apply",
                width=80,
                callback=apply_bc_callback,
            )
            dpg.add_spacer(width=4)
            dpg.add_button(
                label="Close",
                width=80,
                callback=lambda: dpg.configure_item("bc_window", show=False),
            )

    refresh_bc_callback()

def refresh_bc_callback():
    """
    Read physical names from the currently loaded mesh and rebuild the BC rows.

    Clears ``bc_widgets_dialog_container`` and repopulates it with one
    :func:`create_bc_widget` row per physical name.
    """
    from src.utils.gmsh_bc_manager import read_physical_names_from_msh

    mesh_file = dpg.get_value("mesh_file_1") if dpg.does_item_exist("mesh_file_1") else None

    if not mesh_file:
        if dpg.does_item_exist("bc_dialog_instructions"):
            dpg.set_value("bc_dialog_instructions", "No mesh file selected.")
        logger.warning("No mesh file selected. Please select a mesh file first.")
        return

    physical_names = read_physical_names_from_msh(mesh_file)

    if not physical_names:
        if dpg.does_item_exist("bc_dialog_instructions"):
            dpg.set_value("bc_dialog_instructions", "No physical names found in mesh file.")
        logger.warning("No physical names found in mesh file.")
        return

    if dpg.does_item_exist("bc_widgets_dialog_container"):
        dpg.delete_item("bc_widgets_dialog_container", children_only=True)

    if dpg.does_item_exist("bc_dialog_instructions"):
        dpg.set_value(
            "bc_dialog_instructions",
            f"Configuring boundary conditions for {len(physical_names)} physical entities:",
        )

    for i, physical_name in enumerate(physical_names):
        create_bc_widget(physical_name)
        # Small gap between rows, but not after the last one
        if i < len(physical_names) - 1:
            dpg.add_spacer(height=3, parent="bc_widgets_dialog_container")


def create_bc_widget(physical_name: str):
    """
    Add one aligned row of BC controls for *physical_name* to the container.

    Layout (horizontal group, explicit spacers for column alignment):
        [entity label] [BC type combo] [velocity mode combo*] [var inputs*]

    The velocity-mode combo is only shown when the selected BC type exposes
    velocity variables.  Variable inputs update in-memory state on every
    value change — no Enter required.

    Args:
        physical_name: Name of the physical entity (e.g. ``"inlet"``).
    """
    bc_data          = get_boundary_condition(physical_name)
    default_type     = bc_data["type"]      if bc_data else BC_TYPES[0]
    default_vars     = bc_data["variables"] if bc_data else {}
    has_vel          = _bc_type_has_velocity(default_type)
    default_vel_mode = _infer_velocity_mode(default_vars)

    with dpg.group(horizontal=True, parent="bc_widgets_dialog_container"):
        dpg.add_text(physical_name, color=COLORS["subheader"], indent=4)
        pad = max(4, _COL_NAME_W - len(physical_name) * 7 - 4)
        dpg.add_spacer(width=pad)

        dpg.add_combo(
            items=BC_TYPES,
            default_value=default_type,
            tag=f"bc_type_{physical_name}",
            width=_COL_BCTYPE_W,
            user_data=physical_name,                                          # <-- pass here
            callback=lambda s, a, u: _on_bc_type_changed(u, a),              # <-- u is user_data
        )
        dpg.add_spacer(width=6)

        dpg.add_combo(
            items=VELOCITY_INPUT_MODES,
            default_value=default_vel_mode,
            tag=f"bc_velmode_{physical_name}",
            width=_COL_VELMODE_W,
            show=has_vel,
            user_data=physical_name,                                          # <-- pass here
            callback=lambda s, a, u: _on_velocity_mode_changed(u, a),        # <-- u is user_data
        )
        dpg.add_spacer(
            width=6,
            tag=f"bc_velmode_spacer_{physical_name}",
            show=has_vel,
        )

        with dpg.group(horizontal=True, tag=f"bc_vars_{physical_name}"):
            _create_variable_inputs(physical_name, default_type, default_vel_mode, default_vars)


def _create_variable_inputs(
    physical_name: str,
    bc_type: str,
    velocity_mode: str,
    current_values: dict,
):
    """
    Populate the ``bc_vars_<n>`` group with labelled float inputs.

    For BC types with velocity variables the names come from
    :data:`VELOCITY_MODE_VARS`[*velocity_mode*]; for all others from
    :data:`BC_VARIABLES`[*bc_type*].

    Inputs call :func:`_on_var_changed` on every value change — no Enter needed.

    Args:
        physical_name:  Physical entity name.
        bc_type:        Selected BC type string.
        velocity_mode:  Selected velocity input mode string.
        current_values: Dict of ``{var_name: value}`` to pre-fill.
    """
    var_names = _resolve_variable_names(bc_type, velocity_mode)

    for var_name in var_names:
        label       = VELOCITY_VAR_LABELS.get(var_name, var_name)
        default_val = float(current_values.get(var_name, 0.0))

        dpg.add_text(f"{label}=")
        dpg.add_input_float(
            tag=f"bc_var_{physical_name}_{var_name}",
            default_value=default_val,
            width=_COL_VAR_INPUT_W,
            format="%.3f",
            step=0,
            callback=lambda s, a, u=(physical_name, var_name): _on_var_changed(*u, a),
        )
        dpg.add_spacer(width=4)

def _on_bc_type_changed(physical_name: str, new_type: str):
    """
    Rebuild the velocity-mode combo visibility and variable inputs when the
    BC type combo changes.

    Args:
        physical_name: Physical entity name.
        new_type:      Newly selected BC type string.
    """
    has_vel     = _bc_type_has_velocity(new_type)
    velmode_tag = f"bc_velmode_{physical_name}"
    spacer_tag  = f"bc_velmode_spacer_{physical_name}"

    if dpg.does_item_exist(velmode_tag):
        dpg.configure_item(velmode_tag, show=has_vel)
    if dpg.does_item_exist(spacer_tag):
        dpg.configure_item(spacer_tag, show=has_vel)

    vel_mode = (
        dpg.get_value(velmode_tag)
        if has_vel and dpg.does_item_exist(velmode_tag)
        else VELOCITY_INPUT_MODES[0]
    )
    _rebuild_variable_inputs(physical_name, new_type, vel_mode, current_values={})
    update_bc_from_widgets_dialog(physical_name)

def _on_velocity_mode_changed(physical_name: str, new_mode: str):
    """
    Rebuild variable inputs when the velocity decomposition mode changes.

    Variable values reset to 0 because names change between modes
    (e.g. ``u/v/w`` → ``v_n/v_t``).

    Args:
        physical_name: Physical entity name.
        new_mode:      Newly selected velocity input mode string.
    """
    bc_type = dpg.get_value(f"bc_type_{physical_name}")
    _rebuild_variable_inputs(physical_name, bc_type, new_mode, current_values={})
    update_bc_from_widgets_dialog(physical_name)


def _on_var_changed(physical_name: str, var_name: str, new_value: float):
    """
    Called on every value change in a variable input field (no Enter needed).

    Args:
        physical_name: Physical entity name.
        var_name:      Variable name (e.g. ``"u"``, ``"v_n"``).
        new_value:     Current float value from the widget.
    """
    update_bc_from_widgets_dialog(physical_name)

def apply_bc_callback(sender, app_data, user_data):
    """
    Flush all BC widget values to memory, validate, and write ``bc.csv``.

    Args:
        sender:    Button tag.
        app_data:  Unused.
        user_data: Unused.
    """
    from src.utils.gmsh_bc_manager import write_bc_csv, validate_bc_assignment

    for name in get_physical_names():
        update_bc_from_widgets_dialog(name)

    all_assigned, missing = validate_bc_assignment()
    if not all_assigned:
        logger.warning(
            f"Not all physical names have BCs assigned. Missing: {', '.join(missing)}"
        )
        logger.info("Writing bc.csv anyway with the currently assigned BCs.")

    if write_bc_csv():
        logger.success("Boundary conditions saved to bc.csv")


def load_bc_csv_dialog_callback():
    """
    Load boundary conditions from ``bc.csv`` and refresh the dialog rows.
    """
    from src.utils.gmsh_bc_manager import read_bc_csv

    if read_bc_csv():
        logger.success("Boundary conditions loaded from bc.csv")
        refresh_bc_callback()
    else:
        logger.warning("Could not load bc.csv")


def update_bc_from_widgets_dialog(physical_name: str):
    """
    Read all widget values for *physical_name* and push them into the
    in-memory BC store via :func:`set_boundary_condition`.

    Args:
        physical_name: Physical entity name.
    """
    combo_tag = f"bc_type_{physical_name}"
    if not dpg.does_item_exist(combo_tag):
        return

    bc_type  = dpg.get_value(combo_tag)
    vel_mode = (
        dpg.get_value(f"bc_velmode_{physical_name}")
        if dpg.does_item_exist(f"bc_velmode_{physical_name}")
        else VELOCITY_INPUT_MODES[0]
    )

    var_names = _resolve_variable_names(bc_type, vel_mode)
    variables = {}
    for var_name in var_names:
        tag = f"bc_var_{physical_name}_{var_name}"
        if dpg.does_item_exist(tag):
            variables[var_name] = dpg.get_value(tag)

    set_boundary_condition(physical_name, bc_type, variables)


def _bc_type_has_velocity(bc_type: str) -> bool:
    """
    Return ``True`` when *bc_type* exposes velocity variables that should
    be controlled through the velocity-mode selector.

    Checks whether any variable in :data:`BC_VARIABLES` for this type
    overlaps with the known velocity variable names.  Replace the body with
    an explicit allow-list once your BC_TYPES are finalised, e.g.:
        ``return bc_type in {"inlet", "outlet_velocity"}``

    Args:
        bc_type: BC type string (e.g. ``"inlet"``, ``"wall"``).

    Returns:
        ``True`` if the velocity-mode combo should be shown.
    """
    # TODO: replace with an explicit set once BC_TYPES are finalised
    velocity_var_names = {"u", "v", "w", "v_n", "v_t"}
    vars_for_type = BC_VARIABLES.get(bc_type, [])
    return bool(set(vars_for_type) & velocity_var_names)


def _resolve_variable_names(bc_type: str, velocity_mode: str) -> list[str]:
    """
    Return the ordered list of variable names for the given BC type and
    velocity input mode.

    If *bc_type* has velocity variables, names come from
    :data:`VELOCITY_MODE_VARS`[*velocity_mode*]; otherwise from
    :data:`BC_VARIABLES`[*bc_type*].

    Args:
        bc_type:       BC type string.
        velocity_mode: Velocity input mode string.

    Returns:
        List of variable name strings, e.g. ``["u", "v", "w"]``.
    """
    if _bc_type_has_velocity(bc_type):
        return list(VELOCITY_MODE_VARS.get(velocity_mode, VELOCITY_MODE_VARS[VELOCITY_INPUT_MODES[0]]))
    return list(BC_VARIABLES.get(bc_type, []))


def _infer_velocity_mode(stored_vars: dict) -> str:
    """
    Guess the velocity input mode from a previously stored variable dict.

    Checks for the ``"v_n"`` key to identify normal/tangential mode; falls
    back to the first entry of :data:`VELOCITY_INPUT_MODES`.

    Args:
        stored_vars: Dict of ``{var_name: value}`` from the BC store.

    Returns:
        One of the strings in :data:`VELOCITY_INPUT_MODES`.
    """
    if "v_n" in stored_vars:
        return "Normal + Tangential (n, t)"
    return VELOCITY_INPUT_MODES[0]

def _create_variable_inputs(
    physical_name: str,
    bc_type: str,
    velocity_mode: str,
    current_values: dict,
    parent: str | None = None,
):
    var_names = _resolve_variable_names(bc_type, velocity_mode)

    # Only resolve to UUID when a parent is explicitly given
    parent_id = dpg.get_alias_id(parent) if isinstance(parent, str) else parent

    for var_name in var_names:
        label       = VELOCITY_VAR_LABELS.get(var_name, var_name)
        default_val = float(current_values.get(var_name, 0.0))

        kwargs = {"parent": parent_id} if parent_id is not None else {}

        dpg.add_text(f"{label}=", **kwargs)
        dpg.add_input_float(
            tag=f"bc_var_{physical_name}_{var_name}",
            default_value=default_val,
            width=_COL_VAR_INPUT_W,
            format="%.3f",
            step=0,
            callback=lambda s, a, u=(physical_name, var_name): _on_var_changed(*u, a),
            **kwargs,
        )
        dpg.add_spacer(width=4, **kwargs)
        
def _on_velocity_mode_changed(physical_name: str, new_mode: str):
    bc_type = dpg.get_value(f"bc_type_{physical_name}")
    _rebuild_variable_inputs(physical_name, bc_type, new_mode, current_values={})
    update_bc_from_widgets_dialog(physical_name)


def _rebuild_variable_inputs(
    physical_name: str,
    bc_type: str,
    velocity_mode: str,
    current_values: dict,
):
    container_tag = f"bc_vars_{physical_name}"
    if not dpg.does_item_exist(container_tag):
        return

    for mode_vars in VELOCITY_MODE_VARS.values():
        for var_name in mode_vars:
            tag = f"bc_var_{physical_name}_{var_name}"
            if dpg.does_item_exist(tag):
                print(f"[DEBUG] deleting tag: {tag}")
                dpg.delete_item(tag)
    for var_name in BC_VARIABLES.get(bc_type, []):
        tag = f"bc_var_{physical_name}_{var_name}"
        if dpg.does_item_exist(tag):
            print(f"[DEBUG] deleting BC var tag: {tag}")
            dpg.delete_item(tag)

    dpg.delete_item(container_tag, children_only=True)
    _create_variable_inputs(physical_name, bc_type, velocity_mode, current_values, parent=container_tag)