"""
DearPyGUI theme definitions for MeMPhyS GUI
Contains all visual styling and theme configurations
"""

import dearpygui.dearpygui as dpg

def create_button_theme():
    """
    Create primary button theme with rounded corners and hover effects
    Used for main action buttons like "Compile and Run Solver"
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            # Colors
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))          # Normal state
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))   # Hover state
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))    # Pressed state
            
            # Styling
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)               # Rounded corners
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)             # Inner padding (x, y)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)             # Spacing between items
    
    return theme


def create_secondary_button_theme():
    """
    Create secondary button theme for less prominent actions
    Used for utility buttons like "Browse", "Save", "Clear Logs"
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            # Colors - slightly different from primary
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))
            
            # Styling
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)
    
    return theme


def create_menu_bar_theme():
    """
    Create menu bar theme with dark background
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvMenuBar):
            # Dark menu bar background
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg, (28, 28, 32))
            
            # Light text color
            dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 220, 220))
            
            # Padding
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)
    
    return theme


def create_input_theme():
    """
    Create theme for input fields
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvInputText):
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 8, 4)
    
    return theme


def create_window_theme():
    """
    Create theme for child windows and panels
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvChildWindow):
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (30, 30, 35))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (60, 60, 70))
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 10, 10)
            dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, 6)
    
    return theme


def create_plot_theme():
    """
    Create theme for plots
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvPlot):
            dpg.add_theme_color(dpg.mvPlotCol_FrameBg, (25, 25, 30))
            dpg.add_theme_color(dpg.mvPlotCol_PlotBg, (20, 20, 25))
            dpg.add_theme_color(dpg.mvPlotCol_PlotBorder, (60, 60, 70))
            dpg.add_theme_style(dpg.mvPlotStyleVar_PlotPadding, 12, 12)
    
    return theme


def create_combo_theme():
    """
    Create theme for combo boxes (dropdowns)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvCombo):
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 40, 45))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (50, 50, 55))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 8, 4)
    
    return theme


def create_checkbox_theme():
    """
    Create theme for checkboxes
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvCheckbox):
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark, (40, 120, 200))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 3)
    
    return theme


def create_separator_theme():
    """
    Create theme for separators
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvSeparator):
            dpg.add_theme_color(dpg.mvThemeCol_Separator, (80, 80, 90))
    
    return theme


def create_modal_theme():
    """
    Create theme for modal windows (About, Preferences)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvWindowAppItem):
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (35, 35, 40))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (80, 80, 90))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (28, 28, 32))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (40, 120, 200))
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 8)
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 15, 15)
    
    return theme


def create_disabled_theme():
    """
    Create theme for disabled items
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton, enabled_state=False):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_Text, (120, 120, 125))
    
    return theme


def create_success_button_theme():
    """
    Create theme for success state buttons (green)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 180, 80))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (50, 200, 100))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 160, 70))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
    
    return theme


def create_error_button_theme():
    """
    Create theme for error state buttons (red)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (200, 60, 60))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (220, 80, 80))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (180, 50, 50))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
    
    return theme


def initialize_all_themes():
    """
    Initialize all themes and return them in a dictionary for easy access
    
    Returns:
        dict: Dictionary mapping theme names to theme IDs
    """
    themes = {
        "button": create_button_theme(),
        "button_secondary": create_secondary_button_theme(),
        "button_success": create_success_button_theme(),
        "button_error": create_error_button_theme(),
        "menu_bar": create_menu_bar_theme(),
        "input": create_input_theme(),
        "window": create_window_theme(),
        "plot": create_plot_theme(),
        "combo": create_combo_theme(),
        "checkbox": create_checkbox_theme(),
        "separator": create_separator_theme(),
        "modal": create_modal_theme(),
        "disabled": create_disabled_theme(),
    }
    
    return themes

DARK = {
    "window_bg":        (22,  22,  28),
    "child_bg":         (28,  28,  35),
    "frame_bg":         (38,  38,  48),
    "frame_bg_hovered": (50,  50,  62),
    "frame_bg_active":  (60,  60,  75),
    "menubar_bg":       (18,  18,  24),
    "border":           (55,  55,  70),
    "text":             (215, 215, 220),
    "text_disabled":    (120, 120, 130),
    "button":           (55,  90,  160),
    "button_hovered":   (70,  110, 190),
    "button_active":    (45,  75,  140),
    "header":           (55,  90,  160, 120),
    "header_hovered":   (70,  110, 190, 150),
    "header_active":    (45,  75,  140, 200),
    "tab":              (38,  38,  48),
    "tab_hovered":      (70,  110, 190),
    "tab_active":       (55,  90,  160),
    "title_bg":         (18,  18,  24),
    "title_bg_active":  (30,  55,  120),
    "scrollbar_bg":     (18,  18,  24),
    "scrollbar_grab":   (65,  65,  85),
    "scrollbar_hovered":(85,  85,  105),
    "scrollbar_active": (105, 105, 125),
    "separator":        (55,  55,  70),
    "popup_bg":         (28,  28,  35),
    "check_mark":       (100, 160, 255),
    "slider_grab":      (70,  110, 190),
    "slider_grab_active":(90, 135, 220),
}

LIGHT = {
    "window_bg":        (235, 235, 240),
    "child_bg":         (245, 245, 250),
    "frame_bg":         (255, 255, 255),
    "frame_bg_hovered": (230, 235, 250),
    "frame_bg_active":  (210, 220, 245),
    "menubar_bg":       (220, 220, 228),
    "border":           (185, 185, 200),
    "text":             (25,  25,  30),
    "text_disabled":    (150, 150, 160),
    "button":           (100, 140, 220),
    "button_hovered":   (120, 160, 240),
    "button_active":    (80,  120, 200),
    "header":           (170, 195, 240),        # visible blue-grey tint, full opacity
    "header_hovered":   (145, 175, 230),        # slightly deeper on hover
    "header_active":    (120, 155, 220),        # clearly pressed
    "tab":              (210, 215, 230),
    "tab_hovered":      (120, 160, 240),
    "tab_active":       (100, 140, 220),
    "title_bg":         (210, 215, 230),
    "title_bg_active":  (140, 175, 235),
    "scrollbar_bg":     (215, 215, 225),
    "scrollbar_grab":   (170, 175, 195),
    "scrollbar_hovered":(145, 152, 180),
    "scrollbar_active": (120, 130, 165),
    "separator":        (185, 185, 200),
    "popup_bg":         (240, 240, 245),
    "check_mark":       (60,  100, 210),
    "slider_grab":      (100, 140, 220),
    "slider_grab_active":(80, 120, 200),
}


# ── Theme builder ─────────────────────────────────────────────────────────────

def build_theme(palette: dict):
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvAll):

            # Backgrounds
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg,       palette["window_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg,        palette["child_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg,        palette["popup_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg,        palette["frame_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, palette["frame_bg_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive,  palette["frame_bg_active"])
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg,      palette["menubar_bg"])

            # Title bar
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg,        palette["title_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive,  palette["title_bg_active"])
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgCollapsed, palette["title_bg"])

            # Text
            dpg.add_theme_color(dpg.mvThemeCol_Text,           palette["text"])
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled,   palette["text_disabled"])

            # Borders & separators
            dpg.add_theme_color(dpg.mvThemeCol_Border,         palette["border"])
            dpg.add_theme_color(dpg.mvThemeCol_BorderShadow,   (0, 0, 0, 0))
            dpg.add_theme_color(dpg.mvThemeCol_Separator,      palette["separator"])
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorHovered, palette["button_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorActive,  palette["button_active"])

            # Buttons
            dpg.add_theme_color(dpg.mvThemeCol_Button,         palette["button"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered,  palette["button_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,   palette["button_active"])

            # Headers (collapsibles, selectables, tree nodes)
            dpg.add_theme_color(dpg.mvThemeCol_Header,         palette["header"])
            dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered,  palette["header_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_HeaderActive,   palette["header_active"])

            # Tabs
            dpg.add_theme_color(dpg.mvThemeCol_Tab,            palette["tab"])
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered,     palette["tab_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_TabActive,      palette["tab_active"])
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocused,       palette["tab"])
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocusedActive, palette["tab_active"])

            # Scrollbar
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarBg,    palette["scrollbar_bg"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrab,  palette["scrollbar_grab"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabHovered, palette["scrollbar_hovered"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabActive,  palette["scrollbar_active"])

            # Checkboxes & sliders
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark,      palette["check_mark"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab,     palette["slider_grab"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrabActive, palette["slider_grab_active"])

            # Styles (same for both themes)
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding,      10, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding,        6,  4)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing,          8,  4)
            dpg.add_theme_style(dpg.mvStyleVar_ItemInnerSpacing,     4,  4)
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding,       6)
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding,        4)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarSize,       14)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarRounding,    8)
            dpg.add_theme_style(dpg.mvStyleVar_GrabMinSize,         10)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding,         4)
            dpg.add_theme_style(dpg.mvStyleVar_TabRounding,          4)

    return theme


# ── Public API ────────────────────────────────────────────────────────────────

def apply_dark_theme():
    theme = build_theme(DARK)
    dpg.bind_theme(theme)
    return theme


def apply_light_theme():
    theme = build_theme(LIGHT)
    dpg.bind_theme(theme)
    return theme


def toggle_theme(current_dark: bool):
    if current_dark:
        apply_light_theme()
        return False
    else:
        apply_dark_theme()
        return True