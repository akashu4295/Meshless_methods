"""
Menu bar construction for MeMPhyS GUI

Creates the application menu bar with File, Settings, and Help menus.
"""

import dearpygui.dearpygui as dpg

from src.callbacks import (
    open_logs_callback,
    open_help_callback,
    show_about_callback,
    show_preferences_callback,
    exit_application_callback,
    open_config_callback,
    save_config_callback,
    show_options_callback,
    load_config_file_callback,
    save_config_file_callback,
    toggle_theme_callback,
)


def create_menu_bar(themes: dict) -> int:
    """
    Create the main menu bar.
    A small spacer is added above the bar so it doesn't butt against the
    OS title bar on Windows / Linux.
    Args: themes: Dictionary of theme IDs
    Returns: Menu bar tag/ID
    """
    
    _create_config_file_dialogs()
    with dpg.menu_bar(tag="main_menu_bar") as menu_bar:
        # ── File ──────────────────────────────────────────────────────────
        with dpg.menu(label="  File  "):
            dpg.add_menu_item(label="Open Configuration",
                callback=open_config_callback)
            dpg.add_menu_item(label="Open Logs",
                callback=open_logs_callback)
            dpg.add_separator()
            dpg.add_menu_item(label="Save Configuration",
                callback=save_config_callback)
            dpg.add_separator()
            dpg.add_menu_item(label="Exit              ",
                callback=exit_application_callback)

        # ── Settings (consolidates Edit + Options + View) ─────────────────
        with dpg.menu(label="  Settings  "):
            dpg.add_menu_item(label="Preferences",
                callback=show_preferences_callback)
            dpg.add_menu_item(label="Application Options",
                callback=show_options_callback)
            dpg.add_separator()
            dpg.add_menu_item(label="Toggle Dark / Light Theme",
                callback=toggle_theme_callback)

        # ── Help ──────────────────────────────────────────────────────────
        with dpg.menu(label="  Help  "):
            dpg.add_menu_item(label="Documentation ",
                callback=open_help_callback)
            dpg.add_separator()
            dpg.add_menu_item(label="About", callback=show_about_callback,)

    if "menu_bar" in themes:
        dpg.bind_item_theme(menu_bar, themes["menu_bar"])

    return menu_bar


def _create_config_file_dialogs():
    """Create file dialogs for configuration load/save."""
    from src.config import FILE_DIALOG_WIDTH, FILE_DIALOG_HEIGHT

    dpg.add_file_dialog(directory_selector=False, tag="file_dialog_config_load",
        callback=load_config_file_callback, show=False,
        width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT,
        default_path="./config")
    dpg.add_file_extension(".json", parent="file_dialog_config_load", color=(255, 255, 0, 255))

    dpg.add_file_dialog(directory_selector=False, tag="file_dialog_config_save",
        callback=save_config_file_callback, show=False,
        width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT,
        default_path="./config", default_filename="my_config.json",)
    dpg.add_file_extension(".json", parent="file_dialog_config_save", color=(255, 255, 0, 255))