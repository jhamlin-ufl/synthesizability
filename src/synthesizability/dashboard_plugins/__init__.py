"""Dashboard plugin system for synthesizability dashboard."""
import importlib
import pandas as pd
from pathlib import Path


# Registry - comment out to disable, reorder to control section order
PLUGINS = [
    "synthesizability.dashboard_plugins.composition",
    "synthesizability.dashboard_plugins.oqmd",
    "synthesizability.dashboard_plugins.ternary_phases",
    "synthesizability.dashboard_plugins.susceptibility",
    "synthesizability.dashboard_plugins.xrd_rietveld",
    "synthesizability.dashboard_plugins.xrd_comparison",
    "synthesizability.dashboard_plugins.supercon",
]


def load_plugins():
    """Load all registered plugins."""
    plugins = []
    for module_path in PLUGINS:
        try:
            module = importlib.import_module(module_path)
            plugins.append(module)
        except ImportError as e:
            print(f"  WARNING: Could not load plugin {module_path}: {e}")
    return plugins


def collect_summary_cards(plugins, df):
    """Collect summary cards from all plugins."""
    cards = []
    for plugin in plugins:
        try:
            cards.extend(plugin.get_summary_cards(df))
        except Exception as e:
            print(f"  WARNING: {plugin.__name__}.get_summary_cards failed: {e}")
    return cards


def collect_table_columns(plugins, df):
    """Collect ordered list of column names from all plugins."""
    columns = []
    for plugin in plugins:
        try:
            for col in plugin.get_table_columns(df):
                if col in df.columns and col not in columns:
                    columns.append(col)
        except Exception as e:
            print(f"  WARNING: {plugin.__name__}.get_table_columns failed: {e}")
    return columns


def collect_detail_sections(plugins, row, plots_dir, results_dir):
    """Collect detail sections from all plugins for a single sample."""
    sections = []
    for plugin in plugins:
        try:
            section = plugin.get_detail_section(row, plots_dir, results_dir)
            if section is not None:
                sections.append(section)
        except Exception as e:
            print(f"  WARNING: {plugin.__name__}.get_detail_section failed for "
                  f"{row['sample_id']}: {e}")
    return sections


def run_generate(plugins, row, plots_dir, results_dir):
    """Run generate() for all plugins for a single sample."""
    for plugin in plugins:
        try:
            plugin.generate(row, plots_dir, results_dir)
        except Exception as e:
            print(f"  WARNING: {plugin.__name__}.generate failed for "
                  f"{row['sample_id']}: {e}")
