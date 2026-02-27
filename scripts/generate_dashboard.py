"""Generate HTML dashboard for the synthesis dataframe."""
import sys
import pandas as pd
import numpy as np
from pathlib import Path

sys.path.insert(0, "src")
from synthesizability.dashboard_plugins import (
    load_plugins,
    collect_summary_cards,
    collect_table_columns,
    collect_detail_sections,
    run_generate,
)


# Columns always shown first in the table regardless of plugins
CORE_COLUMNS = [
    'sample_number', 'sample_id', 'formula', 'has_summary',
    'superconductivity', 'xrd_type', 'xrd_instrument', 'xrd_result',
    'prediction_list', 'mass_loss_percent', 'price_per_gram', 'arc_meltable',
    'disorder_probability',
]

# Core summary cards always shown
CORE_SUMMARY_CARDS = [
    {'label': 'Total Samples', 'value_fn': lambda df: str(len(df))},
    {'label': 'With XRD Data', 'value_fn': lambda df: str((df['xrd_n_files'] > 0).sum())},
    {'label': 'Arc Meltable', 'value_fn': lambda df: str(df['arc_meltable'].sum())},
    {'label': 'Avg Price ($/g)', 'value_fn': lambda df: f"${df['price_per_gram'].mean():.2f}"},
    {'label': 'Avg Mass Loss (%)', 'value_fn': lambda df: f"{df['mass_loss_percent'].mean():.1f}%"},
]

# Sections always shown on detail pages, in order
CORE_DETAIL_SECTIONS = {
    'Basic Information': ['sample_number', 'sample_id', 'formula', 'files', 'has_summary'],
    'Characterization': ['superconductivity', 'tc_kelvin', 'xrd_type', 'xrd_instrument', 'xrd_result'],
    'Synthesis Details': ['synthesis_content', 'mass_loss_percent', 'initial_mass_g', 'final_mass_g'],
    'Cost & Feasibility': ['price_per_gram', 'arc_meltable', 'disorder_probability', 'prediction_list'],
    'XRD Data': ['xrd_patterns', 'xrd_files', 'xrd_n_files', 'xrd_two_theta_min', 'xrd_two_theta_max'],
    'Susceptibility Data': ['chi_files', 'chi_n_files', 'chi_has_high_field', 'chi_fields'],
    'Status': ['status_content'],
}


def serialize_for_table(val):
    """Convert value to compact HTML for table display."""
    if isinstance(val, np.ndarray):
        return f'<span style="color: #666;">[array: {val.shape}]</span>'
    elif isinstance(val, list):
        if len(val) == 0:
            return '<span style="color: #999;">[]</span>'
        return f'[{len(val)} items]'
    elif isinstance(val, bool):
        color = '#28a745' if val else '#dc3545'
        symbol = '✓' if val else '✗'
        return f'<span style="color: {color}; font-weight: bold;">{symbol}</span>'
    elif isinstance(val, (int, float)):
        if pd.isna(val):
            return '<span style="color: #999;">—</span>'
        return f'{val:.4g}' if isinstance(val, float) else str(val)
    elif isinstance(val, str):
        return val[:47] + '...' if len(val) > 50 else val
    try:
        if pd.isna(val):
            return '<span style="color: #999;">—</span>'
    except (ValueError, TypeError):
        pass
    return str(val)


def serialize_for_detail(val, key=None):
    """Convert value to detailed HTML for detail page."""
    if isinstance(val, np.ndarray):
        return f'<span style="color: #666;">[numpy array with shape {val.shape}]</span>'
    elif isinstance(val, list):
        if len(val) == 0:
            return '<em style="color: #999;">Empty list</em>'
        if key == 'xrd_patterns' and len(val) > 0 and isinstance(val[0], dict):
            html = '<div style="margin-top: 10px;">'
            for i, pattern in enumerate(val, 1):
                html += f'<div style="background: #f8f9fa; padding: 10px; margin: 5px 0; border-radius: 4px;">'
                html += f'<strong>Pattern {i}:</strong><br>'
                html += f'Instrument: {pattern.get("instrument", "N/A")}<br>'
                html += f'Date: {pattern.get("date", "N/A")}<br>'
                html += f'Anode: {pattern.get("anode", "N/A")}<br>'
                html += f'2θ range: {pattern.get("two_theta_min", "?")} - {pattern.get("two_theta_max", "?")}°<br>'
                html += f'Points: {pattern.get("n_points", "?")}<br>'
                html += f'File: {pattern.get("filename", "N/A")}'
                html += '</div>'
            html += '</div>'
            return html
        html = '<ul style="margin: 5px 0; padding-left: 20px;">'
        for item in val:
            html += f'<li>{item}</li>'
        html += '</ul>'
        return html
    elif isinstance(val, bool):
        color = '#28a745' if val else '#dc3545'
        text = 'Yes' if val else 'No'
        return f'<span style="color: {color}; font-weight: bold;">{text}</span>'
    elif isinstance(val, (int, float)):
        if pd.isna(val):
            return '<em style="color: #999;">Not available</em>'
        return f'{val:.4g}' if isinstance(val, float) else str(val)
    elif isinstance(val, str):
        return val.replace('\n', '<br>')
    try:
        if pd.isna(val):
            return '<em style="color: #999;">Not available</em>'
    except (ValueError, TypeError):
        pass
    return str(val)


def _build_index_html(df, plugin_cards, plugin_columns):
    """Build the full index page HTML."""

    # Core summary cards
    core_cards_html = ''
    for card in CORE_SUMMARY_CARDS:
        core_cards_html += f'''
            <div class="stat-card">
                <div class="stat-label">{card["label"]}</div>
                <div class="stat-value">{card["value_fn"](df)}</div>
            </div>'''

    # Plugin summary cards
    plugin_cards_html = ''
    for card in plugin_cards:
        plugin_cards_html += f'''
            <div class="stat-card">
                <div class="stat-label">{card["label"]}</div>
                <div class="stat-value">{card["value"]}</div>
            </div>'''

    # Build ordered column list: core first, then plugin columns not already included
    all_columns = [c for c in CORE_COLUMNS if c in df.columns]
    for col in plugin_columns:
        if col not in all_columns and col in df.columns:
            all_columns.append(col)

    # Column toggle checkboxes
    toggles_html = ''
    for i, col in enumerate(all_columns):
        toggles_html += f'<label class="toggle-label"><input type="checkbox" class="column-toggle" checked onchange="toggleColumn({i+1}, this)"> {col}</label>\n'

    # Table headers
    headers_html = '<th>Detail</th>\n'
    for i, col in enumerate(all_columns):
        headers_html += f'<th onclick="sortTable({i+1})">{col} ▾</th>\n'

    # Table rows
    rows_html = ''
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        rows_html += '<tr>\n'
        rows_html += f'<td><a href="samples/{sample_id}.html" class="sample-link">View</a></td>\n'
        for col in all_columns:
            val = serialize_for_table(row[col])
            css_class = ' class="formula"' if col == 'formula' else ''
            rows_html += f'<td{css_class}>{val}</td>\n'
        rows_html += '</tr>\n'

    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Synthesis Dashboard</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0; padding: 20px; background: #f5f5f5; font-size: 14px;
        }}
        .container {{
            max-width: 100%; margin: 0 auto; background: white;
            padding: 30px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        h1 {{ color: #333; border-bottom: 3px solid #0066cc; padding-bottom: 10px; margin-top: 0; }}
        .summary {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px; margin: 20px 0;
        }}
        .stat-card {{
            background: #f8f9fa; padding: 15px; border-radius: 6px;
            border-left: 4px solid #0066cc;
        }}
        .stat-label {{ color: #666; font-size: 0.85em; margin-bottom: 5px; }}
        .stat-value {{ color: #333; font-size: 1.4em; font-weight: bold; }}
        .controls {{ margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 6px; }}
        .search-box {{
            padding: 8px 12px; width: 300px; font-size: 1em;
            border: 2px solid #ddd; border-radius: 4px; margin-bottom: 15px;
        }}
        .search-box:focus {{ outline: none; border-color: #0066cc; }}
        .column-toggles {{
            display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 8px; margin-top: 10px;
        }}
        .toggle-label {{ display: flex; align-items: center; cursor: pointer; padding: 4px; }}
        .toggle-label input {{ margin-right: 6px; }}
        .table-wrapper {{ overflow-x: auto; margin-top: 20px; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; }}
        th {{
            background: #0066cc; color: white; padding: 10px 6px; text-align: left;
            position: sticky; top: 0; z-index: 10; white-space: nowrap; cursor: pointer;
        }}
        th:hover {{ background: #0052a3; }}
        td {{ padding: 8px 6px; border-bottom: 1px solid #e0e0e0; vertical-align: top; }}
        tr:hover {{ background: #f8f9fa; }}
        .sample-link {{ color: #0066cc; text-decoration: none; font-weight: bold; }}
        .sample-link:hover {{ text-decoration: underline; }}
        .formula {{ font-family: 'Courier New', monospace; font-weight: 500; }}
    </style>
    <script>
        let sortDirection = {{}};
        function filterTable() {{
            const filter = document.getElementById('searchBox').value.toLowerCase();
            const rows = document.getElementById('dataTable').getElementsByTagName('tr');
            for (let i = 1; i < rows.length; i++) {{
                rows[i].style.display = rows[i].textContent.toLowerCase().includes(filter) ? '' : 'none';
            }}
        }}
        function toggleColumn(colIndex, checkbox) {{
            const rows = document.getElementById('dataTable').getElementsByTagName('tr');
            for (let row of rows) {{
                if (row.children[colIndex]) row.children[colIndex].style.display = checkbox.checked ? '' : 'none';
            }}
        }}
        function toggleAllColumns(checked) {{
            document.querySelectorAll('.column-toggle').forEach((cb, idx) => {{
                cb.checked = checked;
                toggleColumn(idx + 1, cb);
            }});
        }}
        function sortTable(colIndex) {{
            const table = document.getElementById('dataTable');
            const rows = Array.from(table.tBodies[0].rows);
            if (!sortDirection[colIndex]) sortDirection[colIndex] = 1;
            sortDirection[colIndex] *= -1;
            const dir = sortDirection[colIndex];
            rows.sort((a, b) => {{
                const aT = a.cells[colIndex].textContent.trim();
                const bT = b.cells[colIndex].textContent.trim();
                const aN = parseFloat(aT), bN = parseFloat(bT);
                return (!isNaN(aN) && !isNaN(bN)) ? dir * (aN - bN) : dir * aT.localeCompare(bT);
            }});
            rows.forEach(row => table.tBodies[0].appendChild(row));
        }}
    </script>
</head>
<body>
    <div class="container">
        <h1>Synthesis Data Dashboard</h1>
        <div class="summary">
            {core_cards_html}
            {plugin_cards_html}
        </div>
        <div class="controls">
            <input type="text" id="searchBox" class="search-box"
                   placeholder="Search table..." onkeyup="filterTable()">
            <div style="margin: 10px 0;">
                <button onclick="toggleAllColumns(true)" style="margin-right: 10px;">Show All</button>
                <button onclick="toggleAllColumns(false)">Hide All</button>
            </div>
            <details open>
                <summary style="cursor: pointer; font-weight: bold; margin-bottom: 10px;">Column Visibility</summary>
                <div class="column-toggles">
                    {toggles_html}
                </div>
            </details>
        </div>
        <div class="table-wrapper">
            <table id="dataTable">
                <thead><tr>{headers_html}</tr></thead>
                <tbody>{rows_html}</tbody>
            </table>
        </div>
    </div>
</body>
</html>"""


def _build_detail_html(row, plugin_sections):
    """Build the full detail page HTML for a single sample."""
    sample_id = row['sample_id']
    formula = row['formula']

    # Core detail sections HTML
    core_sections_html = ''
    for section_name, fields in CORE_DETAIL_SECTIONS.items():
        fields_present = [f for f in fields if f in row.index]
        if not fields_present:
            continue
        core_sections_html += f'<div class="section">\n<div class="section-title">{section_name}</div>\n'
        for field in fields_present:
            val = row[field]
            field_html = serialize_for_detail(val, key=field)
            core_sections_html += f'<div class="field">\n<div class="field-label">{field}</div>\n'
            if field in ['synthesis_content', 'status_content'] and isinstance(val, str) and len(val) > 100:
                core_sections_html += f'<div class="text-block">{field_html}</div>\n'
            else:
                core_sections_html += f'<div class="field-value">{field_html}</div>\n'
            core_sections_html += '</div>\n'
        core_sections_html += '</div>\n'

    # Plugin sections HTML
    plugin_sections_html = ''
    for section in plugin_sections:
        plugin_sections_html += f'''
<div class="section">
    <div class="section-title">{section["title"]}</div>
    {section["html"]}
</div>'''

    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>{sample_id} - {formula}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0; padding: 20px; background: #f5f5f5;
        }}
        .container {{
            max-width: 1200px; margin: 0 auto; background: white;
            padding: 30px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .back-link {{ color: #0066cc; text-decoration: none; margin-bottom: 20px; display: inline-block; }}
        .back-link:hover {{ text-decoration: underline; }}
        h1 {{ color: #333; border-bottom: 3px solid #0066cc; padding-bottom: 10px; margin-top: 0; }}
        .formula {{ font-family: 'Courier New', monospace; font-size: 1.2em; color: #666; margin-bottom: 20px; }}
        .section {{ margin: 30px 0; }}
        .section-title {{
            font-size: 1.3em; font-weight: bold; color: #0066cc;
            border-bottom: 2px solid #e0e0e0; padding-bottom: 5px; margin-bottom: 15px;
        }}
        .field {{ margin: 15px 0; }}
        .field-label {{ font-weight: bold; color: #555; margin-bottom: 5px; }}
        .field-value {{ color: #333; line-height: 1.6; }}
        .text-block {{
            background: #f8f9fa; padding: 15px; border-radius: 4px;
            border-left: 4px solid #0066cc; white-space: pre-wrap;
            font-family: 'Courier New', monospace; font-size: 0.9em;
        }}
        .plot-container {{ margin: 20px 0; text-align: center; }}
        .plot-container img {{ max-width: 100%; height: auto; border-radius: 4px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .fit-table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        .fit-table th, .fit-table td {{ padding: 10px; text-align: left; border-bottom: 1px solid #e0e0e0; }}
        .fit-table th {{ background: #f8f9fa; font-weight: bold; }}
        .oqmd-warning {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; border-radius: 4px; margin: 15px 0; }}
        .oqmd-stable {{ background: #d4edda; border-left: 4px solid #28a745; padding: 15px; border-radius: 4px; margin: 15px 0; }}
        .oqmd-unstable {{ background: #f8d7da; border-left: 4px solid #dc3545; padding: 15px; border-radius: 4px; margin: 15px 0; }}
        .cif-link {{
            display: inline-block; margin: 5px 10px 5px 0; padding: 8px 12px;
            background: #0066cc; color: white; text-decoration: none; border-radius: 4px; font-size: 0.9em;
        }}
        .cif-link:hover {{ background: #0052a3; }}
        .external-link {{ color: #0066cc; text-decoration: none; }}
        .external-link:hover {{ text-decoration: underline; }}
    </style>
</head>
<body>
    <div class="container">
        <a href="../index.html" class="back-link">← Back to Dashboard</a>
        <h1>{sample_id}</h1>
        <div class="formula">{formula}</div>
        {core_sections_html}
        {plugin_sections_html}
    </div>
</body>
</html>"""


def main():
    """Generate all dashboard pages."""
    print("=" * 80)
    print("Dashboard Generation")
    print("=" * 80)

    # Load plugins
    print("\nLoading plugins...")
    plugins = load_plugins()
    print(f"  Loaded {len(plugins)} plugins: {[p.__name__.split('.')[-1] for p in plugins]}")

    # Load dataframe
    df = pd.read_pickle('data/processed/synthesis_data.pkl')
    print(f"\nLoaded {len(df)} samples from dataframe")

    # Collect plugin metadata
    plugin_cards = collect_summary_cards(plugins, df)
    plugin_columns = collect_table_columns(plugins, df)

    # Create output directories
    plots_dir = Path('results/dashboard/plots')
    plots_dir.mkdir(parents=True, exist_ok=True)
    samples_dir = Path('results/dashboard/samples')
    samples_dir.mkdir(parents=True, exist_ok=True)
    results_dir = Path('results')

    # Generate index
    print("\nGenerating index page...")
    index_html = _build_index_html(df, plugin_cards, plugin_columns)
    index_path = Path('results/dashboard/index.html')
    index_path.write_text(index_html)
    print(f"  Written to {index_path}")

    # Generate detail pages
    print("\nGenerating detail pages...")
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        print(f"  {sample_id}...", end='', flush=True)

        # Run plugin generation steps (plots etc.)
        run_generate(plugins, row, plots_dir, results_dir)

        # Collect plugin detail sections
        plugin_sections = collect_detail_sections(plugins, row, plots_dir, results_dir)

        # Build and write detail page
        detail_html = _build_detail_html(row, plugin_sections)
        detail_path = samples_dir / f'{sample_id}.html'
        detail_path.write_text(detail_html)
        print(" done")

    print("\n" + "=" * 80)
    print("Dashboard complete!")
    print(f"View at: {Path('results/dashboard/index.html').absolute()}")
    print("=" * 80)


if __name__ == '__main__':
    main()