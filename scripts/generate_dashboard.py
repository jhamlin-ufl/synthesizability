"""Generate an HTML dashboard for the synthesis dataframe."""
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

sys.path.insert(0, "src")
from synthesizability.susceptibility import (
    load_all_chi_data,
    extract_tc_values,
    plot_single_chi_real,
    plot_single_chi_imaginary,
    plot_single_hc2
)


def serialize_for_table(val):
    """Convert value to compact HTML representation for table display."""
    # Check type before pd.isna() to avoid numpy array issues
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
        if len(val) > 50:
            return val[:47] + '...'
        return val
    # Fallback for any other types including None
    try:
        if pd.isna(val):
            return '<span style="color: #999;">—</span>'
    except (ValueError, TypeError):
        pass
    return str(val)


def serialize_for_detail(val, key=None):
    """Convert value to detailed HTML representation for detail page."""
    # Check types before pd.isna()
    if isinstance(val, np.ndarray):
        return f'<span style="color: #666;">[numpy array with shape {val.shape}]</span>'
    elif isinstance(val, list):
        if len(val) == 0:
            return '<em style="color: #999;">Empty list</em>'
        # Special handling for xrd_patterns
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
        # Regular lists
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
    # Fallback
    try:
        if pd.isna(val):
            return '<em style="color: #999;">Not available</em>'
    except (ValueError, TypeError):
        pass
    return str(val)


def extract_composition(sample_id: str) -> str:
    """Extract composition from sample_id (e.g., '0449_HM_HfTa4Zr' -> 'HfTa4Zr')."""
    parts = sample_id.split('_')
    return parts[2] if len(parts) > 2 else sample_id


def generate_susceptibility_plots(sample_id: str, sample_dir: Path, 
                                  plots_dir: Path) -> dict:
    """
    Generate susceptibility plots for a sample.
    
    Returns:
        dict with keys 'chi_real', 'chi_imag', 'hc2' containing relative paths
        to plot files, and 'fit_results' containing fit parameters
    """
    composition = extract_composition(sample_id)
    
    # Load chi data
    chi_data = load_all_chi_data(sample_dir)
    
    result = {
        'chi_real': None,
        'chi_imag': None,
        'hc2': None,
        'fit_results': None
    }
    
    if len(chi_data) == 0:
        return result
    
    # Generate chi real plot
    fig_real = plot_single_chi_real(chi_data, composition, sample_id)
    real_path = plots_dir / f'{sample_id}_chi_real.png'
    fig_real.savefig(real_path, dpi=150, bbox_inches='tight')
    plt.close(fig_real)
    result['chi_real'] = f'../plots/{real_path.name}'
    
    # Generate chi imaginary plot
    fig_imag = plot_single_chi_imaginary(chi_data, composition, sample_id)
    imag_path = plots_dir / f'{sample_id}_chi_imag.png'
    fig_imag.savefig(imag_path, dpi=150, bbox_inches='tight')
    plt.close(fig_imag)
    result['chi_imag'] = f'../plots/{imag_path.name}'
    
    # Extract Tc values and generate Hc2 plot
    tc_data = extract_tc_values(chi_data)
    if len(tc_data) > 0:
        fig_hc2, fit_results = plot_single_hc2(tc_data, composition)
        hc2_path = plots_dir / f'{sample_id}_hc2.png'
        fig_hc2.savefig(hc2_path, dpi=150, bbox_inches='tight')
        plt.close(fig_hc2)
        result['hc2'] = f'../plots/{hc2_path.name}'
        result['fit_results'] = fit_results
    
    return result


def generate_index(df, fit_params_df, output_path):
    """Generate main dashboard index page."""
    
    # Merge fit parameters into dataframe
    if fit_params_df is not None:
        # Create composition column for merging
        df_with_comp = df.copy()
        df_with_comp['Composition'] = df_with_comp['sample_id'].apply(extract_composition)
        
        # Merge
        df_merged = df_with_comp.merge(
            fit_params_df,
            on='Composition',
            how='left'
        )
    else:
        df_merged = df
    
    # Calculate summary statistics
    n_samples = len(df_merged)
    n_superconducting = df_merged['tc_kelvin'].notna().sum()
    n_with_xrd = (df_merged['xrd_n_files'] > 0).sum()
    n_with_chi = (df_merged['chi_n_files'] > 0).sum()
    n_arc_meltable = df_merged['arc_meltable'].sum()
    avg_price = df_merged['price_per_gram'].mean()
    avg_mass_loss = df_merged['mass_loss_percent'].mean()
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Synthesis Dashboard</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
            font-size: 14px;
        }}
        .container {{
            max-width: 100%;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #0066cc;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #0066cc;
        }}
        .stat-label {{
            color: #666;
            font-size: 0.85em;
            margin-bottom: 5px;
        }}
        .stat-value {{
            color: #333;
            font-size: 1.4em;
            font-weight: bold;
        }}
        .controls {{
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
        }}
        .search-box {{
            padding: 8px 12px;
            width: 300px;
            font-size: 1em;
            border: 2px solid #ddd;
            border-radius: 4px;
            margin-bottom: 15px;
        }}
        .search-box:focus {{
            outline: none;
            border-color: #0066cc;
        }}
        .column-toggles {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 8px;
            margin-top: 10px;
        }}
        .toggle-label {{
            display: flex;
            align-items: center;
            cursor: pointer;
            padding: 4px;
        }}
        .toggle-label input {{
            margin-right: 6px;
        }}
        .table-wrapper {{
            overflow-x: auto;
            margin-top: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.85em;
        }}
        th {{
            background: #0066cc;
            color: white;
            padding: 10px 6px;
            text-align: left;
            position: sticky;
            top: 0;
            z-index: 10;
            white-space: nowrap;
            cursor: pointer;
        }}
        th:hover {{
            background: #0052a3;
        }}
        td {{
            padding: 8px 6px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .sample-link {{
            color: #0066cc;
            text-decoration: none;
            font-weight: bold;
        }}
        .sample-link:hover {{
            text-decoration: underline;
        }}
        .formula {{
            font-family: 'Courier New', monospace;
            font-weight: 500;
        }}
    </style>
    <script>
        let sortDirection = {{}};
        
        function filterTable() {{
            const input = document.getElementById('searchBox');
            const filter = input.value.toLowerCase();
            const table = document.getElementById('dataTable');
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {{
                const row = rows[i];
                const text = row.textContent.toLowerCase();
                row.style.display = text.includes(filter) ? '' : 'none';
            }}
        }}
        
        function toggleColumn(colIndex, checkbox) {{
            const table = document.getElementById('dataTable');
            const rows = table.getElementsByTagName('tr');
            
            for (let row of rows) {{
                const cell = row.children[colIndex];
                if (cell) {{
                    cell.style.display = checkbox.checked ? '' : 'none';
                }}
            }}
        }}
        
        function toggleAllColumns(checked) {{
            const checkboxes = document.querySelectorAll('.column-toggle');
            checkboxes.forEach((cb, idx) => {{
                cb.checked = checked;
                toggleColumn(idx + 1, cb);
            }});
        }}
        
        function sortTable(colIndex) {{
            const table = document.getElementById('dataTable');
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            
            // Toggle sort direction
            if (!sortDirection[colIndex]) sortDirection[colIndex] = 1;
            sortDirection[colIndex] *= -1;
            const direction = sortDirection[colIndex];
            
            rows.sort((a, b) => {{
                const aText = a.cells[colIndex].textContent.trim();
                const bText = b.cells[colIndex].textContent.trim();
                
                // Try numeric comparison first
                const aNum = parseFloat(aText);
                const bNum = parseFloat(bText);
                
                if (!isNaN(aNum) && !isNaN(bNum)) {{
                    return direction * (aNum - bNum);
                }}
                
                // Fall back to string comparison
                return direction * aText.localeCompare(bText);
            }});
            
            rows.forEach(row => tbody.appendChild(row));
        }}
    </script>
</head>
<body>
    <div class="container">
        <h1>Synthesis Data Dashboard</h1>
        
        <div class="summary">
            <div class="stat-card">
                <div class="stat-label">Total Samples</div>
                <div class="stat-value">{n_samples}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Superconducting</div>
                <div class="stat-value">{n_superconducting}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">With XRD Data</div>
                <div class="stat-value">{n_with_xrd}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">With χ Data</div>
                <div class="stat-value">{n_with_chi}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Arc Meltable</div>
                <div class="stat-value">{n_arc_meltable}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Price ($/g)</div>
                <div class="stat-value">${avg_price:.2f}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Mass Loss (%)</div>
                <div class="stat-value">{avg_mass_loss:.1f}%</div>
            </div>
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
"""
    
    # Add column toggle checkboxes
    for i, col in enumerate(df_merged.columns):
        html += f'                    <label class="toggle-label"><input type="checkbox" class="column-toggle" checked onchange="toggleColumn({i+1}, this)"> {col}</label>\n'
    
    html += """                </div>
            </details>
        </div>
        
        <div class="table-wrapper">
            <table id="dataTable">
                <thead>
                    <tr>
                        <th>Detail</th>
"""
    
    # Add column headers with sorting
    for i, col in enumerate(df_merged.columns):
        html += f'                        <th onclick="sortTable({i+1})">{col} ▾</th>\n'
    
    html += """                    </tr>
                </thead>
                <tbody>
"""
    
    # Add data rows
    for idx, row in df_merged.iterrows():
        sample_id = row['sample_id']
        html += '                    <tr>\n'
        html += f'                        <td><a href="samples/{sample_id}.html" class="sample-link">View</a></td>\n'
        
        for col in df_merged.columns:
            val = serialize_for_table(row[col])
            css_class = ' class="formula"' if col == 'formula' else ''
            html += f'                        <td{css_class}>{val}</td>\n'
        
        html += '                    </tr>\n'
    
    html += """                </tbody>
            </table>
        </div>
    </div>
</body>
</html>
"""
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    print(f"Index written to {output_path}")


def generate_detail_page(row, output_path, plot_info: dict):
    """Generate detail page for a single sample."""
    
    sample_id = row['sample_id']
    formula = row['formula']
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>{sample_id} - {formula}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .back-link {{
            color: #0066cc;
            text-decoration: none;
            margin-bottom: 20px;
            display: inline-block;
        }}
        .back-link:hover {{
            text-decoration: underline;
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #0066cc;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .formula {{
            font-family: 'Courier New', monospace;
            font-size: 1.2em;
            color: #666;
            margin-bottom: 20px;
        }}
        .section {{
            margin: 30px 0;
        }}
        .section-title {{
            font-size: 1.3em;
            font-weight: bold;
            color: #0066cc;
            border-bottom: 2px solid #e0e0e0;
            padding-bottom: 5px;
            margin-bottom: 15px;
        }}
        .field {{
            margin: 15px 0;
        }}
        .field-label {{
            font-weight: bold;
            color: #555;
            margin-bottom: 5px;
        }}
        .field-value {{
            color: #333;
            line-height: 1.6;
        }}
        .text-block {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 4px;
            border-left: 4px solid #0066cc;
            white-space: pre-wrap;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }}
        .plot-container {{
            margin: 20px 0;
            text-align: center;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .fit-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        .fit-table th, .fit-table td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        .fit-table th {{
            background: #f8f9fa;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="container">
        <a href="../index.html" class="back-link">← Back to Dashboard</a>
        
        <h1>{sample_id}</h1>
        <div class="formula">{formula}</div>
"""
    
    # Group fields into sections
    sections = {
        'Basic Information': ['sample_number', 'sample_id', 'formula', 'files', 'has_summary'],
        'Characterization': ['superconductivity', 'tc_kelvin', 'xrd_type', 'xrd_instrument', 'xrd_result'],
        'Synthesis Details': ['synthesis_content', 'mass_loss_percent', 'initial_mass_g', 'final_mass_g'],
        'Cost & Feasibility': ['price_per_gram', 'arc_meltable', 'prediction_list'],
        'XRD Data': ['xrd_patterns', 'xrd_files', 'xrd_n_files', 'xrd_two_theta_min', 'xrd_two_theta_max'],
        'Susceptibility Data': ['chi_files', 'chi_n_files', 'chi_has_high_field', 'chi_fields'],
        'Status': ['status_content']
    }
    
    for section_name, fields in sections.items():
        html += f'        <div class="section">\n'
        html += f'            <div class="section-title">{section_name}</div>\n'
        
        for field in fields:
            if field not in row.index:
                continue
                
            val = row[field]
            field_html = serialize_for_detail(val, key=field)
            
            html += f'            <div class="field">\n'
            html += f'                <div class="field-label">{field}</div>\n'
            
            # Use text-block styling for long text fields
            if field in ['synthesis_content', 'status_content'] and isinstance(val, str) and len(val) > 100:
                html += f'                <div class="text-block">{field_html}</div>\n'
            else:
                html += f'                <div class="field-value">{field_html}</div>\n'
            
            html += f'            </div>\n'
        
        html += f'        </div>\n'
    
    # Add susceptibility plots section if available
    if plot_info['chi_real'] or plot_info['chi_imag'] or plot_info['hc2']:
        html += '        <div class="section">\n'
        html += '            <div class="section-title">Susceptibility Analysis</div>\n'
        
        if plot_info['chi_real']:
            html += '            <div class="plot-container">\n'
            html += f'                <img src="{plot_info["chi_real"]}" alt="Real Susceptibility">\n'
            html += '            </div>\n'
        
        if plot_info['chi_imag']:
            html += '            <div class="plot-container">\n'
            html += f'                <img src="{plot_info["chi_imag"]}" alt="Imaginary Susceptibility">\n'
            html += '            </div>\n'
        
        if plot_info['hc2']:
            html += '            <div class="plot-container">\n'
            html += f'                <img src="{plot_info["hc2"]}" alt="Upper Critical Field">\n'
            html += '            </div>\n'
        
        # Add fit parameters table if available
        if plot_info['fit_results']:
            fit = plot_info['fit_results']
            html += '            <div class="field">\n'
            html += '                <div class="field-label">Fit Parameters</div>\n'
            html += '                <table class="fit-table">\n'
            html += '                    <thead><tr><th>Model</th><th>Hc2(0) (T)</th><th>Tc (K)</th></tr></thead>\n'
            html += '                    <tbody>\n'
            html += f'                        <tr><td>Linear</td><td>{fit["linear"]["Hc2_0"]:.3f}</td><td>{fit["linear"]["Tc"]:.3f}</td></tr>\n'
            html += f'                        <tr><td>Quadratic</td><td>{fit["quadratic"]["Hc2_0"]:.3f}</td><td>{fit["quadratic"]["Tc"]:.3f}</td></tr>\n'
            html += '                    </tbody>\n'
            html += '                </table>\n'
            html += '            </div>\n'
        
        html += '        </div>\n'
    
    html += """    </div>
</body>
</html>
"""
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)


def main():
    """Generate all dashboard pages."""
    print("="*80)
    print("Dashboard Generation")
    print("="*80)
    
    # Load dataframe
    df = pd.read_pickle('data/processed/synthesis_data.pkl')
    print(f"\nLoaded {len(df)} samples from dataframe")
    
    # Load fit parameters if available
    fit_params_path = Path('results/susceptibility/hc2_fit_parameters.csv')
    if fit_params_path.exists():
        fit_params_df = pd.read_csv(fit_params_path)
        print(f"Loaded fit parameters for {len(fit_params_df)} samples")
    else:
        fit_params_df = None
        print("No fit parameters file found - susceptibility columns will be empty")
    
    # Create output directories
    plots_dir = Path('results/dashboard/plots')
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    samples_dir = Path('results/dashboard/samples')
    samples_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate index
    print("\nGenerating index page...")
    index_path = Path('results/dashboard/index.html')
    generate_index(df, fit_params_df, index_path)
    
    # Generate detail pages with plots
    print("\nGenerating detail pages with susceptibility plots...")
    data_dir = Path('data/raw')
    
    for idx, row in df.iterrows():
        sample_id = row['sample_id']
        sample_dir = data_dir / sample_id
        
        print(f"  Processing {sample_id}...", end='')
        
        # Generate plots if chi data exists
        if row['chi_n_files'] > 0:
            plot_info = generate_susceptibility_plots(sample_id, sample_dir, plots_dir)
            print(f" [χ plots generated]")
        else:
            plot_info = {
                'chi_real': None,
                'chi_imag': None,
                'hc2': None,
                'fit_results': None
            }
            print(f" [no χ data]")
        
        # Generate detail page
        detail_path = samples_dir / f'{sample_id}.html'
        generate_detail_page(row, detail_path, plot_info)
    
    print("\n" + "="*80)
    print(f"Dashboard complete!")
    print(f"View at: {index_path.absolute()}")
    print("="*80)


if __name__ == '__main__':
    main()