# src/synthesizability/susceptibility.py
"""
AC susceptibility analysis functions.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit


def load_chi_file(filepath: Path) -> pd.DataFrame:
    """
    Load a single chi file.
    
    Returns:
        DataFrame with columns: temperature, lockin_v, chan2_v, field, timestamp
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find data start
    data_start_idx = None
    for i, line in enumerate(lines):
        if 'Lockin_V' in line or 'Chan2_V' in line:
            data_start_idx = i + 1
            break
    
    if data_start_idx is None:
        raise ValueError(f"Could not find data header in {filepath}")
    
    data = pd.read_csv(
        filepath,
        sep='\t',
        skiprows=data_start_idx,
        names=['lockin_v', 'chan2_v', 'temperature', 'field', 'timestamp'],
        engine='python'
    )
    
    return data


def extract_field_from_chi_filename(filename: str) -> float:
    """
    Extract field value from chi filename.
    
    Example: 'HM449_chiAC_vs_T_B_1.0T.txt' -> 1.0
    """
    match = re.search(r'B_([0-9.]+)T', filename)
    if match:
        return float(match.group(1))
    return None


def load_all_chi_data(sample_dir: Path, max_temp: float = 10.0) -> pd.DataFrame:
    """
    Load all chi files for a sample.
    
    Args:
        sample_dir: Sample directory containing chi files
        max_temp: Maximum temperature to include
        
    Returns:
        DataFrame with all chi data for this sample
    """
    chi_files = [f for f in sample_dir.iterdir()
                 if f.suffix == '.txt' and 'chi' in f.name.lower()]
    
    if not chi_files:
        return pd.DataFrame()
    
    all_data = []
    
    for chi_file in chi_files:
        field = extract_field_from_chi_filename(chi_file.name)
        if field is None:
            continue
            
        try:
            data = load_chi_file(chi_file)
            data = data[data['temperature'] <= max_temp]
            data['field_tesla'] = field
            data['filename'] = chi_file.name
            all_data.append(data)
        except Exception as e:
            print(f"Warning: Could not load {chi_file.name}: {e}")
            continue
    
    if not all_data:
        return pd.DataFrame()
    
    return pd.concat(all_data, ignore_index=True)


def extract_tc_from_chi_imaginary(field_data: pd.DataFrame,
                                   window_length: int = 31,
                                   polyorder: int = 3) -> float:
    """
    Extract Tc from minimum in imaginary part of susceptibility.
    
    Args:
        field_data: Data for a single field value
        window_length: Savitzky-Golay filter window
        polyorder: Savitzky-Golay polynomial order
        
    Returns:
        Tc in Kelvin, or None if cannot be determined
    """
    if len(field_data) < window_length:
        return None
    
    field_data = field_data.sort_values('temperature')
    chi_imag = field_data['chan2_v'].values
    temp = field_data['temperature'].values
    
    chi_imag_smoothed = savgol_filter(chi_imag,
                                      window_length=window_length,
                                      polyorder=polyorder)
    
    min_idx = np.argmin(chi_imag_smoothed)
    tc = temp[min_idx]
    
    return tc


def extract_tc_values(chi_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract Tc vs field for all field values in chi data.
    
    Returns:
        DataFrame with columns: field_tesla, tc_kelvin
    """
    tc_records = []
    
    for field in sorted(chi_data['field_tesla'].unique()):
        field_data = chi_data[chi_data['field_tesla'] == field]
        
        tc = extract_tc_from_chi_imaginary(field_data)
        
        if tc is not None:
            tc_records.append({
                'field_tesla': field,
                'tc_kelvin': tc
            })
    
    return pd.DataFrame(tc_records)


def linear_hc2_model(T, Hc2_0, Tc):
    """Linear Hc2 model: H(T) = Hc2(0) * (1 - T/Tc)"""
    return Hc2_0 * (1 - T / Tc)


def quadratic_hc2_model(T, Hc2_0, Tc):
    """Quadratic Hc2 model: H(T) = Hc2(0) * (1 - (T/Tc)^2)"""
    return Hc2_0 * (1 - (T / Tc)**2)


def fit_hc2_models(tc_data: pd.DataFrame) -> dict:
    """
    Fit linear and quadratic Hc2 models.
    
    Args:
        tc_data: DataFrame with columns field_tesla, tc_kelvin
        
    Returns:
        dict with 'linear' and 'quadratic' fit parameters
    """
    if len(tc_data) < 2:
        return None
    
    T = tc_data['tc_kelvin'].values
    H = tc_data['field_tesla'].values
    
    Tc_guess = T.max()
    Hc2_0_guess = H.max() * 1.5
    
    try:
        popt_linear, _ = curve_fit(linear_hc2_model, T, H,
                                   p0=[Hc2_0_guess, Tc_guess])
        popt_quadratic, _ = curve_fit(quadratic_hc2_model, T, H,
                                      p0=[Hc2_0_guess, Tc_guess])
        
        return {
            'linear': {'Hc2_0': popt_linear[0], 'Tc': popt_linear[1]},
            'quadratic': {'Hc2_0': popt_quadratic[0], 'Tc': popt_quadratic[1]}
        }
    except:
        return None


def composition_to_latex(composition: str) -> str:
    """Convert composition string to LaTeX chemical formula."""
    result = ""
    i = 0
    while i < len(composition):
        if i < len(composition) and composition[i].isupper():
            element = composition[i]
            i += 1
            if i < len(composition) and composition[i].islower():
                element += composition[i]
                i += 1
            
            subscript = ""
            while i < len(composition) and composition[i].isdigit():
                subscript += composition[i]
                i += 1
            
            if subscript:
                result += f"{element}$_{{{subscript}}}$"
            else:
                result += element
    
    return result


def trim_high_temp_outliers(field_data: pd.DataFrame, n_trim: int = 3) -> pd.DataFrame:
    """Remove last n_trim data points to eliminate high-temp outliers."""
    if len(field_data) > n_trim:
        return field_data.iloc[:-n_trim]
    return field_data


def plot_chi_real_grid(samples_data: dict, sample_order: list,
                       figsize=(16, 16)) -> plt.Figure:
    """
    Plot real part of susceptibility in 4x4 grid.
    
    Args:
        samples_data: dict mapping sample_id to chi DataFrame
        sample_order: list of sample_ids in desired order (max 15)
        
    Returns:
        matplotlib Figure
    """
    fig, axes = plt.subplots(4, 4, figsize=figsize)
    axes = axes.flatten()
    
    for idx, sample_id in enumerate(sample_order[:15]):
        ax = axes[idx]
        chi_data = samples_data.get(sample_id)
        
        if chi_data is None or len(chi_data) == 0:
            ax.axis('off')
            continue
        
        # Extract composition from sample_id (e.g., '0449_HM_HfTa4Zr' -> 'HfTa4Zr')
        parts = sample_id.split('_')
        composition = parts[2] if len(parts) > 2 else sample_id
        latex_title = composition_to_latex(composition)
        
        # Determine if needs trimming (samples with known issues)
        needs_trimming = sample_id in ['0455_HM_NbTa2Zr', '0457_HM_MoNbZr2']
        
        for field in sorted(chi_data['field_tesla'].unique(), reverse=True):
            field_data = chi_data[chi_data['field_tesla'] == field].sort_values('temperature')
            
            if len(field_data) == 0:
                continue
            
            if needs_trimming and abs(field - 2.0) < 0.1:
                field_data = trim_high_temp_outliers(field_data, n_trim=3)
            
            # Normalize to zero at highest temperature
            high_temp_value = field_data.iloc[-1]['lockin_v']
            normalized = (field_data['lockin_v'] - high_temp_value) * 1e6
            
            ax.plot(field_data['temperature'], normalized,
                   label=f'{field:.1f} T', alpha=0.7, linewidth=1.5)
        
        ax.set_xlabel('Temperature (K)', fontsize=9)
        ax.set_ylabel(r"$\chi'$ offset ($\mu$V)", fontsize=9)
        ax.set_title(latex_title, fontsize=10, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=8)
    
    # Legend in the last subplot (position 15)
    axes[15].axis('off')
    if len(sample_order) > 0:
        handles, labels = axes[0].get_legend_handles_labels()
        axes[15].legend(handles[::-1], labels[::-1],
                       loc='center', fontsize=12, frameon=False)
    
    # Turn off any remaining empty subplots
    for idx in range(len(sample_order), 15):
        axes[idx].axis('off')
    
    plt.tight_layout()
    return fig


def plot_chi_imaginary_grid(samples_data: dict, sample_order: list,
                            window_length: int = 31, polyorder: int = 3,
                            figsize=(16, 16)) -> plt.Figure:
    """
    Plot imaginary part of susceptibility with smoothing in 4x4 grid.
    
    Args:
        samples_data: dict mapping sample_id to chi DataFrame
        sample_order: list of sample_ids in desired order (max 15)
        window_length: Savitzky-Golay filter window
        polyorder: Savitzky-Golay polynomial order
        
    Returns:
        matplotlib Figure
    """
    fig, axes = plt.subplots(4, 4, figsize=figsize)
    axes = axes.flatten()
    
    for idx, sample_id in enumerate(sample_order[:15]):
        ax = axes[idx]
        chi_data = samples_data.get(sample_id)
        
        if chi_data is None or len(chi_data) == 0:
            ax.axis('off')
            continue
        
        parts = sample_id.split('_')
        composition = parts[2] if len(parts) > 2 else sample_id
        latex_title = composition_to_latex(composition)
        
        needs_trimming = sample_id in ['0455_HM_NbTa2Zr', '0457_HM_MoNbZr2']
        
        for field in sorted(chi_data['field_tesla'].unique(), reverse=True):
            field_data = chi_data[chi_data['field_tesla'] == field].sort_values('temperature')
            
            if len(field_data) == 0:
                continue
            
            if needs_trimming and abs(field - 2.0) < 0.1:
                field_data = trim_high_temp_outliers(field_data, n_trim=3)
            
            if len(field_data) < window_length:
                continue
            
            high_temp_value = field_data.iloc[-1]['chan2_v']
            normalized = (field_data['chan2_v'] - high_temp_value) * 1e6
            smoothed = savgol_filter(normalized,
                                    window_length=window_length,
                                    polyorder=polyorder)
            
            ax.plot(field_data['temperature'], smoothed,
                   label=f'{field:.1f} T', alpha=0.7, linewidth=1.5)
        
        ax.set_xlabel('Temperature (K)', fontsize=9)
        ax.set_ylabel(r"$\chi''$ offset ($\mu$V)", fontsize=9)
        ax.set_title(latex_title, fontsize=10, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=8)
    
    # Legend in the last subplot
    axes[15].axis('off')
    if len(sample_order) > 0:
        handles, labels = axes[0].get_legend_handles_labels()
        axes[15].legend(handles[::-1], labels[::-1],
                       loc='center', fontsize=12, frameon=False)
    
    for idx in range(len(sample_order), 15):
        axes[idx].axis('off')
    
    plt.tight_layout()
    return fig


def plot_hc2_grid(samples_tc_data: dict, sample_order: list,
                  figsize=(16, 16)) -> tuple:
    """
    Plot Hc2(T) with fits in 4x4 grid.
    
    Args:
        samples_tc_data: dict mapping sample_id to Tc DataFrame
        sample_order: list of sample_ids in desired order (max 15)
        
    Returns:
        tuple of (figure, fit_results_dict)
    """
    fig, axes = plt.subplots(4, 4, figsize=figsize)
    axes = axes.flatten()
    
    fit_results_all = {}
    
    for idx, sample_id in enumerate(sample_order[:15]):
        ax = axes[idx]
        tc_data = samples_tc_data.get(sample_id)
        
        if tc_data is None or len(tc_data) == 0:
            ax.axis('off')
            continue
        
        parts = sample_id.split('_')
        composition = parts[2] if len(parts) > 2 else sample_id
        latex_title = composition_to_latex(composition)
        
        fit_results = fit_hc2_models(tc_data)
        
        if fit_results is None:
            ax.axis('off')
            continue
        
        fit_results_all[composition] = fit_results
        
        # Determine plot ranges
        max_tc = np.ceil(tc_data['tc_kelvin'].max())
        max_field_data = tc_data['field_tesla'].max()
        max_field_fit = max(fit_results['linear']['Hc2_0'],
                           fit_results['quadratic']['Hc2_0'])
        max_field = np.ceil(max(max_field_data, max_field_fit))
        
        # Plot data
        ax.plot(tc_data['tc_kelvin'], tc_data['field_tesla'],
               'o', markersize=8, color='black', label='Data')
        
        # Plot fits
        T_fit_linear = np.linspace(0, fit_results['linear']['Tc'], 200)
        H_linear = linear_hc2_model(T_fit_linear,
                                    fit_results['linear']['Hc2_0'],
                                    fit_results['linear']['Tc'])
        H_linear = np.clip(H_linear, 0, None)
        ax.plot(T_fit_linear, H_linear, '-', linewidth=2,
               color='C0', label='Linear', alpha=0.7)
        
        T_fit_quadratic = np.linspace(0, fit_results['quadratic']['Tc'], 200)
        H_quadratic = quadratic_hc2_model(T_fit_quadratic,
                                         fit_results['quadratic']['Hc2_0'],
                                         fit_results['quadratic']['Tc'])
        H_quadratic = np.clip(H_quadratic, 0, None)
        ax.plot(T_fit_quadratic, H_quadratic, '--', linewidth=2,
               color='C1', label='Quadratic', alpha=0.7)
        
        ax.set_xlim(0, max_tc)
        ax.set_ylim(0, max_field)
        ax.set_xlabel('$T_c$ (K)', fontsize=9)
        ax.set_ylabel('Field (T)', fontsize=9)
        ax.set_title(latex_title, fontsize=10, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=8)
    
    # Legend
    axes[15].axis('off')
    if len(sample_order) > 0 and len(fit_results_all) > 0:
        handles, labels = axes[0].get_legend_handles_labels()
        axes[15].legend(handles, labels,
                       loc='center', fontsize=12, frameon=False)
    
    for idx in range(len(sample_order), 15):
        axes[idx].axis('off')
    
    plt.tight_layout()
    return fig, fit_results_all


# Single-sample plotting functions for dashboard

def plot_single_chi_real(chi_data: pd.DataFrame, composition: str,
                         sample_id: str = None,
                         figsize=(8, 6)) -> plt.Figure:
    """
    Plot real part of susceptibility for a single sample.
    
    Args:
        chi_data: DataFrame with chi data for one sample
        composition: Chemical composition string
        sample_id: Optional sample_id for trimming logic
        figsize: Figure size
        
    Returns:
        matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if chi_data is None or len(chi_data) == 0:
        ax.text(0.5, 0.5, 'No susceptibility data available',
               ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig
    
    latex_title = composition_to_latex(composition)
    needs_trimming = sample_id in ['0455_HM_NbTa2Zr', '0457_HM_MoNbZr2']
    
    for field in sorted(chi_data['field_tesla'].unique(), reverse=True):
        field_data = chi_data[chi_data['field_tesla'] == field].sort_values('temperature')
        
        if len(field_data) == 0:
            continue
        
        if needs_trimming and abs(field - 2.0) < 0.1:
            field_data = trim_high_temp_outliers(field_data, n_trim=3)
        
        high_temp_value = field_data.iloc[-1]['lockin_v']
        normalized = (field_data['lockin_v'] - high_temp_value) * 1e6
        
        ax.plot(field_data['temperature'], normalized,
               label=f'{field:.1f} T', alpha=0.7, linewidth=2)
    
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel(r"$\chi'$ offset ($\mu$V)", fontsize=12)
    ax.set_title(f"{latex_title} - Real Susceptibility", fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    
    plt.tight_layout()
    return fig


def plot_single_chi_imaginary(chi_data: pd.DataFrame, composition: str,
                              sample_id: str = None,
                              window_length: int = 31, polyorder: int = 3,
                              figsize=(8, 6)) -> plt.Figure:
    """
    Plot imaginary part of susceptibility for a single sample.
    
    Args:
        chi_data: DataFrame with chi data for one sample
        composition: Chemical composition string
        sample_id: Optional sample_id for trimming logic
        window_length: Savitzky-Golay filter window
        polyorder: Savitzky-Golay polynomial order
        figsize: Figure size
        
    Returns:
        matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if chi_data is None or len(chi_data) == 0:
        ax.text(0.5, 0.5, 'No susceptibility data available',
               ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig
    
    latex_title = composition_to_latex(composition)
    needs_trimming = sample_id in ['0455_HM_NbTa2Zr', '0457_HM_MoNbZr2']
    
    for field in sorted(chi_data['field_tesla'].unique(), reverse=True):
        field_data = chi_data[chi_data['field_tesla'] == field].sort_values('temperature')
        
        if len(field_data) == 0:
            continue
        
        if needs_trimming and abs(field - 2.0) < 0.1:
            field_data = trim_high_temp_outliers(field_data, n_trim=3)
        
        if len(field_data) < window_length:
            continue
        
        high_temp_value = field_data.iloc[-1]['chan2_v']
        normalized = (field_data['chan2_v'] - high_temp_value) * 1e6
        smoothed = savgol_filter(normalized,
                                window_length=window_length,
                                polyorder=polyorder)
        
        ax.plot(field_data['temperature'], smoothed,
               label=f'{field:.1f} T', alpha=0.7, linewidth=2)
    
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel(r"$\chi''$ offset ($\mu$V)", fontsize=12)
    ax.set_title(f"{latex_title} - Imaginary Susceptibility", fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    
    plt.tight_layout()
    return fig


def plot_single_hc2(tc_data: pd.DataFrame, composition: str,
                    figsize=(8, 6)) -> tuple:
    """
    Plot Hc2(T) with fits for a single sample.
    
    Args:
        tc_data: DataFrame with Tc vs field data
        composition: Chemical composition string
        figsize: Figure size
        
    Returns:
        tuple of (figure, fit_results_dict or None)
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if tc_data is None or len(tc_data) == 0:
        ax.text(0.5, 0.5, 'No Hc2 data available',
               ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig, None
    
    latex_title = composition_to_latex(composition)
    fit_results = fit_hc2_models(tc_data)
    
    if fit_results is None:
        ax.text(0.5, 0.5, 'Insufficient data for fitting',
               ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig, None
    
    # Determine plot ranges
    max_tc = np.ceil(tc_data['tc_kelvin'].max())
    max_field_data = tc_data['field_tesla'].max()
    max_field_fit = max(fit_results['linear']['Hc2_0'],
                       fit_results['quadratic']['Hc2_0'])
    max_field = np.ceil(max(max_field_data, max_field_fit))
    
    # Plot data
    ax.plot(tc_data['tc_kelvin'], tc_data['field_tesla'],
           'o', markersize=10, color='black', label='Data', zorder=3)
    
    # Plot fits
    T_fit_linear = np.linspace(0, fit_results['linear']['Tc'], 200)
    H_linear = linear_hc2_model(T_fit_linear,
                                fit_results['linear']['Hc2_0'],
                                fit_results['linear']['Tc'])
    H_linear = np.clip(H_linear, 0, None)
    ax.plot(T_fit_linear, H_linear, '-', linewidth=2.5,
           color='C0', label='Linear', alpha=0.7)
    
    T_fit_quadratic = np.linspace(0, fit_results['quadratic']['Tc'], 200)
    H_quadratic = quadratic_hc2_model(T_fit_quadratic,
                                     fit_results['quadratic']['Hc2_0'],
                                     fit_results['quadratic']['Tc'])
    H_quadratic = np.clip(H_quadratic, 0, None)
    ax.plot(T_fit_quadratic, H_quadratic, '--', linewidth=2.5,
           color='C1', label='Quadratic', alpha=0.7)
    
    ax.set_xlim(0, max_tc)
    ax.set_ylim(0, max_field)
    ax.set_xlabel('$T_c$ (K)', fontsize=12)
    ax.set_ylabel('Field (T)', fontsize=12)
    ax.set_title(f"{latex_title} - Upper Critical Field", fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=11)
    
    plt.tight_layout()
    return fig, fit_results