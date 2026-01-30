#!/usr/bin/env python3
"""
Quick dataframe inspection tool.
"""

import pandas as pd
from pathlib import Path
import sys

def main():
    pkl_path = Path("data/processed/synthesis_data.pkl")
    
    if not pkl_path.exists():
        print(f"Error: {pkl_path} not found")
        return
    
    df = pd.read_pickle(pkl_path)
    
    # Get format from command line (default: markdown)
    format_type = sys.argv[1] if len(sys.argv) > 1 else "markdown"
    
    # Drop large columns for cleaner view
    df_view = df.drop(columns=['files', 'status_content', 'synthesis_content', 'xrd_patterns'], errors='ignore')
    
    if format_type == "markdown":
        output = df_view.to_markdown(index=False)
        print(output)
    elif format_type == "html":
        output_path = Path("data/temp/dataframe_view.html")
        output_path.parent.mkdir(exist_ok=True)
        df_view.to_html(output_path, index=False)
        print(f"✓ Saved to {output_path}")
    elif format_type == "excel":
        output_path = Path("data/temp/dataframe_view.xlsx")
        output_path.parent.mkdir(exist_ok=True)
        df_view.to_excel(output_path, index=False)
        print(f"✓ Saved to {output_path}")
    elif format_type == "info":
        print(df.info())
        print("\n" + "="*80)
        print("First few rows:")
        print(df_view.head())
    else:
        print("Usage: python scripts/inspect_dataframe.py [markdown|html|excel|info]")

if __name__ == "__main__":
    main()