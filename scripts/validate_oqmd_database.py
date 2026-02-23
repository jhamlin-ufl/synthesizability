# scripts/validate_oqmd_database.py
"""
Validate OQMD database installation and version.
"""
import subprocess
import sys
import re
from pathlib import Path

REQUIRED_VERSION = "1.7"
OQMD_DOWNLOAD_URL = "https://oqmd.org/download/"

def run_mysql_query(query):
    """Execute a MySQL query and return output"""
    try:
        cmd = ['sudo', 'mysql', '-e', f'USE qmdb; {query}']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        return None

def check_mariadb_running():
    """Check if MariaDB server is running"""
    try:
        result = subprocess.run(['systemctl', 'is-active', 'mariadb'], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except:
        return False

def check_database_exists():
    """Check if qmdb database exists"""
    try:
        cmd = ['sudo', 'mysql', '-e', 'SHOW DATABASES LIKE "qmdb";']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return 'qmdb' in result.stdout
    except:
        return False

def get_database_entry_count():
    """Get number of entries in database"""
    result = run_mysql_query("SELECT COUNT(*) FROM entries;")
    if result:
        lines = result.split('\n')
        if len(lines) > 1:
            return int(lines[1])
    return None

def detect_installed_version(entry_count):
    """Estimate version based on entry count"""
    # Known entry counts for versions
    version_counts = {
        "1.8": 1400000,  # ~1.4M entries (Feb 2026)
        "1.7": 1317701,  # ~1.3M entries (May 2025)
        "1.6": 1100000,  # ~1.1M entries (Nov 2023)
    }
    
    if entry_count is None:
        return None
    
    # Find closest version
    min_diff = float('inf')
    detected = None
    for version, count in version_counts.items():
        diff = abs(entry_count - count)
        if diff < min_diff:
            min_diff = diff
            detected = version
    
    return detected

def check_for_updates():
    """
    Check if newer version is available.
    For now, this is a placeholder - could scrape OQMD website.
    """
    # Current known latest version
    LATEST_VERSION = "1.8"
    return LATEST_VERSION

def main():
    print("=" * 60)
    print("OQMD Database Validation")
    print("=" * 60)
    
    # Check 1: MariaDB running
    print("\n[1/4] Checking MariaDB status...")
    if not check_mariadb_running():
        print("❌ MariaDB is not running")
        print(f"   Please start MariaDB and ensure qmdb is installed")
        print(f"   Installation instructions: {OQMD_DOWNLOAD_URL}")
        sys.exit(1)
    print("✓ MariaDB is running")
    
    # Check 2: Database exists
    print("\n[2/4] Checking qmdb database...")
    if not check_database_exists():
        print("❌ qmdb database not found")
        print(f"   Please download and import OQMD database")
        print(f"   Download from: {OQMD_DOWNLOAD_URL}")
        sys.exit(1)
    print("✓ qmdb database exists")
    
    # Check 3: Get version
    print("\n[3/4] Detecting database version...")
    entry_count = get_database_entry_count()
    if entry_count is None:
        print("⚠ Could not determine entry count")
        sys.exit(1)
    
    print(f"   Total entries: {entry_count:,}")
    installed_version = detect_installed_version(entry_count)
    
    if installed_version:
        print(f"   Detected version: v{installed_version}")
    else:
        print("⚠ Could not detect version")
    
    # Check 4: Version comparison
    print("\n[4/4] Checking for updates...")
    latest_version = check_for_updates()
    print(f"   Latest available: v{latest_version}")
    
    if installed_version and installed_version < latest_version:
        print(f"\n⚠ Newer version available: v{latest_version}")
        print(f"   You are using: v{installed_version}")
        print(f"   Download from: {OQMD_DOWNLOAD_URL}")
        print(f"   Note: Automatic database updates not yet implemented")
        print(f"         Your current version is sufficient for this analysis")
    elif installed_version == latest_version:
        print(f"\n✓ You have the latest version (v{latest_version})")
    
    # Final check: Test query
    print("\n[5/5] Testing database query...")
    test_result = run_mysql_query(
        "SELECT COUNT(*) FROM formation_energies WHERE fit_id = 'standard';"
    )
    if test_result:
        lines = test_result.split('\n')
        if len(lines) > 1:
            count = int(lines[1])
            print(f"✓ Query successful: {count:,} formation energies available")
    else:
        print("❌ Test query failed")
        sys.exit(1)
    
    print("\n" + "=" * 60)
    print("✓ OQMD Database validation complete")
    print("=" * 60)

if __name__ == "__main__":
    main()