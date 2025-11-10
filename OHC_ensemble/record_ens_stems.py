import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Dict
import pdb

# ============================================================================
# CONFIGURATION
# ============================================================================

RECORDS = [
    'IAP', 'NCEILevitus', 'JMA', 'CORA', 'Miniere', 'GCOS', 
    'GREPrean', 'CISRODomingues', 'MOHeaCAN', 'EN4']

EN4_SUBTYPES = ['c14', 'c14', 'c14', 'c13', 'c13', 'c13', 'g10', 'g10', 'l09', 'l09']

YEARS = np.arange(1850, 2025)  # 175 years
N_ENSEMBLE = 100
N_COPIES_PER_RECORD = 10

# Starting year for forward sweep (maximum alignment)
START_YEAR = 2010

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_domain(domain_str: str) -> Tuple[int, int, str]:
    """Parse domain string like '0-2000m_60Lat' into (depth_min, depth_max, lat_coverage)"""
    if ':' in domain_str:
        # Handle already-combined domains, extract leftmost
        domain_str = domain_str.split(':')[0].strip()
    
    parts = domain_str.split('_')
    depth_part = parts[0].replace('m', '')
    depth_min, depth_max = map(int, depth_part.split('-'))
    lat_coverage = parts[1] if len(parts) > 1 else 'allLat'
    
    return depth_min, depth_max, lat_coverage

rt="observ_late"

def load_record_data(record_name: str, subtype: str = None) -> pd.DataFrame:
    """Load main OHCA data for a record"""
    if subtype:
        filename = f"{rt}/{record_name}/{subtype}/ohca_{subtype}.csv"
    else:
        filename = f"{rt}/{record_name}/ohca_{record_name}.csv"  
    
    df = pd.read_csv(filename).apply(pd.to_numeric)
    df['year']=df['year'].apply(np.floor)
    return df


def load_supplemental_data(record_name: str, subtype: str = None) -> pd.DataFrame:
    """Load supplemental OHCA chunks for a record"""
    if subtype:
        filename = f"{rt}/{record_name}/{subtype}/ohca_supp_{subtype}.csv"
    else:
        filename = f"{rt}/{record_name}/ohca_supp_{record_name}.csv"
    
    try:
        df = pd.read_csv(filename)
       # print(filename)
        df=df.apply(pd.to_numeric)
        df['year']=df['year'].apply(np.floor)
        return df
    except FileNotFoundError:
        return None



def compute_annual_differences(df: pd.DataFrame, value_cols: List[str]) -> pd.DataFrame:
    """Compute year-to-year differences for OHCA columns"""
    df_diff = df.copy()
    
    for col in value_cols:
        if col in df.columns and not col.endswith('_se'):
            df_diff[col] = df[col].diff()
    
    return df_diff


def get_record_subtype_cycle() -> List[Tuple[str, str]]:
    """Generate the cycling pattern for records and subtypes"""
    cycle = []
    
    for record in RECORDS:
        if record == 'EN4':
            for subtype in EN4_SUBTYPES:
                cycle.append(('EN4', subtype))
        else:
            for _ in range(N_COPIES_PER_RECORD):
                cycle.append((record, None))
    
    return cycle


def combine_ohca_values(val1: float, val2: float) -> float:
    """Combine two OHCA values (simple addition)"""
    if np.isnan(val1):
        return val2
    if np.isnan(val2):
        return val1
    return val1 + val2


def encode_asymmetric_uncertainty(pos_error: float, neg_error: float) -> complex:
    """
    Encode asymmetric uncertainties as complex number.
    pos_error * (1 + 1j) + neg_error * (1 - 1j)
    
    Args:
        pos_error: Positive standard error
        neg_error: Negative standard error (should be positive value)
    
    Returns:
        Complex number encoding both errors
    """
    return pos_error * (1 + 1j) + neg_error * (1 - 1j)


def decode_asymmetric_uncertainty(encoded: complex) -> Tuple[float, float]:
    """
    Decode complex number back to asymmetric uncertainties.
    
    Args:
        encoded: Complex number from encode_asymmetric_uncertainty
    
    Returns:
        (pos_error, neg_error) tuple
    """
    # From encoded = pos*(1+1j) + neg*(1-1j) = (pos+neg) + 1j*(pos-neg)
    # We get: pos = (real + imag)/2, neg = (real - imag)/2
    pos_error = (encoded.real + encoded.imag) / 2
    neg_error = (encoded.real - encoded.imag) / 2
    return pos_error, neg_error


def is_symmetric_uncertainty(encoded: complex) -> bool:
    """Check if encoded uncertainty is symmetric (imag part ~= 0)"""
    return abs(encoded.imag) < 1e-10

def combine_uncertainties(unc1: float, unc2: float) -> float:
    """Combine two uncertainties (root sum of squares)"""

    #something funny going on with unexpected inputs
    erra = False
    try:
        erra = erra or np.isnan(unc1)
    except:
        breakpoint()
    
    try:
        erra = erra or np.isnan(unc2)
    except:
        breakpoint()

    if erra:
        breakpoint()

        
    return np.sqrt(unc1**2 + unc2**2)

##def combine_uncertainties(unc1: float, unc2: float) -> complex:
##    """
##    Combine two uncertainties (root sum of squares).
##    Handles both float (symmetric) and complex (asymmetric) uncertainties.
##    
##    Args:
##        unc1, unc2: Can be float or complex
##    
##    Returns:
##        Combined uncertainty (float if both symmetric, complex otherwise)
##    """
##    # Handle NaN cases
##    if np.isnan(unc1) or (isinstance(unc1, complex) and np.isnan(unc1.real)):
##        return unc2
##    if np.isnan(unc2) or (isinstance(unc2, complex) and np.isnan(unc2.real)):
##        return unc1
##    
##    # Convert to complex if needed
##    u1 = complex(unc1) if not isinstance(unc1, complex) else unc1
##    u2 = complex(unc2) if not isinstance(unc2, complex) else unc2
##    
##    # For symmetric uncertainties: combine as usual
##    if is_symmetric_uncertainty(u1) and is_symmetric_uncertainty(u2):
##        combined = np.sqrt(u1.real**2 + u2.real**2)
##        return combined
##    
##    # For asymmetric: decode, combine pos/neg separately, re-encode
##    pos1, neg1 = decode_asymmetric_uncertainty(u1)
##    pos2, neg2 = decode_asymmetric_uncertainty(u2)
##    
##    combined_pos = np.sqrt(pos1**2 + pos2**2)
##    combined_neg = np.sqrt(neg1**2 + neg2**2)
##    
##    return encode_asymmetric_uncertainty(combined_pos, combined_neg)


def combine_domain_strings(base_domain: str, base_string: str, 
                           infill_domain: str, infill_record: str) -> str:
    """Combine domain strings with proper notation"""
    return f"{base_domain} : ({base_string} + {infill_domain}_{infill_record})"


# ============================================================================
# INITIALIZATION
# ============================================================================

def initialize_arrays(record_data: Dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Initialize the three main arrays:
    - ohca_change: Annual ΔOHCA (differenced)
    - ohca_uncertainty: Standard errors (NOT differenced)
    - coverage_string: Domain notation strings
    """
    ohca_change = np.full((N_ENSEMBLE, len(YEARS)), np.nan)
    ohca_uncertainty = np.full((N_ENSEMBLE, len(YEARS)), np.nan, dtype=complex)
    coverage_string = np.empty((N_ENSEMBLE, len(YEARS)), dtype=object)
    
    # Initialize with records (10 copies each, grouped by rows)
    row_idx = 0
    
    for record, subtype in get_record_subtype_cycle():
        # Load the record's main data
        df_main = record_data[record]['main'][subtype] if subtype else record_data[record]['main']
        
        # Get primary domain column (first non-year, non-_se column)
        value_cols = [c for c in df_main.columns if c != 'year' and not c.endswith('_se') and not c.endswith('_pse') and not c.endswith('_nse')]
        primary_domain = value_cols[0] if value_cols else None
        if primary_domain:
            for year_idx, year in enumerate(YEARS):
                year_data = df_main[df_main['year'] == year]
                
                if not year_data.empty:
                    if len(value_cols)>1:
                        for i in range(len(value_cols)):
                            ei = len(value_cols)-1 -i
                            #go through them from right to left
                            if not np.isnan(year_data[value_cols[ei]].values[0]):
                                primary_domain = value_cols[ei] #reset primary_domain to the right-most one with valid change
                                break
                    # OHCA change (already differenced)
                    ohca_change[row_idx, year_idx] = year_data[primary_domain].values[0]
                    
                    # Uncertainty (NOT differenced)
                    # Uncertainty - check if asymmetric (MOHeaCAN case)
                    unc_col_pos = f"{primary_domain}_pse"
                    unc_col_neg = f"{primary_domain}_nse"
                    unc_col_sym = f"{primary_domain}_se"
                    
                    if unc_col_pos in df_main.columns and unc_col_neg in df_main.columns:
                        # Asymmetric uncertainty (MOHeaCAN)
                        pos_err = year_data[unc_col_pos].values[0]
                        neg_err = year_data[unc_col_neg].values[0]
                        ohca_uncertainty[row_idx, year_idx] = encode_asymmetric_uncertainty(pos_err, neg_err)
                    elif unc_col_sym in df_main.columns:
                        # Symmetric uncertainty (most records)
                        sym_err = year_data[unc_col_sym].values[0]
                        ohca_uncertainty[row_idx, year_idx] = complex(sym_err, 0)

                    # Coverage string: domain_recordname
                    record_str = f"{record}.{subtype}" if subtype else record
                    coverage_string[row_idx, year_idx] = f"{primary_domain}_{record_str}"
        
        row_idx += 1
    
    return ohca_change, ohca_uncertainty, coverage_string


def load_all_record_data() -> Dict:
    """Load all record data (main and supplemental) into memory"""
    record_data = {}
    
    for record in RECORDS:
        record_data[record] = {'main': {}, 'supp': {}}
        
        if record == 'EN4':
            # Load all EN4 subtypes
            for subtype in set(EN4_SUBTYPES):
                df_main = load_record_data('EN4', subtype)
                # Compute differences for OHCA columns
                value_cols = [c for c in df_main.columns if c != 'year' and not c.endswith('_se')]
                df_main = compute_annual_differences(df_main, value_cols)
                record_data[record]['main'][subtype] = df_main
                
                df_supp = load_supplemental_data('EN4', subtype)
                if df_supp is not None:
                    df_supp = compute_annual_differences(df_supp, 
                        [c for c in df_supp.columns if c != 'year' and not c.endswith('_se')])
                record_data[record]['supp'][subtype] = df_supp
        else:
            # Load regular record
            df_main = load_record_data(record)
            value_cols = [c for c in df_main.columns if c != 'year' and not c.endswith('_se')]
            df_main = compute_annual_differences(df_main, value_cols)
            record_data[record]['main'] = df_main
            
            df_supp = load_supplemental_data(record)
            if df_supp is not None:
                df_supp = compute_annual_differences(df_supp,
                    [c for c in df_supp.columns if c != 'year' and not c.endswith('_se')])
            record_data[record]['supp'] = df_supp
    
    # Load Zanna (for backsweep only, pre-padded to 1850)
    df_zanna = load_record_data('zanna')
    value_cols = [c for c in df_zanna.columns if c != 'year' and not c.endswith('_se')]
    df_zanna = compute_annual_differences(df_zanna, value_cols)

    df_zanna_supp = load_supplemental_data('zanna')
    value_cols_supp = [c for c in df_zanna_supp.columns if c != 'year' and not c.endswith('_se')]
    df_zanna_supp = compute_annual_differences(df_zanna_supp, value_cols_supp)
    
    record_data['Zanna'] = {'main': df_zanna, 'supp': df_zanna_supp}
    
    return record_data


# ============================================================================
# INFILLING FUNCTIONS
# ============================================================================

def get_infill_cycle(available_records: List[str], 
                     exclude_records: List[str] = None) -> List[Tuple[str, str]]:
    """
    Generate cycling pattern from available records.
    First cycles through unique base records, then replaces EN4 with its subtypes.
    
    Args:
        available_records: List of (record_name, subtype) tuples
        exclude_records: List of record names to exclude from cycle
    
    Returns:
        List of (record_name, subtype) tuples in cycle order
    """
    if exclude_records is None:
        exclude_records = []
    
    # Get unique base records (collapse EN4 subtypes to just 'EN4')
    unique_base_records = []
    seen = set()
    
    for record in available_records:
        if record in exclude_records:
            continue
            
        if record not in seen: #removes potential duplicates
            unique_base_records.append(record)
            seen.add(record)
    
    # Build cycle: alternate through base records
    cycle = []
    if 'EN4' in unique_base_records:
        for en4_subtype in ['c14', 'c13', 'g10', 'l09']:
            for base_record in unique_base_records:
                if base_record == 'EN4':
                    # Replace EN4 with its subtypes
                    cycle.append(('EN4', en4_subtype))
                else:
                    cycle.append((base_record, None))
    else:
        for base_record in unique_base_records:
            cycle.append((base_record, None))
    
    return cycle


def get_available_records_for_year(year: int, record_data: Dict, 
                                   include_zanna: bool = False) -> List[Tuple[str, str]]:
    """
    Get list of all available records for a given year.
    
    Returns:
        List of (record_name, subtype) tuples that have data for this year
    """
    available = []
    
    for record in RECORDS:
        if record == 'EN4':
            for subtype in set(EN4_SUBTYPES):
                df = record_data[record]['main'][subtype]
                if not df[df['year'] == year].empty:
                    available.append(('EN4', subtype))
        else:
            df = record_data[record]['main']
            if not df[df['year'] == year].empty:
                available.append((record, None))
    
    if include_zanna:
        df_zanna = record_data['Zanna']['main']
        if not df_zanna[df_zanna['year'] == year].empty:
            available.append(('Zanna', None))
    
    return available


def infill_latitude(ohca_change: np.ndarray, ohca_uncertainty: np.ndarray,
                    coverage_string: np.ndarray, year_idx: int,
                    record_data: Dict) -> None:
    """
    Step 1: Infill latitude coverage (60Lat → allLat)
    Applies to GCOS and GREPrean
    Uses cycling through available records (CORA and EN4 subtypes)
    """
    year = YEARS[year_idx]
    
    # Get available records for latitude infill (CORA and EN4 have supplemental data)
    available_for_lat_infill = ['CORA','EN4']
    
    # Generate cycle: CORA, EN4, CORA, EN4, ... then EN4 → c14, c13, g10, l09
    lat_infill_cycle = get_infill_cycle(available_for_lat_infill)

    cycle_idx = 0
    
    for row_idx in range(N_ENSEMBLE):
        domain_str = coverage_string[row_idx, year_idx]
        
        if domain_str is None or pd.isna(domain_str):
            continue
        
        # Parse current domain
        depth_min, depth_max, lat_cov = parse_domain(domain_str)
        
        # Check if latitude infill needed
        if lat_cov == '60Lat':
            # Need to add polar latitudes
            infill_record, infill_subtype = lat_infill_cycle[cycle_idx % len(lat_infill_cycle)]
            cycle_idx += 1
            
            # Get supplemental data
            if infill_subtype:
                df_supp = record_data[infill_record]['supp'][infill_subtype]
            else:
                df_supp = record_data[infill_record]['supp']
            
            if df_supp is not None:
                # Look for polar complement
                polar_domain = f"{depth_min}-{depth_max}m_polarLat"
                
                if polar_domain in df_supp.columns:
                    year_data = df_supp[df_supp['year'] == year]
                    
                    if not year_data.empty:
                        polar_val = year_data[polar_domain].values[0]
                        polar_unc_col = f"{polar_domain}_se"
                        polar_unc = year_data[polar_unc_col].values[0] if polar_unc_col in df_supp.columns else 0
                        
                        # Combine values
                        ohca_change[row_idx, year_idx] = combine_ohca_values(
                            ohca_change[row_idx, year_idx], polar_val)
                        ohca_uncertainty[row_idx, year_idx] = combine_uncertainties(
                            ohca_uncertainty[row_idx, year_idx], polar_unc)
                        
                        # Update coverage string
                        infill_str = f"{infill_record}.{infill_subtype}" if infill_subtype else infill_record
                        new_domain = f"{depth_min}-{depth_max}m_allLat"
                        coverage_string[row_idx, year_idx] = combine_domain_strings(
                            new_domain, domain_str, polar_domain, infill_str)


def infill_depth(ohca_change: np.ndarray, ohca_uncertainty: np.ndarray,
                coverage_string: np.ndarray, year_idx: int,
                record_data: Dict, include_zanna: bool = False) -> None:
    """
    Step 2: Infill depth coverage (partial → 0-6000m)
    Searches for complementary depth chunks in supplemental files
    If include_zanna=True (backsweep), adds Zanna to available records
    """
    year = YEARS[year_idx]
    for row_idx in range(N_ENSEMBLE):
        domain_str = coverage_string[row_idx, year_idx]
        
        if domain_str is None or pd.isna(domain_str):
            continue
        
        # Parse current domain
        depth_min, depth_max, lat_cov = parse_domain(domain_str)
        
        # Check if depth infill needed (not already 0-6000m)
        if depth_max < 6000:
            # Need deeper data
            target_depth_min = depth_max
            target_depth_max = 6000
            target_domain = f"{target_depth_min}-{target_depth_max}m_allLat"
            
            # Search all records for this depth chunk
            available_chunks = []
            
            # Search regular records
            for record in RECORDS:
                if record == 'EN4':
                    for subtype in ['c14', 'c13', 'g10', 'l09']:
                        df_supp = record_data[record]['supp'][subtype]
                        if df_supp is not None and target_domain in df_supp.columns:
                            year_data = df_supp[df_supp['year'] == year]
                            if not year_data.empty and not pd.isna(year_data[target_domain].values[0]):
                                available_chunks.append((record, subtype, df_supp, year_data))
                else:
                    df_supp = record_data[record]['supp']
                    if df_supp is not None and target_domain in df_supp.columns:
                        year_data = df_supp[df_supp['year'] == year]
                        if not year_data.empty and not pd.isna(year_data[target_domain].values[0]):
                            available_chunks.append((record, None, df_supp, year_data))
            
            # Add Zanna if in backsweep
            if include_zanna:
                df_zanna = record_data['Zanna']['supp']
                if target_domain in df_zanna.columns:
                    year_data = df_zanna[df_zanna['year'] == year]
                    if not year_data.empty and not pd.isna(year_data[target_domain].values[0]):
                        available_chunks.append(('Zanna', None, df_zanna, year_data))
            
            if available_chunks:
                # Convert to (record, subtype) tuples for cycling
                chunk_records = [rec for rec, _, _, _ in available_chunks]
                depth_cycle = get_infill_cycle(chunk_records)
                # Pick one chunk using consistent ordering
                chosen_idx = row_idx % len(depth_cycle)
                chosen_record, chosen_subtype = depth_cycle[chosen_idx]
                
                # Find the matching chunk data
                for rec, sub, df, data in available_chunks:
                    if rec == chosen_record and sub == chosen_subtype:
                        chosen_df = df
                        chosen_data = data
                        break
                
                depth_val = chosen_data[target_domain].values[0]
                depth_unc_col = f"{target_domain}_se" #pse and nse only apply to full-depth MOHeaCAN
                depth_unc = chosen_data[depth_unc_col].values[0] if depth_unc_col in chosen_df.columns else 0
                
                # Combine values
                ohca_change[row_idx, year_idx] = combine_ohca_values(
                    ohca_change[row_idx, year_idx], depth_val)
                ohca_uncertainty[row_idx, year_idx] = combine_uncertainties(
                    ohca_uncertainty[row_idx, year_idx], depth_unc)
                
                # Update coverage string
                infill_str = f"{chosen_record}.{chosen_subtype}" if chosen_subtype else chosen_record
                new_domain = f"{depth_min}-6000m_allLat"
                coverage_string[row_idx, year_idx] = combine_domain_strings(
                    new_domain, domain_str, target_domain, infill_str)


def infill_missing_years(ohca_change: np.ndarray, ohca_uncertainty: np.ndarray,
                        coverage_string: np.ndarray, year_idx: int, record_zanna: pd.DataFrame,
                        include_zanna: bool = False) -> None:
    """
    Step 3: Infill completely empty cells using completed records from same year
    Alternates between all complete records (0-6000m_allLat) available that year
    """
    # Find all complete records in this year-column
    complete_records = []
    
    for row_idx in range(N_ENSEMBLE):
        domain_str = coverage_string[row_idx, year_idx]
        
        if domain_str is not None and not pd.isna(domain_str):
            depth_min, depth_max, lat_cov = parse_domain(domain_str)
            
            if depth_min == 0 and depth_max == 6000 and lat_cov == 'allLat':
                complete_records.append((row_idx, domain_str))

    if include_zanna:
        complete_records.append((-1, "0-6000m_allLat_Zanna"))
    
    if not complete_records:
        return
    
    # Extract unique record combinations from complete records (likely containing duplicates)
    unique_complete = []
    unique_complete_withEN4 = []
    for row_idx, domain_str in complete_records:
        if "EN4." in domain_str:
            if domain_str not in [d for _, d in unique_complete_withEN4]:
                unique_complete_withEN4.append((row_idx, domain_str))
        else:
            if domain_str not in [d for _, d in unique_complete]:
                unique_complete.append((row_idx, domain_str))
    
    # (This naturally handles EN4 subtypes being different)
    # Generate supercycle: interleave EN4 variants with non-EN4 records
    temporal_cycle = []

    max_len = max(len(unique_complete_withEN4), 1)  # Avoid division by zero

    for i in range(max_len):
        if i < len(unique_complete_withEN4):
            temporal_cycle.append(unique_complete_withEN4[i][0])
        
        # Add all non-EN4 records after each EN4 variant
        for row_idx, _ in unique_complete:
            temporal_cycle.append(row_idx)

    # Handle case where no EN4 variants exist
    if not temporal_cycle and unique_complete:
        temporal_cycle = [row_idx for row_idx, _ in unique_complete]
    
    
    # Now fill empty cells by cycling through complete records
    cycle_idx = 0
    
    for row_idx in range(N_ENSEMBLE):
        domain_str = coverage_string[row_idx, year_idx]
        
        # Check if cell is empty
        if domain_str is None or pd.isna(domain_str) or np.isnan(ohca_change[row_idx, year_idx]):
            # Copy from complete record
            source_row = temporal_cycle[cycle_idx % len(temporal_cycle)]
            cycle_idx += 1
            
            if source_row == -1: #here completely relying on Zanna (copying from hidden row), also has 1850 starting row

                ohca_change[row_idx, year_idx] = record_zanna.loc[year_idx,"0-6000m_allLat"]
                ohca_uncertainty[row_idx, year_idx] = record_zanna.loc[year_idx,"0-6000m_allLat_se"]
                coverage_string[row_idx, year_idx] = "0-6000m_allLat_Zanna"
            
            else:
                ohca_change[row_idx, year_idx] = ohca_change[source_row, year_idx]
                ohca_uncertainty[row_idx, year_idx] = ohca_uncertainty[source_row, year_idx]
                coverage_string[row_idx, year_idx] = coverage_string[source_row, year_idx]


# ============================================================================
# MAIN ALGORITHM
# ============================================================================

def process_year_column(ohca_change: np.ndarray, ohca_uncertainty: np.ndarray,
                       coverage_string: np.ndarray, year_idx: int,
                       record_data: Dict, include_zanna: bool = False) -> None:
    """Process one year-column: latitude infill → depth infill → temporal infill"""
    
    # Step 1: Infill latitude
    infill_latitude(ohca_change, ohca_uncertainty, coverage_string, year_idx, record_data)
    
    # Step 2: Infill depth
    infill_depth(ohca_change, ohca_uncertainty, coverage_string, year_idx, 
                record_data, include_zanna)
    
    # Step 3: Infill missing years
    infill_missing_years(ohca_change, ohca_uncertainty, coverage_string, year_idx, 
                        record_data['Zanna']['main'],include_zanna )


def create_ensemble():
    """Main function to create 100-member ensemble"""
    
    print("Loading all record data...")
    record_data = load_all_record_data()
    #fix missing Zanna initial point
    record_data['Zanna']['main'].loc[0,"0-6000m_allLat"]=0
    
    print("Initializing arrays...")
    ohca_change, ohca_uncertainty, coverage_string = initialize_arrays(record_data)

    
        
    
    # Find index of starting year
    start_idx = np.where(YEARS == START_YEAR)[0][0]
    
    print(f"Processing year {START_YEAR} (maximum alignment)...")
    process_year_column(ohca_change, ohca_uncertainty, coverage_string, 
                       start_idx, record_data, include_zanna=False)
    
    # Forward sweep: 2011 → 2024
    print("Forward sweep...")
    for year_idx in range(start_idx + 1, len(YEARS)):
        year = YEARS[year_idx]
        print(f"  Processing year {year}...")
        process_year_column(ohca_change, ohca_uncertainty, coverage_string,
                          year_idx, record_data, include_zanna=False)

    
    
    # Backward sweep: 2009 → 1850 (includes Zanna)
    print("Backward sweep (including Zanna)...")
    for year_idx in range(start_idx - 1, -1, -1):
        year = YEARS[year_idx]
        print(f"  Processing year {year}...")
        process_year_column(ohca_change, ohca_uncertainty, coverage_string,
                          year_idx, record_data, include_zanna=True)
    
    print("Ensemble creation complete!")


    
    return ohca_change, ohca_uncertainty, coverage_string, YEARS


# ============================================================================
# VISUALIZATION
# print(ohca_change[:,150:])
# ============================================================================
import sys
np.set_printoptions(threshold=sys.maxsize, precision=2, suppress=True,linewidth=np.inf)

def plot_ensemble_stems(ohca: np.ndarray, coverage_string: np.ndarray,
                        years: np.ndarray, title:str, filename : str, changes=True):
    """Plot the ensemble 'stems' showing coverage patterns"""
    """Plot the ensemble OHC timeseries with colored blocks and interactive datacursor"""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    # Find START_YEAR index
    start_idx = np.where(years == START_YEAR)[0][0]
    
    # Compute cumulative OHCA, anchored at 0 in START_YEAR
    if(changes):
        ohca_cumulative = np.nancumsum(ohca, axis=1)
        ohca_cumulative = ohca_cumulative - ohca_cumulative[:, start_idx:start_idx+1]
    else:
        ohca_cumulative = ohca
    
    # Setup colors: one per 10-row block
    n_blocks = N_ENSEMBLE // N_COPIES_PER_RECORD
    colors = cm.tab10(np.arange(n_blocks) % 10)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    lines = []
    for row_idx in range(N_ENSEMBLE):
        block_idx = row_idx // N_COPIES_PER_RECORD
        color = colors[block_idx]
        if row_idx % N_COPIES_PER_RECORD ==0:
            line, = ax.plot(years, ohca_cumulative[row_idx, :], 
                       color=color, alpha=0.4, linewidth=0.8,label=RECORDS[block_idx])
        else:
            line, = ax.plot(years, ohca_cumulative[row_idx, :], 
                       color=color, alpha=0.4, linewidth=0.8)
        lines.append(line)
    
    # Interactive datacursor
    def on_hover(event):
        if event.inaxes == ax:
            for row_idx, line in enumerate(lines):
                if line.contains(event)[0]:
                    year_idx = np.argmin(np.abs(years - event.xdata))
                    coverage = coverage_string[row_idx, year_idx]
                    ax.set_title(f"Row {row_idx}, Year {int(years[year_idx])}: {coverage}", 
                               fontsize=10)
                    fig.canvas.draw_idle()
                    return
    
    fig.canvas.mpl_connect('motion_notify_event', on_hover)

    ax.axvline(START_YEAR, color='red', linestyle='--', alpha=0.5, label='Anchor year')
    ax.set_xlabel('Year')
    ax.set_ylabel('Cumulative OHCA (J)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()

    plt.savefig(filename,  dpi=500)
    print(f"Saved {filename}")
    plt.show()


# ============================================================================
# RUN
# ============================================================================

if __name__ == "__main__":
    run_ensemble=False
    
    if(run_ensemble):
        ohca_change, ohca_uncertainty, coverage_string, years = create_ensemble()
    # Save outputs
        np.save('ohca_change_ensemble.npy', ohca_change)
        np.save('ohca_uncertainty_ensemble.npy', ohca_uncertainty)
        np.save('coverage_string_ensemble.npy', coverage_string)
    else:
        ohca_change = np.load('ohca_change_ensemble.npy')
        ohca_uncertainty = np.load('ohca_uncertainty_ensemble.npy')
        coverage_string =np.load('coverage_string_ensemble.npy', allow_pickle=True)
        years=YEARS
    
    # Plot
    plot_ensemble_stems(ohca_change, coverage_string, years, 'Ensemble OHCA Timeseries Stems',
                        "ensemble_stems")
    
    print("\nSaved arrays:")
    print("  - ohca_change_ensemble.npy (100 x 175)")
    print("  - ohca_uncertainty_ensemble.npy (100 x 175)")
    print("  - coverage_string_ensemble.npy (100 x 175)")
    print("  - years.npy (175,)")
