import re
import pandas as pd

# --- 1) Optional: bias lookup (fallback when no explicit order is provided) ---
bias_map = {}
try:
    index_mapping_new = pd.read_csv('all_methods_statistics_251110True.csv')  # fixed missing quote
    bias_map = dict(zip(index_mapping_new["method_name"], index_mapping_new["bias50"]))
except Exception as e:
    print("not found newest DATAFRAME SAVED", e)

def bias50(method_name: str) -> float:
    # Fallback bias used only when no explicit order applies
    return bias_map.get(method_name, float('inf'))

# --- 2) Explicit within-class orders you requested ---
ORDER_FAIR = [
    "FaIR_anthro","FaIR_anthro_unB", "FaIR_nonat",  "FaIR_nonat_unB", "FaIR_all",
     "FaIR_all_unB",
    "FaIR_comb_unB"
]


ORDER_CGWL = [
    "CGWL10y_forec", "CGWL10y_for_halfU",
    "CGWL10y_pUKCP", "CGWL10y_sUKCP",
    "CGWL10y_sfUKCP", "CGWL10y_IPCC"
]

# Helper: find the first token from an order list that occurs in the method name
def index_in_named_list(method_name: str, tokens: list[str]) -> tuple[int, int]:
    for idx, tok in enumerate(tokens):
        if tok in method_name:
            return (0, idx)  # (matched, position)
    return (1, len(tokens))  # didn't match any token

# --- 3) Special handling for the LOWESS family inside ST_Fits class ---
LOWESS_FAMILY_PREFIX = "lowess"  # treat anything containing 'lowess' as same family

def lowess_family_key(method_name: str) -> tuple[int, str]:
    """
    Group all LOWESS variants together and sort them consistently.
    If it's LOWESS -> (0, normalized_name)
    else           -> (1, original_name)
    """
    if LOWESS_FAMILY_PREFIX in method_name.lower():
        # Normalize slightly so 'lowess1dg20wnc' and close variants cluster cleanly.
        # (Light touch—feel free to refine this to your exact scheme.)
        norm = re.sub(r'[^a-z0-9_]+', '', method_name.lower())
        return (0, norm)
    return (1, method_name.lower())


ORDER_HUMAN_INDUCED = [
    "eROF_anthro", "KCC_human", "eROF_tot",  "GWI_anthro", "GWI_anthro_orig", "GWI_anthro_CGWL", "GWI_anthro_SR15",
    "GWI_anthro_AR6", "eROF_tot", "KCC_all",
     "GWI_tot","GWI_tot_orig",
    "GWI_tot_CGWL","GWI_tot_SR15", "GWI_tot_AR6", 
]

HUMAN_INDUCED_RANK = {name: i for i, name in enumerate(ORDER_HUMAN_INDUCED)}
FAIR_INDUCED_RANK = {name: i for i, name in enumerate(ORDER_FAIR)}

def order_within_class_key(method_name: str, method_class: str):
    if method_class == "43_Forcing_Based/3_Human_Induced":
        if method_name in HUMAN_INDUCED_RANK:                 # exact match only
            return (0, HUMAN_INDUCED_RANK[method_name])
        # not in the explicit list: push after the explicit block
        return (1, bias50(method_name), method_name.lower())

    # (unchanged) FaIR block uses exact matching too:
    if method_class == "43_Forcing_Based/1_ERF_FaIR":
        if method_name in FAIR_INDUCED_RANK:                 # exact match only
            return (0, FAIR_INDUCED_RANK[method_name])
        elif not hasattr(order_within_class_key, "_fair_rank"):
            fair_list = [
                "FaIR_anthro", "FaIR_nonat", "FaIR_all",
                "FaIR_anthro_unB", "FaIR_nonat_unB", "FaIR_all_unB",
                "FaIR_comb_unB"
            ]
            order_within_class_key._fair_rank = {n: i for i, n in enumerate(fair_list)}
        fair_rank = order_within_class_key._fair_rank
        if method_name in fair_rank:
            return (0, fair_rank[method_name])
        return (1, bias50(method_name), method_name.lower())

    # (unchanged) CGWL exact list
    if method_class == "44_EarthModel_CGWL":
        cgwl_list = [
            "CGWL10y_forec", "CGWL10y_for_halfU",
            "CGWL10y_pUKCP", "CGWL10y_sUKCP",
            "CGWL10y_sfUKCP", "CGWL10y_IPCC"
        ]
        cgwl_rank = {n: i for i, n in enumerate(cgwl_list)}
        if method_name in cgwl_rank:
            return (0, cgwl_rank[method_name])
        return (1, bias50(method_name), method_name.lower())

    # (unchanged) LOWESS grouping for ST_Fits
    if method_class == "42_Temp_Alone/3_ST_Fits":
        fam_group, fam_name = lowess_family_key(method_name)
        return (fam_group, fam_name) if fam_group == 0 else (1, bias50(method_name), method_name.lower())

    # default
    return (bias50(method_name), method_name.lower())

# --- 5) Final sort: primary by method_class, secondary by class-aware order key ---
# results: dict[str, dict], where each value has at least 'method_class'
# Example:
# sorted_results = sorted(results.items(), key=lambda item: (item[1]['method_class'], order_within_class_key(item[0], item[1]['method_class'])))

# If you previously expected a list of (name, payload) pairs:
def sort_results(results: dict[str, dict]):
    return sorted(
        results.items(),
        key=lambda item: (
            item[1]['method_class'],
            order_within_class_key(item[0], item[1]['method_class'])
        )
    )
