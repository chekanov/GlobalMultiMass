#!/usr/bin/env python3
import sys
import time
import argparse
import os
import json
import warnings
import numpy as np
from scipy.special import gammainc
from array import array
import ROOT

# Ensure pyBumpHunter is in path
sys.path.append("./pyBumpHunter/")
import pyBumpHunter as BH

# =================================================================
# 1. Math & ROOT Functions
# =================================================================
class FiveParam2015:
    def __call__(self, x, par):
        xx = x[0] / 13000.0
        if xx <= 0 or xx >= 1: return 0.0
        ff1 = par[0] * ROOT.TMath.Power((1.0 - xx), par[1])
        ff2 = ROOT.TMath.Power(xx, (par[2] + par[3] * ROOT.TMath.Log(xx) + par[4] * ROOT.TMath.Log(xx) * ROOT.TMath.Log(xx)))
        return ff1 * ff2

def FiveParam_NP(Ecm, x_center, p1, p2, p3, p4, p5):
    x = x_center / Ecm
    nlog = np.log(x)
    return p1 * np.power((1.0 - x), p2) * np.power(x, (p3 + p4 * nlog + p5 * nlog * nlog))

def fast_bumphunter_stat_with_loc(data_hist, bkg_hist, max_width=30):
    d_crop = np.asarray(data_hist).flatten()
    b_crop = np.asarray(bkg_hist).flatten()
    
    max_t, best_w, best_idx = 0.0, 0, 0
    for w in range(2, max_width + 1):
        k = np.ones(w)
        D = np.convolve(d_crop, k, mode='valid')
        B = np.convolve(b_crop, k, mode='valid')
        mask = (D > B) & (B > 0)
        if not np.any(mask): continue
        p = np.clip(gammainc(D[mask], B[mask]), 1e-300, 1.0)
        t_array = -np.log(p)
        if np.max(t_array) > max_t:
            max_t = np.max(t_array)
            best_w = w
            best_idx = np.where(mask)[0][np.argmax(t_array)]
    return max_t, best_idx, best_w

# =================================================================
# 2. Setup & Configuration
# =================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--events", type=int, default=100)
parser.add_argument("--zval", type=float, default=5.0)
parser.add_argument("-b", "--batch", action="store_true")
parser.add_argument("--fit", action="store_true")
parser.add_argument("--chimax", type=float, default=2.0)
parser.add_argument("--no-overlap", action="store_true", help="Disable JJ correlation overlap")
args = parser.parse_args()

if args.batch: ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal # Suppress GSL roundoff noise

log_file = open("run_log.txt", "w")
terminal = sys.__stdout__
sys.stdout = log_file

CHANNELS = ["jj", "bb", "jb", "je", "jm", "jg", "be", "bm", "bg"]
cms = 13000.0
MinEntries = 50

DEFALT_OVERLAP1 = {"jj": 0.0, "bb": 0.410, "jb": 0.635, "je": 0.325, "jm": 0.339, "jg": 0.062, "be": 0.204, "bm": 0.205, "bg": 0.028}
DEFALT_OVERLAP2 = {"jj": 0.0, "bb": 0.27, "jb": 0.37, "jm": 0.53, "je": 0.53, "jg": 0.01, "be": 0.24, "bm": 0.25, "bg": 0.01}
DEFALT_OVERLAP3 = {"jj": 0.0, "bb": 0.01, "jb": 0.28, "jm": 0.53, "je": 0.52, "jg": 0.01, "be": 0.16, "bm": 0.19, "bg": 0.01}
DEFALT_OVERLAP4 = {"jj": 0.0, "bb": 0.13, "jb": 0.14, "jm": 0.02, "je": 0.02, "jg": 0.99, "be": 0.01, "bm": 0.01, "bg": 0.24}
DEFALT_OVERLAP5 = {"jj": 0.0, "bb": 0.03, "jb": 0.16, "jm": 0.02, "je": 0.01, "jg": 0.99, "be": 0.01, "bm": 0.01, "bg": 0.16}
DEFALT_OVERLAP6 = {"jj": 0.0, "bb": 0.31, "jb": 0.47, "jm": 0.04, "je": 0.05, "jg": 0.12, "be": 0.02, "bm": 0.02, "bg": 0.02}
DEFALT_OVERLAP7 = {"jj": 0.0, "bb": 0.43, "jb": 0.57, "jm": 0.04, "je": 0.05, "jg": 0.06, "be": 0.03, "bm": 0.03, "bg": 0.02}
DEFALT_OVERLAP_TRIGGER = {1:DEFALT_OVERLAP1, 2:DEFALT_OVERLAP2, 3:DEFALT_OVERLAP3, 4:DEFALT_OVERLAP4, 5:DEFALT_OVERLAP5, 6:DEFALT_OVERLAP6, 7:DEFALT_OVERLAP7}

# Overlap Disablement Logic
if args.no_overlap:
    terminal.write("⚠️  Running WITHOUT channel overlaps (Independent Toys)\n")
    for d in DEFALT_OVERLAP_TRIGGER.values():
        for k in d:
            d[k] = 0.0
else:
    terminal.write("🔗 Running WITH channel overlaps\n")

mjjBinsL = [99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156]
mjjBins = np.array(mjjBinsL)
bin_centers = (mjjBins[1:] + mjjBins[:-1]) / 2.0

# =================================================================
# 3. Pre-Compute Expected Backgrounds & TF1 Templates
# =================================================================
expected_bkg, eval_centers_dict, edges_root_dict, tf1_dict = {}, {}, {}, {}

for TRIG_TYPE in range(1, 8):
    XMIN = {1:600, 2:300, 3:300, 4:500, 5:300, 6:1500, 7:600}.get(TRIG_TYPE, 300)
    mask = (bin_centers >= XMIN) & (bin_centers <= 9000)
    eval_centers_dict[TRIG_TYPE] = bin_centers[mask]
    
    idx = np.where(mask)[0]
    edges_crop = mjjBins[idx[0] : idx[-1] + 2]
    edges_root_dict[TRIG_TYPE] = array("d", edges_crop)
    
    expected_bkg[TRIG_TYPE] = {}
    tf1_dict[TRIG_TYPE] = {}
    
    for ch in CHANNELS:
        f = f"fits/fitme_p5_t{TRIG_TYPE}_{ch}.json"
        if not os.path.exists(f): continue
        with open(f, "r") as j:
            data_fit = json.load(j)
            params = [float(p) for p in data_fit["parameters"][:5]]
            expected_bkg[TRIG_TYPE][ch] = np.clip(FiveParam_NP(cms, eval_centers_dict[TRIG_TYPE], *params), 0, 1e7)
            
            if args.fit:
                back = ROOT.TF1(f"back_{TRIG_TYPE}_{ch}", FiveParam2015(), float(data_fit["fmin"]), float(data_fit["fmax"]), 5)
                for i, val in enumerate(params):
                    back.SetParameter(i, val)
                    if val == 0.0: back.FixParameter(i, val)
                tf1_dict[TRIG_TYPE][ch] = back

# =================================================================
# 4. Main Head-to-Head Loop
# =================================================================
terminal.write(f"🚀 Running Comparison (FitAgain={args.fit}, ChiMax={args.chimax})...\n")
terminal.write("-" * 55 + "\n")

time_pybh, time_fast, max_z_diff = 0.0, 0.0, 0.0
NrFound_pybh, NrFound_fast, fit_failures, loc_mismatches = 0, 0, 0, 0
np.random.seed(42)

for event in range(args.events):
    event_has_bump_pybh = False
    event_has_bump_fast = False
    
    if event > 0 and event % max(1, (args.events // 20)) == 0:
        progress = int((event / args.events) * 100)
        terminal.write(f"\rProgress: [{('=' * (progress//5)).ljust(20)}] {progress}% ")
        terminal.flush()

    for TRIG_TYPE in range(1, 8):
        if TRIG_TYPE not in expected_bkg or "jj" not in expected_bkg[TRIG_TYPE]: continue
        OVERLAP = DEFALT_OVERLAP_TRIGGER.get(TRIG_TYPE, DEFALT_OVERLAP1)
        edges_root = edges_root_dict[TRIG_TYPE]
        eval_centers = eval_centers_dict[TRIG_TYPE]
        
        exp_jj = expected_bkg[TRIG_TYPE]["jj"]
        pseudo_jj = np.random.poisson(exp_jj)
        residual_jj = (pseudo_jj - exp_jj) / np.maximum(exp_jj, 1.0)

        for channel in CHANNELS:
            if channel not in expected_bkg[TRIG_TYPE]: continue
            exp_ch = expected_bkg[TRIG_TYPE][channel]
            
            if channel == "jj":
                pseudo_ch = pseudo_jj
            else:
                frac = OVERLAP.get(channel, 0.0)
                overlap_pseudo = (exp_ch * frac) * (1.0 + residual_jj)
                indep_expected = np.clip(exp_ch - (exp_ch * frac), 0.0, 1e7) 
                pseudo_ch = np.clip(overlap_pseudo + np.random.poisson(indep_expected), 0.0, None)

            if np.sum(pseudo_ch) < MinEntries: continue

            # --- DYNAMIC REFITTING LOGIC ---
            active_bkg = exp_ch
            if args.fit:
                h_name = f"h_tmp_{event}_{TRIG_TYPE}_{channel}"
                h_tmp = ROOT.TH1D(h_name, h_name, len(edges_root)-1, edges_root)
                h_tmp.SetDirectory(0)
                
                for i, val in enumerate(pseudo_ch):
                    if val > 0:
                        h_tmp.SetBinContent(i+1, val)
                        h_tmp.SetBinError(i+1, ROOT.TMath.Sqrt(val))
                
                tf1_template = tf1_dict[TRIG_TYPE][channel].Clone()
                fit_result = h_tmp.Fit(tf1_template, "ISMR0Q") 
                
                ndf = tf1_template.GetNDF()
                chi2ndf = tf1_template.GetChisquare() / ndf if ndf > 0 else float('inf')
                
                if not fit_result.IsValid() or chi2ndf > args.chimax:
                    fit_failures += 1
                    continue
                
                active_bkg = np.clip(np.array([tf1_template.Eval(c) for c in eval_centers]), 0.0, 1e7)

            # --- Method 1: pyBumpHunter ---
            t0 = time.time()
            hunter = BH.BumpHunter1D(rang=None, width_min=2, width_max=30, npe=0, bins=edges_root)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hunter.bump_scan(data=pseudo_ch, bkg=active_bkg, is_hist=True, do_pseudo=False)
            time_pybh += (time.time() - t0)
            
            pybh_loc = hunter.min_loc_ar[0] if np.ndim(hunter.min_loc_ar) else hunter.min_loc_ar
            pybh_width = hunter.min_width_ar[0] if np.ndim(hunter.min_width_ar) else hunter.min_width_ar
            
            # --- Method 2: Fast Vectorized BumpHunter ---
            t0 = time.time()
            fast_tstat, fast_loc, fast_width = fast_bumphunter_stat_with_loc(pseudo_ch, active_bkg, max_width=30)
            time_fast += (time.time() - t0)
            
            # --- Comparisons ---
            pybh_pval = np.clip(hunter.min_Pval_ar[0], 1e-300, 1.0)
            fast_pval = np.clip(np.exp(-fast_tstat), 1e-300, 1.0)
            
            z_pybh = ROOT.RooStats.PValueToSignificance(pybh_pval) if pybh_pval < 1.0 else 0.0
            z_fast = ROOT.RooStats.PValueToSignificance(fast_pval) if fast_pval < 1.0 else 0.0
            
            max_z_diff = max(max_z_diff, abs(z_pybh - z_fast))
            
            if z_fast > 0 and (pybh_loc != fast_loc or pybh_width != fast_width): loc_mismatches += 1
            if z_pybh > args.zval: event_has_bump_pybh = True
            if z_fast > args.zval: event_has_bump_fast = True

    if event_has_bump_pybh: NrFound_pybh += 1
    if event_has_bump_fast: NrFound_fast += 1

# =================================================================
# 5. Final Output Comparison
# =================================================================
terminal.write(f"\rProgress: [{'=' * 20}] 100%\n\n")
terminal.write("=" * 55 + "\n")
terminal.write(" BUMPHUNTER METHOD COMPARISON RESULTS\n")
terminal.write("=" * 55 + "\n")
terminal.write(f"Total Pseudo-experiments : {args.events}\n")
terminal.write(f"Overlap Status           : {'Disabled' if args.no_overlap else 'Enabled'}\n")
terminal.write("-" * 55 + "\n")
terminal.write(f"Total Time (pyBumpHunter): {time_pybh:.4f} seconds\n")
terminal.write(f"Total Time (Fast NumPy)  : {time_fast:.4f} seconds\n")

if time_fast > 0:
    terminal.write(f"SPEEDUP MULTIPLIER       : {time_pybh / time_fast:.1f}x FASTER\n")

terminal.write("-" * 55 + "\n")
terminal.write(f"Toys w/ >{args.zval}σ Bumps (pyBH): {NrFound_pybh}\n")
terminal.write(f"Toys w/ >{args.zval}σ Bumps (Fast): {NrFound_fast}\n")
terminal.write(f"Rejected Fits (chi2 > {args.chimax}) : {fit_failures}\n")
terminal.write(f"Maximum Z-Score Divergence: {max_z_diff:.5e} σ\n")
terminal.write(f"Total Location Mismatches : {loc_mismatches}\n")

if max_z_diff < 1e-4 and loc_mismatches == 0:
    terminal.write("-> Mathematical equivalence & Coordinate matching VERIFIED.\n")
else:
    terminal.write("-> WARNING: Divergence or mismatch detected.\n")
terminal.write("=" * 55 + "\n")

log_file.close()
sys.stdout = terminal
