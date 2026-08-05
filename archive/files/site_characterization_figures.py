"""
=============================================================================
SITE CHARACTERIZATION FIGURES — APV-AD Study
Publication-quality figures for Section 3.1 / 3.2 of the article
=============================================================================

Generates seven figures describing the three semi-arid study sites:
  Fig S1 — Geographic map with site locations and climate zone polygons
  Fig S2 — Monthly climate profiles (GHI, T2m, ET0, precipitation) × 3 sites
  Fig S3 — Wind roses × 3 sites
  Fig S4 — Diurnal GHI profiles (seasonal composites) × 3 sites
  Fig S5 — Temperature exceedance frequency × 3 sites
  Fig S6 — Solar resource summary radar chart
  Fig S7 — Growing season calendar + crop-relevant climate metrics

All figures use a unified color palette and publication typography suitable
for Renewable & Sustainable Energy Reviews (RSER) / Applied Energy.

=============================================================================
USAGE
=============================================================================
1. Place PVGIS TMY CSV files in the same folder or adjust CLIMATE_FILES paths.
2. Run:  python site_characterization_figures.py
3. Figures saved to ./site_figures/ at 300 DPI (PDF + PNG).

If TMY files are not yet available, the script uses built-in synthetic TMY
data derived from the published PVGIS mean values in Table 4, so all figures
render immediately and can be replaced when real TMY data are loaded.
=============================================================================
DEPENDENCIES
=============================================================================
  numpy, pandas, matplotlib, scipy, pathlib
  Optional: cartopy (for geographic map Fig S1)
    Install: pip install cartopy
=============================================================================
"""

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────────────────────────────────────
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
from datetime import datetime
from scipy.ndimage import uniform_filter1d

# Optional cartopy
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _CARTOPY = True
except ImportError:
    _CARTOPY = False
    print("[INFO] cartopy not found — Fig S1 (map) will use a simplified version.")
    print("       Install with: pip install cartopy")


# =============================================================================
# SECTION 0 — GLOBAL STYLE & CONFIGURATION
# =============================================================================

# ── Publication color palette ─────────────────────────────────────────────────
# Three site colors: consistent across ALL figures
SITE_COLORS = {
    "Konya"       : "#1B6CA8",   # deep blue  — cold semi-arid (BSk)
    "Almeria"     : "#E05C2A",   # burnt orange — hot Mediterranean (BSh)
    "Ouagadougou" : "#2E9B5A",   # forest green — tropical Sahel (BSh)
}

SITE_MARKERS = {"Konya": "o", "Almeria": "s", "Ouagadougou": "^"}
SITE_LABELS  = {"Konya": "Konya (BSk, Turkey)",
                "Almeria": "Almería (BSh, Spain)",
                "Ouagadougou": "Ouagadougou (BSh, Burkina Faso)"}

# ── Secondary palette for climate variables ───────────────────────────────────
CLR_GHI    = "#F0A500"   # golden yellow — irradiance
CLR_TEMP   = "#C0392B"   # red — temperature
CLR_ET0    = "#8E44AD"   # purple — evapotranspiration
CLR_PRECIP = "#2980B9"   # blue — precipitation
CLR_WIND   = "#27AE60"   # green — wind

# ── Typography & figure style ─────────────────────────────────────────────────
mpl.rcParams.update({
    # Font — use Helvetica-like sans-serif stack for journal compatibility
    "font.family"        : "sans-serif",
    "font.sans-serif"    : ["DejaVu Sans", "Liberation Sans", "Arial"],
    "font.size"          : 9,
    "axes.titlesize"     : 10,
    "axes.labelsize"     : 9,
    "xtick.labelsize"    : 8,
    "ytick.labelsize"    : 8,
    "legend.fontsize"    : 8,
    "legend.framealpha"  : 0.92,
    "legend.edgecolor"   : "#CCCCCC",
    # Lines & spines
    "axes.linewidth"     : 0.7,
    "axes.spines.top"    : False,
    "axes.spines.right"  : False,
    "grid.linewidth"     : 0.4,
    "grid.alpha"         : 0.4,
    "grid.color"         : "#AAAAAA",
    # Figure
    "figure.dpi"         : 150,
    "savefig.dpi"        : 300,
    "savefig.bbox"       : "tight",
    "savefig.facecolor"  : "white",
    # Patches
    "patch.linewidth"    : 0.5,
})

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = Path("site_figures")
OUT_DIR.mkdir(exist_ok=True)

# ── Month labels ──────────────────────────────────────────────────────────────
MONTHS      = np.arange(1, 13)
MONTH_SHORT = ["J","F","M","A","M","J","J","A","S","O","N","D"]
MONTH_FULL  = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]
DAYS_PER_MONTH = np.array([31,28,31,30,31,30,31,31,30,31,30,31])


# =============================================================================
# SECTION 1 — SITE METADATA
# =============================================================================

SITES = {
    "Konya": {
        "lat": 37.87, "lon": 32.49, "elev": 1016,
        "koppen": "BSk", "country": "Turkey",
        "GHI_annual": 1650,
        "T_mean": 11.5, "T_summer_max": 38, "T_winter_min": -10,
        "ET0_annual": 1200, "precip_annual": 320,
        "irrig_req": 880,
        # Monthly means — from PVGIS-SARAH3 representative values
        "GHI_monthly": np.array([65,92,142,175,200,218,228,210,163,115,72,54]),
        "T_monthly"  : np.array([-1.5, 0.5, 5.5, 11.0, 16.5, 21.0,
                                  23.5, 23.0, 18.5, 12.5, 6.0, 1.0]),
        "ET0_monthly": np.array([25,32,65,90,120,148,162,155,110,70,35,22]),
        "precip_monthly": np.array([35,30,30,35,40,25,10,8,15,25,30,37]),
        "WS_monthly" : np.array([3.2,3.5,3.8,3.6,3.2,2.9,2.7,2.6,2.8,3.0,3.2,3.1]),
        "growing_start": (4, 25), "growing_end": (10, 30),
        "climate_file": "C:/Users/AMIDOU MAIGA/OneDrive - Institut 2IE/MAIGA Amidou/Erciyes University/Master/Thesis/Data/konya_climate.csv",
    },
    "Almeria": {
        "lat": 36.83, "lon": -2.46, "elev": 20,
        "koppen": "BSh", "country": "Spain",
        "GHI_annual": 1850,
        "T_mean": 18.5, "T_summer_max": 35, "T_winter_min": 5,
        "ET0_annual": 1350, "precip_annual": 200,
        "irrig_req": 1150,
        "GHI_monthly": np.array([80,110,158,195,218,238,248,232,178,130,88,68]),
        "T_monthly"  : np.array([11.0,12.0,14.0,16.5,19.5,23.5,
                                  26.5,27.0,24.0,19.5,15.0,12.0]),
        "ET0_monthly": np.array([40,52,82,108,135,162,180,168,120,85,50,38]),
        "precip_monthly": np.array([22,18,22,18,12,5,2,2,15,28,28,28]),
        "WS_monthly" : np.array([3.8,4.0,4.2,3.9,3.5,3.1,2.8,2.7,3.0,3.4,3.6,3.7]),
        "growing_start": (3, 1), "growing_end": (6, 30),
        "climate_file": "C:/Users/AMIDOU MAIGA/OneDrive - Institut 2IE/MAIGA Amidou/Erciyes University/Master/Thesis/Data/almeria_climate.csv",
    },
    "Ouagadougou": {
        "lat": 12.37, "lon": -1.52, "elev": 303,
        "koppen": "BSh", "country": "Burkina Faso",
        "GHI_annual": 2050,
        "T_mean": 28.5, "T_summer_max": 41, "T_winter_min": 15,
        "ET0_annual": 1800, "precip_annual": 800,
        "irrig_req": 1000,
        "GHI_monthly": np.array([190,200,210,205,188,168,155,158,170,185,192,188]),
        "T_monthly"  : np.array([25.0,28.0,31.5,33.5,32.5,29.0,
                                  27.5,27.0,27.5,28.5,27.5,25.5]),
        "ET0_monthly": np.array([155,160,190,190,168,140,128,130,138,150,145,148]),
        "precip_monthly": np.array([2, 4, 12, 28, 80, 130, 155, 175, 90, 42, 10, 2]),
        "WS_monthly" : np.array([2.8,3.2,3.5,3.2,2.9,2.5,2.2,2.0,2.1,2.3,2.5,2.6]),
        "growing_start": (7, 1), "growing_end": (10, 30),
        "climate_file": "C:/Users/AMIDOU MAIGA/OneDrive - Institut 2IE/MAIGA Amidou/Erciyes University/Master/Thesis/Data/ouagadougou_climate.csv",
    },
}


def load_tmy(site_name: str) -> pd.DataFrame | None:
    """
    Attempt to load PVGIS TMY CSV. Returns None if file not found.
    """
    fpath = Path(SITES[site_name]["climate_file"])
    if not fpath.exists():
        return None
    try:
        with open(fpath, "r", encoding="utf-8") as fh:
            lines = fh.readlines()
        header_idx = next(
            (i for i, ln in enumerate(lines) if "time(UTC)" in ln), 17)
        df = pd.read_csv(fpath, skiprows=header_idx)
        df.columns = df.columns.str.strip()
        rename = {}
        for c in df.columns:
            cl = c.lower()
            if "time"  in cl: rename[c] = "time"
            elif "t2m"  in cl: rename[c] = "T2m"
            elif "rh"   in cl: rename[c] = "RH"
            elif "g(h)" in cl: rename[c] = "GHI"
            elif "gb(n)"in cl: rename[c] = "DNI"
            elif "gd(h)"in cl: rename[c] = "DHI"
            elif "ws10m"in cl: rename[c] = "WS10m"
            elif "wd"   in cl: rename[c] = "WD"
        df = df.rename(columns=rename)
        df = df[df["time"].astype(str).str.match(r"\d{8}:\d{4}", na=False)]
        def _parse(t):
            t = str(t).strip()
            try:
                return datetime(int(t[0:4]),int(t[4:6]),int(t[6:8]),int(t[9:11]))
            except Exception: return None
        df["datetime"] = df["time"].apply(_parse)
        df = df.dropna(subset=["datetime"])
        for col in ["GHI","DNI","DHI","T2m","RH","WS10m"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
            else:
                df[col] = 0.0
        df["month"] = df["datetime"].dt.month
        df["hour"]  = df["datetime"].dt.hour
        return df.sort_values("datetime").reset_index(drop=True)
    except Exception as e:
        print(f"  [WARNING] Could not load {fpath}: {e}")
        return None


def _save(fig, name, tight=True):
    """Save figure in both PNG and PDF."""
    for ext in ["png", "pdf"]:
        p = OUT_DIR / f"{name}.{ext}"
        fig.savefig(p, format=ext)
    print(f"  Saved: {OUT_DIR / name}.png / .pdf")


# =============================================================================
# FIGURE S1 — Geographic location map
# =============================================================================

def fig_s1_map():
    """
    World map highlighting the three study sites with climate zone context.
    Uses cartopy if available, falls back to a clean matplotlib scatter map.
    """
    print("[Fig S1] Geographic map...")

    if _CARTOPY:
        fig = plt.figure(figsize=(9, 5))
        ax  = fig.add_subplot(1, 1, 1,
              projection=ccrs.PlateCarree(central_longitude=15))

        ax.set_extent([-20, 55, -5, 50], crs=ccrs.PlateCarree())

        ax.add_feature(cfeature.LAND,       facecolor="#F5F0E8", zorder=1)
        ax.add_feature(cfeature.OCEAN,      facecolor="#D6EAF8", zorder=1)
        ax.add_feature(cfeature.COASTLINE,  linewidth=0.4, edgecolor="#999999", zorder=2)
        ax.add_feature(cfeature.BORDERS,    linewidth=0.3, edgecolor="#BBBBBB",
                       linestyle="--", zorder=2)
        ax.add_feature(cfeature.RIVERS,     linewidth=0.3, edgecolor="#AED6F1",
                       alpha=0.5, zorder=2)

        # Semi-arid zone shading (approx 15–40°N strip)
        ax.add_patch(mpatches.Rectangle(
            (-20, 15), 75, 25,
            transform=ccrs.PlateCarree(),
            facecolor="#F9E79F", alpha=0.25, zorder=1,
            label="Semi-arid belt (BSk/BSh)"))

        # Site markers
        for name, sc in SITES.items():
            ax.plot(sc["lon"], sc["lat"],
                    marker=SITE_MARKERS[name], color=SITE_COLORS[name],
                    markersize=12, zorder=5, transform=ccrs.PlateCarree(),
                    markeredgecolor="white", markeredgewidth=1.2)
            offset = {"Konya": (2, 1.5), "Almeria": (-4, 2),
                      "Ouagadougou": (2, -3)}[name]
            ax.text(sc["lon"] + offset[0], sc["lat"] + offset[1],
                    name, fontsize=9, fontweight="bold",
                    color=SITE_COLORS[name],
                    transform=ccrs.PlateCarree(), zorder=6)

        gl = ax.gridlines(draw_labels=True, linewidth=0.3,
                          color="grey", alpha=0.5, linestyle="--")
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {"size": 7}

    else:
        # Fallback: simple lat-lon scatter
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.set_facecolor("#D6EAF8")
        ax.set_xlim(-20, 55); ax.set_ylim(-5, 50)

        # Continent outline simplified
        import matplotlib.patches as _mp
        ax.add_patch(_mp.FancyBboxPatch(
            (-5, 5), 45, 40, boxstyle="round,pad=2",
            facecolor="#F5F0E8", edgecolor="#999999", lw=0.5, zorder=1))

        ax.axhspan(15, 40, alpha=0.15, color="#F9E79F",
                   label="Semi-arid belt (BSk/BSh)", zorder=0)

        for name, sc in SITES.items():
            ax.scatter(sc["lon"], sc["lat"],
                       s=180, marker=SITE_MARKERS[name],
                       color=SITE_COLORS[name], zorder=5,
                       edgecolors="white", linewidths=1.5)
            offset = {"Konya": (2, 1.5), "Almeria": (-5, 2.5),
                      "Ouagadougou": (2.5, -3)}[name]
            ax.annotate(name,
                        xy=(sc["lon"], sc["lat"]),
                        xytext=(sc["lon"] + offset[0], sc["lat"] + offset[1]),
                        fontsize=9, fontweight="bold", color=SITE_COLORS[name],
                        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7))

        ax.set_xlabel("Longitude (°)", fontsize=9)
        ax.set_ylabel("Latitude (°)", fontsize=9)
        ax.grid(True, lw=0.3, alpha=0.5)
        ax.legend(loc="lower left", fontsize=8)

    # ── Legend panel ──────────────────────────────────────────────────────────
    handles = [
        mpatches.Patch(color=SITE_COLORS[n], label=SITE_LABELS[n])
        for n in SITES
    ]
    if _CARTOPY:
        ax.legend(handles=handles, loc="lower left", fontsize=8,
                  framealpha=0.9, title="Study sites", title_fontsize=8)

    fig.suptitle(
        "Figure S1. Location of the three APV-AD study sites across the semi-arid belt\n"
        "of Europe, the Middle East, and West Africa",
        fontsize=10, fontweight="bold", y=0.98)

    _save(fig, "FigS1_site_map")
    plt.close(fig)


# =============================================================================
# FIGURE S2 — Monthly climate profiles (4-panel per site)
# =============================================================================

def fig_s2_monthly_climate():
    """
    3 × 4 panel grid. Each column = one site. Each row = one variable:
      Row 1: Monthly GHI (kWh/m²) — bar chart
      Row 2: Monthly T2m (°C) — area chart with min/max band
      Row 3: Monthly ET₀ (mm) vs. Precipitation (mm) — water balance
      Row 4: Monthly wind speed (m/s) — line
    """
    print("[Fig S2] Monthly climate profiles...")

    fig = plt.figure(figsize=(13, 12))
    fig.suptitle(
        "Figure S2. Monthly climate profiles for the three study sites\n"
        "(PVGIS-SARAH3 TMY; horizontal bars = growing season)",
        fontsize=10, fontweight="bold", y=0.99)

    gs = gridspec.GridSpec(4, 3, figure=fig,
                           hspace=0.55, wspace=0.35,
                           top=0.94, bottom=0.06, left=0.07, right=0.97)

    row_labels = ["GHI (kWh/m²)", "Temperature (°C)",
                  "ET₀ & Precipitation (mm)", "Wind Speed (m/s)"]
    row_ylims  = [None, None, None, None]   # auto

    for col, (site_name, sc) in enumerate(SITES.items()):
        color = SITE_COLORS[site_name]

        # ── Determine growing season bar positions ─────────────────────────
        gs_start_m = sc["growing_start"][0]
        gs_end_m   = sc["growing_end"][0]
        if gs_end_m < gs_start_m:  # cross-year
            gs_months = list(range(gs_start_m, 13)) + list(range(1, gs_end_m + 1))
        else:
            gs_months = list(range(gs_start_m, gs_end_m + 1))

        # ── Row 0: GHI ────────────────────────────────────────────────────
        ax = fig.add_subplot(gs[0, col])
        bars = ax.bar(MONTHS, sc["GHI_monthly"],
                      color=[color if m in gs_months else "#CCCCCC" for m in MONTHS],
                      width=0.8, zorder=3, edgecolor="white", linewidth=0.3)
        ax.set_xticks(MONTHS)
        ax.set_xticklabels(MONTH_SHORT)
        ax.set_ylabel(row_labels[0] if col == 0 else "", fontsize=8)
        ax.set_title(f"{site_name}\n({sc['koppen']}, {sc['country']})",
                     fontsize=9, fontweight="bold", color=color)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(4))
        ann_y = sc["GHI_monthly"].max() * 0.95
        ax.text(0.97, 0.95, f"∑ {sc['GHI_annual']} kWh/m²",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color="#555555")
        ax.grid(axis="y", lw=0.3, alpha=0.5)
        # Growing season bracket
        if gs_months:
            ax.annotate("", xy=(max(gs_months) + 0.45, -16),
                        xytext=(min(gs_months) - 0.45, -16),
                        xycoords=("data", "axes fraction"),
                        textcoords=("data", "axes fraction"),
                        arrowprops=dict(
                            arrowstyle="<->", color=color, lw=1.2))

        # ── Row 1: Temperature ───────────────────────────────────────────
        ax = fig.add_subplot(gs[1, col])
        T   = sc["T_monthly"]
        ax.plot(MONTHS, T, color=color, lw=2, marker="o",
                markersize=4, zorder=4)
        ax.fill_between(MONTHS, T, 0, where=(T >= 0),
                        alpha=0.15, color=CLR_TEMP, zorder=2)
        ax.fill_between(MONTHS, T, 0, where=(T < 0),
                        alpha=0.15, color="#5DADE2", zorder=2)
        ax.axhline(0, color="#AAAAAA", lw=0.7, ls="--", zorder=1)
        # Threshold lines
        ax.axhline(10,  color="#F39C12", lw=0.5, ls=":", alpha=0.7, zorder=1)
        ax.axhline(26,  color="#E74C3C", lw=0.5, ls=":", alpha=0.7, zorder=1)
        ax.set_xticks(MONTHS); ax.set_xticklabels(MONTH_SHORT)
        ax.set_ylabel(row_labels[1] if col == 0 else "", fontsize=8)
        ax.text(0.97, 0.95, f"T̄ = {sc['T_mean']}°C",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color="#555555")
        ax.grid(axis="y", lw=0.3, alpha=0.5)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(5))

        # ── Row 2: Water balance (ET0 vs precip) ─────────────────────────
        ax = fig.add_subplot(gs[2, col])
        w   = 0.38
        ax.bar(MONTHS - w/2, sc["ET0_monthly"],    width=w, color=CLR_ET0,
               alpha=0.85, label="ET₀", zorder=3, edgecolor="white", lw=0.2)
        ax.bar(MONTHS + w/2, sc["precip_monthly"], width=w, color=CLR_PRECIP,
               alpha=0.85, label="Precip.", zorder=3, edgecolor="white", lw=0.2)
        # Deficit shading
        deficit = sc["ET0_monthly"] - sc["precip_monthly"]
        ax2 = ax.twinx()
        ax2.fill_between(MONTHS, deficit, 0, where=(deficit > 0),
                         alpha=0.12, color=CLR_TEMP, label="Water deficit")
        ax2.set_ylabel("Deficit (mm)" if col == 2 else "", fontsize=7,
                       color=CLR_TEMP)
        ax2.tick_params(axis="y", labelcolor=CLR_TEMP, labelsize=7)
        ax2.yaxis.set_major_locator(mticker.MaxNLocator(4))
        ax.set_xticks(MONTHS); ax.set_xticklabels(MONTH_SHORT)
        ax.set_ylabel(row_labels[2] if col == 0 else "", fontsize=8)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(4))
        ax.grid(axis="y", lw=0.3, alpha=0.5)
        if col == 0:
            ax.legend(fontsize=7, loc="upper left", ncol=2,
                      framealpha=0.85, handlelength=1.2)
        ax.text(0.97, 0.95,
                f"P={sc['precip_annual']} | ET₀={sc['ET0_annual']} mm",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7, color="#555555")

        # ── Row 3: Wind speed ────────────────────────────────────────────
        ax = fig.add_subplot(gs[3, col])
        WS = sc["WS_monthly"]
        ax.plot(MONTHS, WS, color=CLR_WIND, lw=1.8, marker="^",
                markersize=4, zorder=4)
        ax.fill_between(MONTHS, WS, WS.min() - 0.3,
                        alpha=0.12, color=CLR_WIND, zorder=2)
        ax.set_xticks(MONTHS); ax.set_xticklabels(MONTH_SHORT)
        ax.set_ylabel(row_labels[3] if col == 0 else "", fontsize=8)
        ax.set_ylim(bottom=0)
        ax.text(0.97, 0.95, f"WS̄ = {WS.mean():.1f} m/s",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color="#555555")
        ax.grid(axis="y", lw=0.3, alpha=0.5)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(4))
        ax.set_xlabel("Month", fontsize=8)

    # ── Row labels on the left ─────────────────────────────────────────────
    row_titles = ["Solar\nIrradiance", "Air\nTemperature",
                  "Water\nBalance", "Wind\nSpeed"]
    for i, rt in enumerate(row_titles):
        fig.text(0.005, 0.87 - i * 0.222, rt,
                 fontsize=8, fontweight="bold", color="#444444",
                 va="center", rotation=90)

    _save(fig, "FigS2_monthly_climate")
    plt.close(fig)


# =============================================================================
# FIGURE S3 — Diurnal GHI profiles (seasonal composites)
# =============================================================================

def fig_s3_diurnal_ghi():
    """
    For each site: diurnal GHI curve for 4 seasons (DJF, MAM, JJA, SON).
    Uses actual TMY data if available, otherwise uses a theoretical Gaussian.
    """
    print("[Fig S3] Diurnal GHI profiles...")

    SEASON_COLORS = {
        "DJF": "#5DADE2",  # blue
        "MAM": "#82E0AA",  # light green
        "JJA": "#F39C12",  # orange
        "SON": "#C0392B",  # red
    }
    SEASON_MONTHS = {
        "DJF": [12, 1, 2], "MAM": [3, 4, 5],
        "JJA": [6, 7, 8],  "SON": [9, 10, 11]
    }
    SEASON_LABELS = {
        "DJF": "Dec–Feb (Winter)", "MAM": "Mar–May (Spring)",
        "JJA": "Jun–Aug (Summer)", "SON": "Sep–Nov (Autumn)"
    }

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5),
                             sharey=False, gridspec_kw={"wspace": 0.35})
    fig.suptitle(
        "Figure S3. Mean diurnal GHI profiles by season for the three study sites\n"
        "(shaded band = ±1 standard deviation)",
        fontsize=10, fontweight="bold")

    hours = np.arange(0, 24)

    for ax, (site_name, sc) in zip(axes, SITES.items()):
        color = SITE_COLORS[site_name]
        tmy   = load_tmy(site_name)

        for season, s_months in SEASON_MONTHS.items():
            if tmy is not None:
                mask = tmy["month"].isin(s_months)
                grp  = tmy[mask].groupby("hour")["GHI"]
                mean_ghi = grp.mean().reindex(hours, fill_value=0).values
                std_ghi  = grp.std().reindex(hours, fill_value=0).values
            else:
                # Synthetic: bell curve scaled by monthly GHI mean
                peak_month_idx = [m - 1 for m in s_months]
                peak_ghi = np.mean([sc["GHI_monthly"][i] for i in peak_month_idx])
                # Convert monthly kWh/m² to peak W/m² (approx 10 daylight hrs)
                peak_W   = peak_ghi * 1000 / (10 * 30)
                mean_ghi = np.zeros(24)
                std_ghi  = np.zeros(24)
                for h in range(6, 19):
                    x = (h - 12.0) / 3.5
                    mean_ghi[h] = peak_W * np.exp(-0.5 * x**2)
                    std_ghi[h]  = mean_ghi[h] * 0.15

            sc_col = SEASON_COLORS[season]
            ax.plot(hours, mean_ghi, color=sc_col, lw=1.8,
                    label=SEASON_LABELS[season], zorder=4)
            ax.fill_between(hours,
                            np.maximum(0, mean_ghi - std_ghi),
                            mean_ghi + std_ghi,
                            color=sc_col, alpha=0.12, zorder=2)

        # Peak annotation
        all_peaks = []
        for s_months in SEASON_MONTHS.values():
            if tmy is not None:
                mask = tmy["month"].isin(s_months)
                grp  = tmy[mask].groupby("hour")["GHI"]
                mean_g = grp.mean().reindex(hours, fill_value=0).values
            else:
                peak_month_idx = [m - 1 for m in s_months]
                peak_ghi = np.mean([sc["GHI_monthly"][i] for i in peak_month_idx])
                peak_W   = peak_ghi * 1000 / (10 * 30)
                mean_g   = np.array([peak_W * np.exp(-0.5 * ((h-12)/3.5)**2)
                                     if 6 <= h < 19 else 0 for h in hours])
            all_peaks.append(mean_g.max())

        ax.set_title(f"{site_name}\n({sc['koppen']})",
                     fontsize=9, fontweight="bold", color=color)
        ax.set_xlabel("Hour of day (UTC+local)", fontsize=8)
        ax.set_ylabel("GHI (W/m²)" if ax == axes[0] else "", fontsize=8)
        ax.set_xlim(0, 23)
        ax.set_ylim(0)
        ax.set_xticks([0, 4, 8, 12, 16, 20])
        ax.axvline(12, color="#AAAAAA", lw=0.5, ls="--", alpha=0.7)
        ax.text(12.3, ax.get_ylim()[1] * 0.95, "solar noon",
                fontsize=7, color="#AAAAAA", va="top")
        ax.grid(axis="y", lw=0.3, alpha=0.5)
        if ax == axes[0]:
            ax.legend(fontsize=7.5, loc="upper left",
                      framealpha=0.9, handlelength=1.5)

    _save(fig, "FigS3_diurnal_GHI")
    plt.close(fig)


# =============================================================================
# FIGURE S4 — Temperature frequency distribution & extremes
# =============================================================================

def fig_s4_temperature_analysis():
    """
    For each site:
      Left: Monthly T distribution (violin + box)
      Right: Cumulative exceedance curve for T2m showing critical thresholds
    """
    print("[Fig S4] Temperature analysis...")

    THRESH_COLORS = {
        "Tomato optimum\n(18–26°C)": ("#27AE60", 18, 26),
        "Heat stress\n(>32°C)":      ("#E74C3C", 32, 50),
        "Frost risk\n(<2°C)":        ("#5DADE2", -20, 2),
        "Digester min\n(26°C)":      ("#8E44AD", 26, None),
    }

    fig = plt.figure(figsize=(13, 5.5))
    fig.suptitle(
        "Figure S4. Air temperature characterization: monthly distributions "
        "and threshold exceedance frequencies",
        fontsize=10, fontweight="bold")

    gs_outer = gridspec.GridSpec(1, 3, figure=fig,
                                 wspace=0.38, left=0.06, right=0.97,
                                 top=0.88, bottom=0.10)

    for col, (site_name, sc) in enumerate(SITES.items()):
        color = SITE_COLORS[site_name]
        tmy   = load_tmy(site_name)

        gs_inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=gs_outer[col], hspace=0.45)

        # ── Panel A: Monthly T ─────────────────────────────────────────────
        ax_top = fig.add_subplot(gs_inner[0])
        T_m = sc["T_monthly"]

        if tmy is not None:
            # Violin for real TMY
            data_by_month = [tmy[tmy["month"] == m]["T2m"].values
                             for m in range(1, 13)]
            parts = ax_top.violinplot(data_by_month, positions=MONTHS,
                                      widths=0.7, showmedians=True,
                                      showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(color); pc.set_alpha(0.35)
            parts["cmedians"].set_color(color); parts["cmedians"].set_lw(1.5)
        else:
            # Approximate from monthly means ± seasonal std
            ax_top.bar(MONTHS, T_m,
                       color=[color if t > 10 else "#5DADE2" for t in T_m],
                       alpha=0.6, width=0.7, zorder=3)

        ax_top.plot(MONTHS, T_m, color=color, lw=1.5, zorder=5,
                    marker="o", markersize=3.5)
        ax_top.axhline(0,  color="#AAAAAA", lw=0.6, ls="--", zorder=1)
        ax_top.axhline(26, color="#8E44AD", lw=0.6, ls=":", alpha=0.7, zorder=1)
        ax_top.axhspan(18, 26, alpha=0.07, color="#27AE60", zorder=0)
        ax_top.set_xticks(MONTHS); ax_top.set_xticklabels(MONTH_SHORT)
        ax_top.set_ylabel("T2m (°C)" if col == 0 else "", fontsize=8)
        ax_top.set_title(f"{site_name} ({sc['koppen']})",
                         fontsize=9, fontweight="bold", color=color)
        ax_top.grid(axis="y", lw=0.3, alpha=0.5)
        ax_top.yaxis.set_major_locator(mticker.MaxNLocator(5))

        # ── Panel B: Exceedance curve ─────────────────────────────────────
        ax_bot = fig.add_subplot(gs_inner[1])

        if tmy is not None:
            T_all = tmy["T2m"].values
        else:
            # Synthesize hourly from monthly means + diurnal ±7°C
            T_all = np.concatenate([
                np.random.normal(t, 5, d * 24)
                for t, d in zip(sc["T_monthly"], DAYS_PER_MONTH)
            ])

        T_sorted = np.sort(T_all)
        exceed_pct = 100 * (1 - np.arange(len(T_sorted)) / len(T_sorted))

        ax_bot.plot(T_sorted, exceed_pct, color=color, lw=1.8, zorder=4)

        # Threshold shading
        thresh = [
            (18, 26, "#27AE60", "Tomato optimum\n18–26°C"),
            (32, 50, "#E74C3C", "Heat stress >32°C"),
            (-25, 2, "#5DADE2", "Frost risk <2°C"),
        ]
        for t_lo, t_hi, tc, tl in thresh:
            mask = (T_sorted >= t_lo) & (T_sorted <= t_hi)
            if mask.any():
                ax_bot.fill_betweenx(exceed_pct[mask],
                                     T_sorted[mask].min(), T_sorted[mask].max(),
                                     alpha=0.12, color=tc, zorder=2)
        ax_bot.axvline(26, color="#8E44AD", lw=0.7, ls=":", alpha=0.8, zorder=3)
        ax_bot.set_xlabel("T2m (°C)", fontsize=8)
        ax_bot.set_ylabel("Exceedance (%)" if col == 0 else "", fontsize=8)
        ax_bot.set_xlim(T_sorted.min() - 2, T_sorted.max() + 2)
        ax_bot.set_ylim(0, 100)
        ax_bot.grid(lw=0.3, alpha=0.5)

        # Fraction annotations
        frost_pct  = 100 * np.mean(T_all < 2)
        optim_pct  = 100 * np.mean((T_all >= 18) & (T_all <= 26))
        stress_pct = 100 * np.mean(T_all > 32)
        ax_bot.text(0.03, 0.97,
                    f"<2°C: {frost_pct:.0f}%\n18–26°C: {optim_pct:.0f}%\n>32°C: {stress_pct:.0f}%",
                    transform=ax_bot.transAxes, va="top",
                    fontsize=7, color="#444444",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec="#CCCCCC", alpha=0.9))

    # ── Shared legend ─────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color="#27AE60", alpha=0.5, label="Tomato optimum (18–26°C)"),
        mpatches.Patch(color="#E74C3C", alpha=0.5, label="Heat stress (>32°C)"),
        mpatches.Patch(color="#5DADE2", alpha=0.5, label="Frost risk (<2°C)"),
        mpatches.Patch(color="#8E44AD", alpha=0.5, label="Digester min. T = 26°C"),
    ]
    fig.legend(handles=legend_patches, loc="lower center",
               ncol=4, fontsize=7.5, framealpha=0.9,
               bbox_to_anchor=(0.5, 0.01))

    _save(fig, "FigS4_temperature_analysis")
    plt.close(fig)


# =============================================================================
# FIGURE S5 — Solar resource summary: GHI heatmap (day × hour)
# =============================================================================

def fig_s5_ghi_heatmap():
    """
    Day-of-year × hour-of-day GHI heatmap for each site.
    Clearly shows seasonal solar availability and growing season window.
    """
    print("[Fig S5] GHI heatmap (day × hour)...")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5),
                             gridspec_kw={"wspace": 0.35})
    fig.suptitle(
        "Figure S5. Day-of-year × hour-of-day GHI heatmap (W/m²)\n"
        "Horizontal lines = growing season boundaries; "
        "white contour = 400 W/m² threshold",
        fontsize=10, fontweight="bold")

    cmap = LinearSegmentedColormap.from_list(
        "solar", ["#0D1B2A", "#1A4A8A", "#F0A500", "#FFEE58", "#FFFFFF"], N=256)

    for ax, (site_name, sc) in zip(axes, SITES.items()):
        color = SITE_COLORS[site_name]
        tmy   = load_tmy(site_name)

        if tmy is not None:
            tmy["doy"] = tmy["datetime"].dt.dayofyear
            pivot = tmy.pivot_table(
                values="GHI", index="doy", columns="hour",
                aggfunc="mean").reindex(range(1, 366), fill_value=0)
            Z = pivot.values.T   # shape (24, 365)
        else:
            # Synthesize from monthly GHI data
            Z = np.zeros((24, 365))
            day = 0
            for m_idx, (n_days, ghi_m) in enumerate(
                    zip(DAYS_PER_MONTH, sc["GHI_monthly"])):
                peak_W = ghi_m * 1000 / (10 * n_days)
                for d in range(n_days):
                    doy = day + d
                    if doy >= 365: break
                    for h in range(24):
                        if 5 <= h <= 19:
                            x = (h - 12.0) / 3.2
                            Z[h, doy] = max(0, peak_W * np.exp(-0.5 * x**2)
                                           * (0.9 + 0.1 * np.random.rand()))
                day += n_days

        im = ax.imshow(Z, aspect="auto", origin="lower",
                       extent=[1, 365, 0, 23], cmap=cmap, vmin=0, vmax=900)
        ax.contour(np.arange(1, Z.shape[1] + 1), np.arange(Z.shape[0]),
                   Z, levels=[400], colors="white", linewidths=0.5, alpha=0.7)

        # Growing season lines
        gs_start_m, gs_start_d = sc["growing_start"]
        gs_end_m,   gs_end_d   = sc["growing_end"]
        gs_start_doy = sum(DAYS_PER_MONTH[:gs_start_m - 1]) + gs_start_d
        gs_end_doy   = sum(DAYS_PER_MONTH[:gs_end_m   - 1]) + gs_end_d
        ax.axvline(gs_start_doy, color=color, lw=1.5, ls="--", alpha=0.9)
        ax.axvline(gs_end_doy,   color=color, lw=1.5, ls="--", alpha=0.9)
        ax.text(gs_start_doy + 3, 21, "Planting", fontsize=6.5, color=color,
                rotation=90, va="top")
        ax.text(gs_end_doy + 3,   21, "Harvest",  fontsize=6.5, color=color,
                rotation=90, va="top")

        ax.set_title(f"{site_name}\n({sc['koppen']})",
                     fontsize=9, fontweight="bold", color=color)
        ax.set_xlabel("Day of year", fontsize=8)
        ax.set_ylabel("Hour of day" if ax == axes[0] else "", fontsize=8)
        ax.set_yticks([0, 6, 12, 18, 23])
        ax.set_yticklabels(["00:00","06:00","12:00","18:00","23:00"],
                           fontsize=7)

        # Month tick marks on x-axis
        month_doys = np.cumsum([0] + list(DAYS_PER_MONTH[:-1])) + 15
        ax.set_xticks(month_doys)
        ax.set_xticklabels(MONTH_SHORT, fontsize=7)

    # ── Shared colorbar ───────────────────────────────────────────────────
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.65])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("GHI (W/m²)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.subplots_adjust(right=0.90)
    _save(fig, "FigS5_GHI_heatmap")
    plt.close(fig)


# =============================================================================
# FIGURE S6 — Multi-variable radar / spider chart
# =============================================================================

def fig_s6_radar():
    """
    Radar chart comparing the three sites across 8 normalized dimensions:
    GHI, T_mean, T_variability, ET₀, Precipitation, Irrigation requirement,
    PV potential (GHI), and Digester heating challenge (1/T_winter_min scaled).
    """
    print("[Fig S6] Radar chart...")

    categories = [
        "Solar resource\n(GHI)",
        "Mean temperature\n(T̄)",
        "Summer heat\nstress",
        "Reference ET₀",
        "Precipitation",
        "Irrigation\nrequirement",
        "Winter heating\nchallenge",
        "Growing\nseason length",
    ]
    N = len(categories)

    # Raw values — normalize to [0, 1] across sites
    raw = {
        "Konya"       : [1650, 11.5, 38, 1200, 320,  880, 10, 189],
        "Almeria"     : [1850, 18.5, 35, 1350, 200, 1150,  5, 122],
        "Ouagadougou" : [2050, 28.5, 41, 1800, 800, 1000, 15, 122],
    }
    raw_arr = np.array(list(raw.values()), dtype=float)
    # For "winter heating challenge": higher T_winter_min → LESS challenge
    # So invert: challenge = max - T_winter_min → [15-(-10)=25, 15-5=10, 15-15=0]
    winter_idx = 6
    raw_arr[:, winter_idx] = raw_arr[:, winter_idx].max() - raw_arr[:, winter_idx] + 1

    # Normalize to [0.2, 1.0] for each dimension
    col_min = raw_arr.min(axis=0)
    col_max = raw_arr.max(axis=0)
    normed  = 0.2 + 0.8 * (raw_arr - col_min) / (col_max - col_min + 1e-9)

    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]   # close the loop

    fig, ax = plt.subplots(figsize=(7, 6.5),
                           subplot_kw={"polar": True})
    fig.suptitle(
        "Figure S6. Normalized site comparison radar chart\n"
        "(outer = higher value; 8 climatic and agricultural dimensions)",
        fontsize=10, fontweight="bold")

    ax.set_facecolor("#F8F9FA")
    fig.patch.set_facecolor("white")

    # Background grid circles
    for r in [0.2, 0.4, 0.6, 0.8, 1.0]:
        ax.plot(angles, [r] * (N + 1),
                color="#DDDDDD", lw=0.5, ls="-", zorder=0)

    # Plot each site
    for idx, (site_name, sc) in enumerate(SITES.items()):
        color  = SITE_COLORS[site_name]
        values = normed[idx].tolist()
        values += values[:1]  # close loop
        ax.plot(angles, values, color=color, lw=2.0, zorder=4,
                marker="o", markersize=6)
        ax.fill(angles, values, color=color, alpha=0.12, zorder=3)

    # Category labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=7.5, color="#333333")
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2","0.4","0.6","0.8","1.0"], fontsize=6.5,
                       color="#AAAAAA")
    ax.set_rlim(0, 1.1)
    ax.spines["polar"].set_visible(False)
    ax.grid(False)

    # Legend
    legend_handles = [
        mpatches.Patch(color=SITE_COLORS[n], label=SITE_LABELS[n])
        for n in SITES
    ]
    ax.legend(handles=legend_handles,
              loc="lower center", bbox_to_anchor=(0.5, -0.18),
              ncol=1, fontsize=8, framealpha=0.9,
              handlelength=1.5, handleheight=1.2)

    _save(fig, "FigS6_radar_comparison")
    plt.close(fig)


# =============================================================================
# FIGURE S7 — Growing season calendar & crop-climate compatibility
# =============================================================================

def fig_s7_growing_season():
    """
    Gantt-style calendar showing growing season and key climate constraints.
    Overlaid: monthly T bar (colored by stress level) + ET₀-precip balance.
    """
    print("[Fig S7] Growing season calendar...")

    fig, axes = plt.subplots(3, 1, figsize=(12, 8),
                             gridspec_kw={"hspace": 0.55,
                                          "top": 0.90, "bottom": 0.08,
                                          "left": 0.12, "right": 0.96})
    fig.suptitle(
        "Figure S7. Growing season calendar, monthly temperature stress levels,\n"
        "and water balance for tomato cultivation at the three study sites",
        fontsize=10, fontweight="bold")

    for ax, (site_name, sc) in zip(axes, SITES.items()):
        color = SITE_COLORS[site_name]
        T_m   = sc["T_monthly"]
        P_m   = sc["precip_monthly"]
        ET_m  = sc["ET0_monthly"]

        # ── Color bars by temperature stress ─────────────────────────────
        bar_colors = []
        for t in T_m:
            if t < 5:     bar_colors.append("#5DADE2")  # frost/very cold
            elif t < 10:  bar_colors.append("#AED6F1")  # cool
            elif t < 18:  bar_colors.append("#82E0AA")  # sub-optimal
            elif t <= 26: bar_colors.append("#27AE60")  # optimal
            elif t <= 32: bar_colors.append("#F39C12")  # warm stress
            else:         bar_colors.append("#E74C3C")  # heat stress

        ax.bar(MONTHS, T_m, color=bar_colors, width=0.7,
               zorder=3, alpha=0.75, edgecolor="white", lw=0.3)
        ax.plot(MONTHS, T_m, color="#333333", lw=1.2,
                marker="o", markersize=3.5, zorder=5)

        # ── Water balance (thin stacked bars below T axis) ────────────────
        ax2 = ax.twinx()
        deficit = np.maximum(0, ET_m - P_m)
        surplus = np.maximum(0, P_m - ET_m)
        ax2.bar(MONTHS, -deficit, width=0.7, color=CLR_TEMP,
                alpha=0.18, zorder=2, bottom=0)
        ax2.bar(MONTHS,  surplus, width=0.7, color=CLR_PRECIP,
                alpha=0.18, zorder=2, bottom=0)
        ax2.axhline(0, color="#AAAAAA", lw=0.5, ls="-", zorder=1)
        ax2.set_ylabel("Water balance (mm)", fontsize=7, color="#777777")
        ax2.tick_params(axis="y", labelcolor="#777777", labelsize=7)
        ax2.yaxis.set_major_locator(mticker.MaxNLocator(4))

        # ── Growing season band ───────────────────────────────────────────
        gs_m = sc["growing_start"][0]
        ge_m = sc["growing_end"][0]
        if ge_m >= gs_m:
            ax.axvspan(gs_m - 0.4, ge_m + 0.4, alpha=0.10,
                       color=color, zorder=0)
            ax.annotate("", xy=(ge_m + 0.45, ax.get_ylim()[1] * 0.85),
                        xytext=(gs_m - 0.45, ax.get_ylim()[1] * 0.85),
                        arrowprops=dict(arrowstyle="<->",
                                        color=color, lw=1.5))
            gs_days = sc["growing_end"][1] - sc["growing_start"][1] + \
                      sum(DAYS_PER_MONTH[gs_m - 1:ge_m])
            mid_m   = (gs_m + ge_m) / 2
            ax.text(mid_m, ax.get_ylim()[1] * 0.90,
                    f"Growing season\n({abs(gs_days)} days)",
                    ha="center", fontsize=7.5, color=color, fontweight="bold",
                    va="bottom")

        # ── APV shading benefit annotation ───────────────────────────────
        stress_months = [i + 1 for i, t in enumerate(T_m) if t > 30]
        for sm in stress_months:
            if gs_m <= sm <= ge_m or (ge_m < gs_m and (sm >= gs_m or sm <= ge_m)):
                ax.annotate("APV\nshading\nbenefit",
                            xy=(sm, T_m[sm - 1]),
                            xytext=(sm, T_m[sm - 1] + 5),
                            fontsize=5.5, color="#E74C3C", ha="center",
                            arrowprops=dict(arrowstyle="->",
                                            color="#E74C3C", lw=0.7))

        ax.set_title(f"{site_name} ({sc['koppen']}, {sc['country']})",
                     fontsize=9, fontweight="bold", color=color)
        ax.set_ylabel("Mean T (°C)", fontsize=8)
        ax.set_xticks(MONTHS)
        ax.set_xticklabels(MONTH_FULL, fontsize=7.5)
        ax.grid(axis="y", lw=0.3, alpha=0.5, zorder=0)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(5))

        # ── Stress legend ─────────────────────────────────────────────────
        stress_legend = [
            mpatches.Patch(color="#5DADE2", label="<5°C (frost)"),
            mpatches.Patch(color="#82E0AA", label="10–18°C (sub-opt.)"),
            mpatches.Patch(color="#27AE60", label="18–26°C (optimal)"),
            mpatches.Patch(color="#F39C12", label="26–32°C (warm)"),
            mpatches.Patch(color="#E74C3C", label=">32°C (heat stress)"),
        ]
        if ax == axes[0]:
            ax.legend(handles=stress_legend,
                      loc="upper right", fontsize=6.5, ncol=5,
                      framealpha=0.9, handlelength=1.2,
                      title="Temperature stress class", title_fontsize=7)

        ax.set_xlabel("Month" if ax == axes[-1] else "", fontsize=8)

    _save(fig, "FigS7_growing_season_calendar")
    plt.close(fig)


# =============================================================================
# FIGURE S8 — Summary statistics table figure (publication table)
# =============================================================================

def fig_s8_summary_table():
    """
    Renders the key site parameters as a formatted matplotlib table figure
    suitable for direct inclusion in the paper (matches Table 4 content).
    """
    print("[Fig S8] Summary statistics table...")

    rows = [
        ["Parameter",              "Konya (BSk)",  "Almería (BSh)", "Ouagadougou (BSh)"],
        ["Latitude / Longitude",   "37.87°N, 32.49°E", "36.83°N, 2.46°W", "12.37°N, 1.52°W"],
        ["Elevation (m a.s.l.)",   "1,016",         "20",              "303"],
        ["Annual GHI (kWh/m²)",    "1,650",         "1,850",           "2,050"],
        ["Mean annual T (°C)",     "11.5",          "18.5",            "28.5"],
        ["Summer T_max (°C)",      "38",            "35",              "41"],
        ["Winter T_min (°C)",      "−10",           "+5",              "+15"],
        ["Annual ET₀ (mm)",        "1,200",         "1,350",           "1,800"],
        ["Annual precipitation (mm)", "320",        "200",             "800"],
        ["Irrigation requirement (mm)", "880",      "1,150",           "1,000"],
        ["Growing season",         "25 Apr – 30 Oct", "1 Mar – 30 Jun", "1 Jul – 30 Oct"],
        ["Growing season (days)",  "189",           "122",             "122"],
        ["Tomato price (USD/t)",   "300",           "400",             "250"],
        ["Electricity price (USD/kWh)", "0.06",     "0.09",            "0.12"],
    ]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_axis_off()
    fig.suptitle(
        "Figure S8. Summary of site characteristics and model parameters\n"
        "(reproduced from Table 4–5 of the article for figure-panel compatibility)",
        fontsize=10, fontweight="bold")

    col_widths = [0.30, 0.22, 0.22, 0.26]
    header_color = "#1B4F72"
    row_colors_alt = ["#EBF5FB", "white"]
    site_colors_header = ["white",
                          SITE_COLORS["Konya"],
                          SITE_COLORS["Almeria"],
                          SITE_COLORS["Ouagadougou"]]

    table = ax.table(
        cellText   = rows,
        colWidths  = col_widths,
        loc        = "center",
        cellLoc    = "center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.55)

    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#CCCCCC")
        if row == 0:
            cell.set_facecolor(site_colors_header[col])
            cell.set_text_props(color="white" if col > 0 else "white",
                                fontweight="bold", fontsize=9)
            if col == 0:
                cell.set_facecolor(header_color)
        else:
            if col == 0:
                cell.set_facecolor("#EBF5FB")
                cell.set_text_props(fontweight="bold", color="#1B4F72",
                                    ha="left")
                cell.PAD = 0.04
            else:
                cell.set_facecolor(row_colors_alt[row % 2])
                cell.set_text_props(
                    color=site_colors_header[col],
                    fontweight="normal")

    _save(fig, "FigS8_summary_table")
    plt.close(fig)


# =============================================================================
# MAIN — Run all figures
# =============================================================================

def main():
    print("=" * 60)
    print("  APV-AD Site Characterization Figures")
    print("  Output directory:", OUT_DIR.resolve())
    print("=" * 60)

    fig_s1_map()
    fig_s2_monthly_climate()
    fig_s3_diurnal_ghi()
    fig_s4_temperature_analysis()
    fig_s5_ghi_heatmap()
    fig_s6_radar()
    fig_s7_growing_season()
    fig_s8_summary_table()

    print("\n" + "=" * 60)
    print(f"  All figures saved to: {OUT_DIR.resolve()}")
    print("  Formats: PNG (300 DPI) + PDF (vector)")
    print("=" * 60)


if __name__ == "__main__":
    main()
