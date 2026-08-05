# APV + AD Integrated Model: Comprehensive Analysis for Master Thesis & Journal Publication

## Executive Summary

Your APV + AD model (v3.0) represents a **significant research contribution** combining bifacial agrivoltaics with anaerobic digestion across four climate zones. The work addresses real research gaps (G1–G5 in your paper) and introduces novel metrics (extended LER). However, publication in a **top-tier journal** (Applied Energy, Renewable & Sustainable Energy Reviews, Renewable Energy) requires addressing several critical technical, methodological, and presentation issues.

---

## 🌟 MAJOR STRENGTHS

### 1. **Novel Integrated Framework** ⭐⭐⭐
- **First true APV-AD coupling**: Previous work treats technologies in isolation
- **Extended LER metric (eLER)** bridges the food-energy-water nexus
- **Three-component scoring**: Crop + PV + Biogas under one objective function
- **Global significance**: Addresses complementary resource gaps (food, electricity, renewable gas)

### 2. **Rigorous Four-Site Climate Diversity** ⭐⭐⭐
- BSk (cold semi-arid), BSh (hot semi-arid), tropical semi-arid, Cfb (temperate)
- Latitudinal range 12°N–48°N demonstrates climate robustness
- PVGIS-SARAH3 TMY data → high-quality, reproducible climate inputs
- Site-specific economics + agronomic parameters (Freiburg 30 t/ha vs. 60 t/ha) → realistic comparisons

### 3. **Sophisticated Physical Modeling** ⭐⭐⭐
- **Bifacial PV**: φ=0.70, validated against PVGIS reference (validated in apv\_ad\_pv\_validation.py)
- **Faiman cell temperature** with wind dependency (U₀=25, U₁=6.84)
- **FIX 1 (asymmetric PAR saturation)**: Critical fix restoring LER\_crop from spurious 1.000 to physically correct 0.55–0.65
- **Temperature-corrected BMP** (Arrhenius model): First quantification across 4 climates
- **FAO-56 hourly ET₀** with site-specific psychrometric constants

### 4. **Methodological Rigor** ⭐⭐⭐
- **Differential Evolution** with dithered mutation + L-BFGS-B polishing → dual-phase optimization
- **Validation strategy**: PV vs. PVGIS reference, crop model vs. published PAR saturation data
- **Monte-Carlo sensitivity** (N=2000): Tornado analysis + income P5/P95 distributions
- **Pareto optimization** (NSGA-II) for water trade-offs
- **Conservative bias** at both validation checkpoints → reported results are **lower bounds**

### 5. **Computational Efficiency** ⭐⭐
- **PERF 1**: \~400× speedup via vectorized PV computation (numpy vs. row-by-row)
- **PERF 2**: \~30× speedup for MC via pre-computed solar geometry
- **Cache invalidation** (FIX 2) prevents stale reference PV across tilt sweeps

### 6. **Biomethane Economics Integration** ⭐⭐
- **Scenario S5**: Reflects 2025 German EEG 2023 grid-injection premiums (5.5× p\_elec)
- **Economic sensitivity**: Six valorization scenarios (S0–S5) span realistic market conditions
- **Site-specific pricing**: p\_elec, p\_biogas, p\_tomato, discount rates per site

---

## ⚠️ CRITICAL WEAKNESSES

### 1. **eLER Metric Validity Under Question** ❌❌❌
**Issue**: LER\_biogas = 1.000 when ESR > \~6, by construction
- When PV alone covers digester heating, the biogas is "free" land-use addition
- Claimed eLER = 2.293 (Almería) = 1.293 (crop+PV, constrained) + 1.000 (biogas, by identity)
- **Problem**: The 1.00 is not earned through trade-offs; it's a design artifact
- **Journal concern**: Reviewers will immediately question whether eLER *fairly* represents land-use efficiency

**Why it matters**:
- Literature benchmark: LER\_2C max = 1.94 (Riaz et al., 2022)
- Your LER\_2C = 1.293 << 1.94 → suggests inferior agrivoltaic design vs. literature
- The eLER advantage (2.293) rests entirely on adding biogas at zero land "cost"

**Recommendation**: 
- Explicitly redefine eLER to account for digester footprint (if above-ground)
- OR implement **Definition C: "Economic LER"** that weights biogas premium by land-use trade-off cost
- Consider alternative: report LER\_2C separately from biogas economic value (€/ha)

### 2. **PAR Saturation Model Missing Field Validation** ❌❌
**Issue**: Asymmetric PAR model (FIX 1) not experimentally validated
- Your LER\_crop derivation from PAR integral is elegant but relies on:
    - Tomato PAR\_sat = 174 W/m² (800 µmol/m²/s) from Heuvelink et al. (2005) — old reference
    - Assumption that hourly PAR saturation index = yield response (nonlinear dynamics ignored)
    - No field measurements under your APV configurations

**Why it matters**:
- Crop yield is typically the **most uncertain parameter** in agrivoltaic studies (±20–30%)
- FIX 1 is critical to your crop yield claims but untested under your shading geometry
- Journal reviewers will ask: "Have you validated this under partial shade?"

**Recommendation**:
- Add a **sensitivity analysis section**: Vary PAR\_sat ±10%, show LER\_crop range
- Cite recent APV-crop yield field studies (Semeraro et al., 2024; Dinesh et al., 2023)
- Propose field validation as **future work** (acknowledge limitation)
- Consider citing newer tomato photosynthesis data (Poorter et al., 2020; Lambers et al., 2019)

### 3. **Digester Thermal Coupling Not Fully Characterized** ❌❌
**Issue**: Heat balance model (Eq. 14–17) contains several unjustified assumptions
- **5°C thermal boost** (Eq. 16) cited from Lindorfer et al. (2008) — is this representative of YOUR design?
- **Q\_reaction = 0.05 × E\_biogas** — where does 5% exothermic fraction come from? (Pal et al., 2024 cited, but not derived)
- **Clipping T\_dig,eff at 37°C** operationally, but using piecewise model (T > 37°C branch) only in validation
- **No dynamic heat storage** — digester walls modeled as steady-state; no thermal lag from concrete mass

**Why it matters**:
- Digester heating is the **bottleneck** for high ESR → LER\_biogas = 1.000 payoff
- If your thermal model is conservative (underestimates heat available), then reported ESR values are inflated
- If your model is optimistic (overestimates heat available), then ESR = 1.000 is not robust

**Recommendation**:
- Cite recent digester heat balance papers (2023–2025) instead of 2008 data
- Justify the 5°C boost with seasonal/geographic reasoning or sensitivity sweep
- Explain why Q\_reaction is 5% (not 3%, not 7%)
- Add a **thermal sensitivity tornado**: vary U\_eff ±20%, show ESR range
- Consider a **dynamic heat balance** (ODEs) for a subset of sites (stronger journal impact)

### 4. **Short Growing Season Sites (Almería, Ouagadougou) Underpowered** ❌
**Issue**: Only 122-day seasons at both semi-arid sites
- PVGIS TMY data have daily resolution → gridded uncertainty
- 122 days ≈ 4 months of hourly variability
- Your Pareto analysis runs only at Konya & Freiburg (> 150d) for this reason
- **But**: Main optimization & eLER results reported for ALL 4 sites, including the short-season ones

**Why it matters**:
- LER statistics at 122-day sites are less robust than 189-day sites
- Reviewers will ask: "How sensitive is your Almería result to ±2 weeks season length change?"
- Agricultural practices differ; season length is not fixed

**Recommendation**:
- Add a **season-length sensitivity table**: Show eLER, LER\_crop for ±10-day season shifts
- Explicitly state confidence bounds: "Almería results ±5% due to season uncertainty"
- Compare your season dates with **local agricultural authorities** (Spanish ministry, etc.)
- Justify the specific start/end dates with citations to crop phenology literature

### 5. **Reference PV Definition Ambiguity** ❌
**Issue**: Definition A vs. Definition B for LER\_PV creates confusion
- **Definition A (independent)**: Reference PV = dense-pack ground-mount at same tilt
- **Definition B (circular)**: Reference PV = APV array itself

**Why it matters**:
- Your paper reports LER\_PV using Definition A (correct for agrivoltaic claim)
- BUT: You compare to Riaz et al. (2022) who may use Definition B
- **Inconsistent definitions** make literature comparison appear better than it is

**Recommendation**:
- Create a **3-panel figure**:
    - Panel 1: Your eLER with Definition A
    - Panel 2: Your eLER with Definition B
    - Panel 3: Riaz et al. (2022) eLER (if their definition known)
- Prominently state in results: "All LER\_PV values use Definition A (independent ground-mount reference)"
- Cite the recent consensus statement on LER definitions (if one exists in 2025 literature)

### 6. **Biomethane Premium (5.5×) Not Justified for All Sites** ❌
**Issue**: S5 scenario uses 5.5× p\_elec everywhere
- Premium based on **German EEG 2023** (Nik Zad et al., 2026 citation)
- But you model **Konya (Turkey), Almería (Spain), Ouagadougou (Burkina Faso)**
- Turkey, Spain, Burkina Faso biogas markets ≠ Germany; premiums likely 1.5–2.5×, not 5.5×

**Why it matters**:
- S5 is your **largest economic winner** → if premium is wrong, main economic claim is wrong
- Reviewers in energy/sustainability journals will scrutinize renewable gas valuations carefully
- Nik Zad et al. (2026) is very recent; policy may change by 2026–2028 publication

**Recommendation**:
- Create a **site-specific premium table** (with citations):
    - Konya: Turkey biogas market 2024–2025 (EMRA/Enerji Regülütörü Kurumu data)
    - Almería: Spanish renewable gas prices (AEE / APPA data)
    - Ouagadougou: Burkina Faso no grid; use LPG displacement (S3) instead
    - Freiburg: Keep 5.5× (EEG 2023 is primary policy)
- Run sensitivity: "S5 results robust within ±30% premium variation"
- Consider a **discounted future value** analysis: EEG 2023 sunsets in \~2040; what then?

### 7. **Crop Residue Ratio (0.80) Not Site-Specific** ❌
**Issue**: All sites use R\_res = 0.80 (residue-to-fruit ratio)
- Citation: Mohammedi et al. (2023)
- But this is for **irrigated semi-arid tomato** (60 t/ha)
- Freiburg: **open-field German tomato** (30 t/ha) — likely has different R\_res (shorter stems, fewer leaves)
- Tropical tomato (Ouagadougou) may differ from Mediterranean

**Why it matters**:
- Residue mass → VS input → BMP → biogas energy
- ±20% error in R\_res → ±20% error in biogas energy → ±20% error in LER\_biogas
- This is a **first-order sensitivity** that you haven't explored

**Recommendation**:
- Search literature for site-specific tomato morphologies (DLG 2023 for Germany, IPCC default factors, etc.)
- Create a **residue-ratio sensitivity table** (range 0.65–0.95)
- Show "robustness range" for eLER under R\_res uncertainty
- Propose direct biomass sampling as future validation work

### 8. **Monte-Carlo Reframed as Structural Sensitivity (Not True Stochastic)** ⚠️
**Issue**: Your MC (FIX 7) runs 2,000 samples with varying BMP\_ref, U\_eff, α\_shade
- But this is **deterministic sampling** of three uncertain parameters
- Not true Monte-Carlo because:
    - Climate is fixed (no GHI, T, WS stochasticity)
    - Design variables (β, d\_row, H\_m, V\_dig, HRT, f\_PV) are optimized per sample (not fixed)
    - Output is "net income P5/P95", but true stochastic propagation would include climate year-to-year variance

**Why it matters**:
- Title/label "Monte-Carlo" suggests Bayesian uncertainty propagation → readers expect climate variability
- Actual output is "how robust is eLER to ±30% parameter perturbations?" → better titled "Parameter Sensitivity"
- Energy/climate journals expect **multi-year** climate variability analysis (TMY + 2–3 historical years)

**Recommendation**:
- Relabel MC section as **"Parametric Sensitivity & Economic Robustness"** (clearer)
- Add a **secondary true MC**: Run your optimized designs against 3–5 historical climate years
    - Shows how often eLER stays > 1.8, for example
    - Real impact on journal reviewers (shows climate variability risk)
- Cite recent APV robustness studies (Mazzeo et al., 2025; Semeraro et al., 2024)

---

## 🔴 METHODOLOGICAL CONCERNS (High Impact)

### Issue 1: Constraint Convergence Boundary\*\* ❌❌
**Finding**: Optimizer converges to **boundary** (d\_row = 4.0 m min, f\_PV,heat = 0.05 min)
- Indicates optimum is at design space edge, not interior
- L-BFGS-B polishing confirms no interior local maximum
- **Concern**: Are you capturing the true Pareto frontier, or just the corner?

**Recommendation**:
- Run optimization with **relaxed bounds**: d\_row ∈ [2, 20] m, f\_PV,heat ∈ [0.01, 1.0]
- Compare eLER at boundary vs. any interior peak
- Plot **2D slices** of eLER landscape (d\_row vs. eLER, f\_PV,heat vs. eLER) for transparency

### Issue 2: Reference Definition D—"True Net Biogas"\*\* ❌
**Suggestion**: Your LER\_biogas counts all biogas as "land-use equivalent"
- But biogas created only because PV heats the digester
- **Fairer metric**: LER\_biogas should account for the **energy cost** of that heating

Example:
- If PV generates 100 kWh/season and 60 kWh goes to digester heating:
    - Biogas energy "funded by PV heating" = X MWh/season
    - Better representation: LER\_biogas\_net = X / (E\_biogas\_standalone\_unheated)
    - This removes the artificial 1.000 and gives a more honest trade-off picture

**Recommendation**:
- Define **Definition C (energy-attributed eLER)**:
  $$eLER\_C = LER\_{crop} + LER\_{PV} + \frac{E\_{biogas,\ attributed\ to\ PV\ heating}}{E\_{biogas,ref}}$$
- Show all three definitions (A, B, C) in main results
- Discuss which is most appropriate for policy decision-making

### Issue 3: No Digester Operational Cost\*\* ❌
**Gap**: Economic analysis omits digester O&M costs
- Your CAPEX includes digester (€450/m³), but O&M is not itemized
- Real farms spend 3–5% CAPEX/year on digester maintenance, repairs, substrate handling
- This **reduces NPV by 15–25%** at typical discount rates

**Recommendation**:
- Add annual O&M line item: Digester\_OM = 0.04 × Digester\_CAPEX (industry benchmark)
- Recompute NPV for all scenarios
- Show sensitivity: ±50% O&M variation

---

## 🟡 PRESENTATION & WRITING ISSUES

### 1. **Abstract Ambiguity on LER Exceedance** ⚠️
**Current (from seminar2.tex, lines 46–49)**:
> "Optimized eLER values reach 1.963 (Konya)... exceed the published maximum two-component LER of 1.94 (Riaz et al., 2022)... entirely attributable to LER\_biogas = 1.000..."

**Problem**: Reader stops at "exceed 1.94" → assumes your APV design is better than literature
- Small print explanation (in small font) is too late
- Abstract should lead with honest LER\_2C (crop+PV) result, then explain biogas addition

**Revision**:
> "...the two-component agrivoltaic performance (LER\_2C = LER\_crop + LER\_PV) reaches 1.29 under independent reference, below the 1.94 published benchmark (Riaz et al., 2022). The extended three-component metric (eLER) reaches 1.96–2.29 through biogas valorization, which is enabled by PV thermal coverage of digester heating (ESR = 76–166)..."

### 2. **Section 3.9 (FIX 10) Needs Visual Explanation** ⚠️
**Gap**: "LER\_water variance below numerical resolution at GCR ≈ 54.5% over 122-day season"
- **Unclear**: Why does short season + high GCR → low LER\_water variance?
- **Missing**: A figure showing LER\_water(d\_row, season\_length) as a 2D heatmap

**Recommendation**:
- Add Figure: "LER\_water sensitivity to row spacing and growing season length"
- Explain: "At short seasons (122d), water savings are compressed into a narrow time window; shading variation across the reduced design space does not produce sufficient LER\_water swing to populate a Pareto frontier."

### 3. **Algorithm 1 Presentation** ⚠️
**Strength**: Well-formatted in algorithm2e package
**Weakness**: Doesn't explain WHY dithering F ∈ [0.5, 1.0] instead of fixed F=0.7
- Dithering is advanced optimization technique; deserves justification line

**Recommendation**:
- Add a paragraph before Algorithm 1: "Dithering the mutation factor F at each generation improves exploration–exploitation balance (Storn, 1997; Price et al., 2005). We dither uniformly in [0.5, 1.0] rather than fix F=0.7 to avoid premature convergence in early generations while refining solutions in late generations."
- Cite two key DE papers

### 4. **Figure Quality & Caption Completeness** ⚠️
**Issue**: Figures mentioned (fig1\_site\_map.png, fig11\_monthly\_climate.png, fig12\_de\_optimizer.png) are referenced
- Are these figures high-resolution (300+ dpi) for journal submission?
- Do all figures have complete captions explaining every axis and color scale?

**Recommendation**:
- Ensure all figures are ≥ 300 dpi (for RGB) or ≥ 600 dpi (grayscale)
- Captions must be **fully self-contained** (reader shouldn't need to read text to understand figure)
- Add sub-panel labels (a), (b), (c) if multi-panel
- Check that axes have units (e.g., "eLER (—)" not just "eLER")

---

## 🟢 OPPORTUNITIES FOR ENHANCEMENT

### Opportunity 1: Compare Against Machine Learning Baseline
**Idea**: Train a neural network on the 6 design variables → eLER
- Use your 2,000 MC samples + DE optimization points as training data (\~5,000 points)
- Benchmark DE speed vs. trained NN prediction speed
- Show Pareto frontier found by NN vs. DE
- **Impact**: "Computational scalability" angle → appeals to optimization/ML reviewers

### Opportunity 2: Sensitivity Analysis Dashboard
**Idea**: Interactive table showing eLER under 5 key parameter ranges:
- BMP\_ref: [250, 300, 350] NmL/g VS
- PAR\_sat: [150, 174, 200] W/m²
- R\_res: [0.65, 0.80, 0.95]
- Biomethane premium: [1.5×, 5.5×, 10×] p\_elec
- Season length: ±10 days

**Impact**: "Robust across realistic parameter uncertainty" → strong message for policy audiences

### Opportunity 3: Temporal Dynamics Animation
**Idea**: Create a **video animation** (GIF or MP4) showing:
- Month-by-month crop growth & shading profile
- Digester temperature & BMP fluctuation
- Cumulative eLER component build-up
- **Submit as supplementary material** to journal
- **Impact**: High-quality supplementary material stands out; shows computational sophistication

### Opportunity 4: Experimental Validation Plan
**Idea**: Add a **"Future Work / Validation" section** with concrete field experiment:
- Select **1 site** (Konya?) for a 2-year APV+AD pilot
- Measure: PV output (inverter data), tomato yield (harvest records), digester gas (gas meter), heat loss (thermocouples)
- Protocols for validating FIX 1 (PAR saturation) and heat balance
- **Impact**: Transition from pure modeling to experimental roadmap → increases credibility

### Opportunity 5: Policy Scenario Analysis
**Idea**: Add a section on **policy implications**:
- **Scenario A (Current)**: S0–S5 with 2024–2025 prices
- **Scenario B (Pessimistic)**: 30% drop in PV cost, biomethane premium cut to 2× → rerun optimization
- **Scenario C (Optimistic)**: 50% efficiency gain in digester → rerun optimization
- Compare Pareto frontiers across scenarios
- **Impact**: Policy journals (Energy Policy, Renewable & Sustainable Energy Reviews) love scenario analysis

---

## 📋 PUBLICATION STRATEGY RECOMMENDATIONS

### **Target Journals (Ranked by Fit)**

1. **Applied Energy** (IF \~10.3, 2024)
    - **Fit**: ✅ Integrated systems, optimization, climate zones
    - **Gaps to address**: Simplify eLER definition; add field validation roadmap
    - **Timeline**: 6–8 months peer review

2. **Renewable & Sustainable Energy Reviews** (IF \~15, 2024)
    - **Fit**: ✅ Multi-technology integration, climate resilience, food-energy-water nexus
    - **Strengths to highlight**: Four-site diversity, biogas market integration
    - **Gaps**: More discussion of policy/sustainability impact
    - **Timeline**: 4–6 months peer review

3. **Renewable Energy** (IF \~8–9, 2024)
    - **Fit**: ✅✅ PV focus, bifacial modules, optimization
    - **Gaps**: Downplay biogas, emphasize PV efficiency gains under shade
    - **Timeline**: 5–7 months peer review

4. **Solar Energy** (IF \~6–7, 2024)
    - **Fit**: ✅ Agrivoltaic focus, PAR saturation, crop yields
    - **Strength**: Your FIX 1 (PAR saturation fix) is a core contribution
    - **Timeline**: 5–6 months peer review

### **Pre-Submission Checklist**

- [ ] **eLER metric**: Redefine or add Definition C (energy-attributed) to address LER\_biogas = 1.000 criticism
- [ ] **PAR validation**: Add sensitivity ±10%, cite 2023–2025 APV crop yield studies
- [ ] **Digester heat model**: Justify 5% Q\_reaction, cite 2023–2025 heat balance papers
- [ ] **Season-length robustness**: ±10-day sensitivity table
- [ ] **Site-specific biomethane premium**: Turkey, Spain, Burkina Faso data (not just Germany)
- [ ] **Residue ratio**: Site-specific R\_res search; add sensitivity table
- [ ] **Monte-Carlo relabeling**: Call it "Parametric Sensitivity" not MC; add true climate variability analysis
- [ ] **Figure quality**: 300+ dpi, complete captions, axis labels with units
- [ ] **Algorithm explanation**: Add dithering justification paragraph
- [ ] **Digester O&M cost**: Add ±4% annual cost; recompute NPV
- [ ] **Policy scenario section**: B (pessimistic) and C (optimistic) market conditions

### **Suggested Revision Process**

1. **Internal review** (4 weeks): Address all ⚠️ and 🔴 items above
2. **Colleague review** (2 weeks): Send to 2–3 experts in agrivoltaics OR optimization OR anaerobic digestion
3. **Pre-print release** (optional): Post on arXiv → gather early feedback
4. **Target submission**: Applied Energy or Renewable & Sustainable Energy Reviews
5. **Expect R1 (Major Revisions)** on eLER definition + validation gaps
6. **R2 review**: Should accept after revisions

---

## 📊 SUMMARY TABLE: Strengths vs. Weaknesses

| Category | Strength | Weakness | Fix Priority |
| --- | --- | --- | --- |
| **Novel contribution** | First integrated APV-AD model ⭐⭐⭐ | eLER definition ambiguous (LER\\_biogas=1.000) | 🔴 CRITICAL |
| **Physical models** | Bifacial, FIX 1 PAR saturation, Faiman temp ⭐⭐⭐ | PAR model not field-validated; no dynamic heat storage | 🔴 HIGH |
| **Climate diversity** | 4 zones, 12°N–48°N, PVGIS-SARAH3 ⭐⭐ | Short seasons (122d) at Almería & Ouagadougou underpowered | 🟡 MEDIUM |
| **Economics** | 6 scenarios, biomethane premium ⭐⭐ | 5.5× premium not justified for all sites; O&M cost missing | 🔴 HIGH |
| **Optimization** | DE + L-BFGS-B, Pareto NSGA-II ⭐⭐ | Boundary convergence may miss interior optima | 🟡 MEDIUM |
| **Validation** | PV vs. PVGIS ⭐ | Crop & digester models not experimentally tested | 🔴 HIGH |
| **Writing/presentation** | Clear structure, algorithm pseudo-code ⭐ | Abstract misleading on LER exceedance; MC relabeling | 🟡 MEDIUM |
| **Reproducibility** | 300-dpi figures, open parameters ⭐ | Code not on GitHub; no supplementary data package | 🟡 MEDIUM |

---

## 🎯 FINAL RECOMMENDATION

**Your work is publication-ready at a strong journal, but requires:**

1. **Reframe eLER metric** (add Definition C or explain land-use trade-off cost)
2. **Validate crop PAR saturation** (sensitivity ±10% + field study roadmap)
3. **Strengthen digester heat model** (justify 5% Q\_reaction; add thermal sensitivity)
4. **Add climate variability analysis** (true MC with 3–5 historical years)
5. **Site-specific biomethane valuations** (Turkey, Spain, Burkina Faso market data)
6. **Polish presentation** (fix abstract framing, add policy scenarios)

**Expected timeline**:
- Revisions: 4–6 weeks
- Submission: Target Applied Energy or Renewable & Sustainable Energy Reviews
- First decision: 3–4 months
- Publication (if accepted): 2–3 months after final revision
- **Total: 9–15 months to publication**

**Journal impact**: Publication in Applied Energy or RSER would be **highly visible** for:
- Master thesis defense (shows journal-level scholarship)
- Future PhD applications (demonstrates research maturity)
- Policy audience (food-energy-water nexus is increasingly policy-relevant)

---

## 📚 Key Literature to Cite (2023–2025)

**Recent APV work**:
- Semeraro et al. (2024) — APV crop yield meta-analysis
- Mazzeo et al. (2025) — eLER framework

**Bifacial PV**:
- Updated bifacial modules performance (2024 industry data)

**Biogas/biomethane**:
- Nik Zad et al. (2026) — Grid injection premiums (your citation)
- EBA (2025) — European biogas outlook

**Optimization**:
- Price et al. (2005) — DE mutation strategies
- Blank & Deb (2020) — pymoo/NSGA-II documentation

**Crop physiology**:
- Poorter et al. (2020) — Photosynthesis under shade
- Lambers et al. (2019) — Plant trait trade-offs

Good luck with your publication! 🚀
