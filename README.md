# RGC_Neuroprotection

## Magnesium Neuroprotection in Retinal Ganglion Cells: A Computational Study

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![Status](https://img.shields.io/badge/Status-Accepted-brightgreen.svg)]()
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19846309.svg)](https://doi.org/10.5281/zenodo.19846309)


This repository contains the MATLAB code for investigating magnesium (Mg²⁺)-mediated neuroprotection in retinal ganglion cells (RGCs) during glutamatergic excitotoxicity, as described in our manuscript accepted for publication in PLOS ONE.

---
## Overview

Retinal ganglion cells are vulnerable to excitotoxic damage in conditions such as glaucoma and retinal ischemia. This computational study systematically investigates how extracellular Mg²⁺ concentration and intervention timing affect neuroprotection through voltage-dependent NMDA receptor block.

### Key Findings

- **Frequency-dependent therapeutic windows**: The effective Mg²⁺ range narrows from 2.0 mM width at physiological frequencies (10 Hz) to 0.4 mM at excitotoxic frequencies (80 Hz), closing entirely at 90–100 Hz
- **Optimal therapeutic range**: 1.6–2.0 mM at 80 Hz satisfies both neuroprotection (Ca²⁺ < 1.0 μM) and function preservation (spike loss ≤ 20%)
- **Spike loss plateau**: An optimal zone where spike reduction plateaus at 20% while calcium reduction continues to improve (1.4–2.0 mM)
- **Critical timing**: Protection efficacy is maximal with pre-treatment or immediate intervention (100%), declines steeply with delay (> 50% protection requires intervention within 0.2 s in the model protocol), reflecting relative phase sensitivity rather than a literal clinical deadline
- **State-based threshold**: ≥ 50% protection requires intervention before ~35% of peak Ca²⁺ accumulation, generalizing the timing result across timescales
- **Numerical validation**: Forward Euler at dt = 0.02 ms agrees with RK4 to < 0.2% in peak Ca²⁺

---

## Repository Structure

```
RGC_Neuroprotection/
│
├── README.md                         # This file
├── LICENSE                           # MIT License
│
└── RGC_Mg_Revised.m                  # Main analysis script — generates all
                                      # main figures (1–5) and supplementary
                                      # figures (S1–S5), including revised
                                      # analyses addressing reviewer comments
```

> **Note**: This revised version consolidates the main and supplementary analyses into a single script for reproducibility, and adds four new analyses (extended frequency range, Ca²⁺ threshold sensitivity, numerical method validation, expanded parameter table) requested during peer review.

---

## Model Description

The model implements a biophysically detailed, conductance-based RGC with:

### Intrinsic Conductances

| Channel | Conductance | Description |
|---------|-------------|-------------|
| Na⁺ | g_Na = 120 mS/cm² | Fast sodium (spike generation) |
| K⁺ (DR) | g_Kdr = 36 mS/cm² | Delayed rectifier potassium |
| K⁺ (A) | g_KA = 8 mS/cm² | A-type potassium |
| Ca²⁺ (L) | g_CaL = 0.3 mS/cm² | L-type calcium |
| K⁺ (Ca) | g_KCa = 0.3 mS/cm² | Calcium-activated potassium |
| Leak | g_L = 0.35 mS/cm² | Leak conductance |

### Synaptic Receptors

| Receptor | Parameters | Description |
|----------|------------|-------------|
| AMPA | τ_rise = 0.3 ms, τ_decay = 3 ms | Fast glutamatergic transmission |
| NMDA | τ_rise = 5 ms, τ_decay = 80 ms | Slow, Mg²⁺-sensitive transmission |

### Mg²⁺ Block
Voltage-dependent block following the Jahr–Stevens formulation:
```
B(V) = 1 / (1 + η · [Mg²⁺] · exp(−γ · V))
```
where η = 0.28 mM⁻¹ and γ = 0.062 mV⁻¹

### Calcium Dynamics
Intracellular calcium evolves via:
```
d[Ca²⁺]ᵢ/dt = −k_NMDA · f_Ca · I_NMDA − k_CaL · I_CaL − ([Ca²⁺]ᵢ − [Ca²⁺]_rest) / τ_Ca
```
with τ_Ca = 200 ms, f_Ca = 0.15, [Ca²⁺]_rest = 0.05 μM.

---

## Quick Start

### Requirements
- MATLAB R2020b or later
- No additional toolboxes required

### Running the Analysis

```matlab
% Generate all main and supplementary figures
run('RGC_Mg_Revised.m')
```

The script runs five sequential analyses:

1. **Main dose-response simulations** — Mg²⁺ × frequency grid (8 × 6 = 48 conditions)
2. **Intervention timing analysis** — 12 delay conditions at 80 Hz
3. **High-resolution analysis at 80 Hz** — 0.1 mM steps from 1.0 to 2.5 mM
4. **Ca²⁺ toxicity threshold sensitivity** — 17 threshold values (0.6–1.4 μM in 0.05 μM steps)
5. **Numerical method validation** — Forward Euler (4 step sizes) vs. RK4

Estimated runtime: ~5–15 minutes depending on hardware.

---

## Output Files

Running the script generates the following figures (saved as `.tif` at 600 DPI, `.pdf` vector, and `.png` preview):

### Main Figures

| Figure | Filename | Description |
|--------|----------|-------------|
| Figure 1 | `Figure1_ModelOverview` | Model schematic, Mg²⁺ block curve, voltage and Ca²⁺ traces |
| Figure 2 | `Figure2_DoseResponse` | Dose-response curves across 10–100 Hz (functional cost & neuroprotection) |
| Figure 3 | `Figure3_TherapeuticWindows` | Trade-off space at 80 Hz; window width vs. frequency |
| Figure 4 | `Figure4_InterventionTiming` | Ca²⁺ dynamics and protection efficacy vs. intervention delay |
| Figure 5 | `Figure5_Mechanism` | Mechanistic cascade: B(V), NMDA/AMPA currents, Ca²⁺ reduction |

### Supplementary Figures

| Figure | Filename | Description |
|--------|----------|-------------|
| Figure S1 | `FigureS1_HighResolution` | High-resolution dose-response and trade-off space at 80 Hz |
| Figure S2 | `FigureS2_ResultsTable` | Detailed numerical results table at 80 Hz (0.1 mM resolution) |
| Figure S3 | `FigureS3_ThresholdSensitivity` | Ca²⁺ toxicity threshold sensitivity analysis (0.6–1.4 μM) |
| Figure S4 | `FigureS4_NumericalValidation` | Euler convergence and RK4 agreement |
| Figure S5 | `FigureS5_ExpandedParameters` | Complete model parameter table with references |

---

## Key Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `Mg_levels` | 0.2, 0.5, 1.0, 1.4, 1.6, 1.8, 2.0, 2.5 | mM | Coarse Mg²⁺ sweep |
| `Mg_highres` | 1.0 : 0.1 : 2.5 | mM | High-resolution Mg²⁺ sweep |
| `frequencies` | 10, 30, 60, 80, 90, 100 | Hz | Stimulation frequencies |
| `Ca_toxic` | 1.0 | μM | Default Ca²⁺ toxicity threshold |
| `Ca_thresholds` | 0.6 : 0.05 : 1.4 | μM | Sensitivity analysis range |
| `max_spike_loss` | 20 | % | Maximum acceptable spike reduction |
| `Tmax` | 3000 | ms | Simulation duration (main) |
| `dt` | 0.02 | ms | Integration time step |
| `transient` | 500 | ms | Initial transient discarded from analysis |

### Simulation Protocol Details

- **Spike counting window**: T_count = 2500 ms (t = 500–3000 ms)
- **At 80 Hz**: N_pulses = 200 expected pulses in the analysis window
- **Spike loss**: 100 × (1 − N_spikes / N_pulses) — relative to expected pulses, not baseline Mg²⁺
- **Intervention timing**: Stress window t = 0.5–4.5 s (4 s duration), total simulation 6 s

---

## Helper Functions

The script contains four internal simulation functions:

| Function | Description |
|----------|-------------|
| `run_simulation_euler` | Forward Euler integration (main workhorse) |
| `run_simulation_rk4` | 4th-order Runge-Kutta (used for numerical validation only) |
| `run_simulation_detailed` | Euler with full current and Mg²⁺ block recording (Figure 5/6) |
| `derivatives` | State-derivative function called by the RK4 integrator |
| `save_figure` | Saves each figure as TIFF (600 DPI), PDF (vector), and PNG (300 DPI) |

---

## Therapeutic Window Summary

| Frequency | Optimal [Mg²⁺] Range | Window Width |
|-----------|----------------------|--------------|
| 10 Hz | 0.5–2.5 mM | 2.0 mM |
| 30 Hz | 1.0–2.5 mM | 1.5 mM |
| 60 Hz | 1.4–2.5 mM | 1.1 mM |
| 80 Hz | 1.6–2.0 mM | 0.4 mM |
| 90 Hz | None | — |
| 100 Hz | None | — |

---

---

## Data Availability

All simulation code, analysis scripts, and minimal data required to reproduce the manuscript figures and results are included in this repository.

Permanent Zenodo archive: https://doi.org/10.5281/zenodo.19846309

Published article DOI: https://doi.org/10.1371/journal.pone.0348068

## Citation

If you use this code in your research, please cite:



```bibtex
@article{Borjkhani2026,
  title   = {Magnesium neuroprotection in retinal ganglion cells: A computational
             study of frequency-dependent therapeutic windows and intervention timing},
  author  = {Borjkhani, Mehdi and Borjkhani, Hadi and Sharif, Morteza A.},
  journal = {PLOS ONE},
  year    = {2026},
  note    = {Accepted}
}
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Mehdi Borjkhani**  
International Centre for Translational Eye Research (ICTER)  
Institute of Physical Chemistry, Polish Academy of Sciences  
Warsaw, Poland

📧 mborjkhani@ichf.edu.pl  
🔗 [GitHub](https://github.com/borjkhani)

---

*Last updated: April 2026*
