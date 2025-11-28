# RGC_Neuroprotection

# RGC_Neuroprotection

## Magnesium Neuroprotection in Retinal Ganglion Cells: A Computational Study

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)](https://www.mathworks.com/products/matlab.html)
# RGC_Neuroprotection

## Magnesium Neuroprotection in Retinal Ganglion Cells: A Computational Study

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b+-blue.svg)](https://www.mathworks.com/products/matlab.html)

This repository contains the MATLAB code for investigating magnesium (Mg²⁺)-mediated neuroprotection in retinal ganglion cells (RGCs) during glutamatergic excitotoxicity, as described in our manuscript submitted to PLOS ONE.

---

## Overview

Retinal ganglion cells are vulnerable to excitotoxic damage in conditions such as glaucoma and retinal ischemia. This computational study systematically investigates how extracellular Mg²⁺ concentration and intervention timing affect neuroprotection through voltage-dependent NMDA receptor block.

### Key Findings

- **Frequency-dependent therapeutic windows**: The effective Mg²⁺ range narrows from 2.0 mM width at physiological frequencies to 0.4 mM at excitotoxic frequencies (80 Hz)
- **Optimal therapeutic range**: 1.6–2.0 mM at 80 Hz satisfies both neuroprotection (Ca²⁺ < 1.0 μM) and function preservation (spike loss ≤20%)
- **Spike loss plateau**: An "optimal zone" where spike reduction plateaus at 20% while calcium reduction continues
- **Critical timing**: Mg²⁺ must be applied within 0.5 seconds of stress onset for meaningful protection

---

## Repository Structure

```
RGC_Neuroprotection/
│
├── README.md                         # This file
├── LICENSE                           # MIT License
│
├── RGC_Mg_Publication_Final.m        # Main analysis script (generates Figures 1-5)
└── RGC_Mg_Supplementary_Final.m      # Supplementary analyses (generates Figures S1-S2)
```

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
Voltage-dependent block following Jahr-Stevens formulation:
```
B(V) = 1 / (1 + η·[Mg²⁺]·exp(-γ·V))
```
where η = 0.28 mM⁻¹ and γ = 0.062 mV⁻¹

---

## Quick Start

### Requirements
- MATLAB R2020b or later
- No additional toolboxes required

### Running the Main Analysis

```matlab
% Generate all main figures (Figures 1-5)
run('RGC_Mg_Publication_Final.m')
```

This script performs:
1. Dose-response simulations across Mg²⁺ concentrations (0.2–2.5 mM)
2. Frequency analysis (10, 30, 60, 80 Hz)
3. Therapeutic window identification
4. Intervention timing analysis
5. Automatic figure generation and saving

### Running Supplementary Analyses

```matlab
% Generate supplementary figures (Figures S1-S2)
run('RGC_Mg_Supplementary_Final.m')
```

This script performs:
1. High-resolution analysis at 80 Hz (0.1 mM steps from 1.0–2.5 mM)
2. Detailed results tables

---

## Output Files

Running the scripts generates the following figures:

### Main Figures
| Figure | Filename | Description |
|--------|----------|-------------|
| Figure 1 | `Figure1_ModelOverview.png/.pdf` | Model schematic and example traces |
| Figure 2 | `Figure2_DoseResponse.png/.pdf` | Dose-response curves |
| Figure 3 | `Figure3_TherapeuticWindows.png/.pdf` | Frequency-dependent therapeutic windows |
| Figure 4 | `Figure4_InterventionTiming.png/.pdf` | Timing analysis |
| Figure 5 | `Figure5_Mechanism.png/.pdf` | Mechanistic demonstration |

### Supplementary Figures
| Figure | Filename | Description |
|--------|----------|-------------|
| Figure S1 | `FigureS1_HighResolution.png/.pdf` | High-resolution analysis at 80 Hz |
| Figure S2 | `FigureS2_ResultsTable.png/.pdf` | Detailed numerical results |

---

## Key Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `Mg_levels` | 0.2–2.5 | mM | Extracellular Mg²⁺ concentrations |
| `frequencies` | 10, 30, 60, 80 | Hz | Stimulation frequencies |
| `Ca_toxic` | 1.0 | μM | Calcium toxicity threshold |
| `max_spike_loss` | 20 | % | Maximum acceptable spike reduction |
| `Tmax` | 3000 | ms | Simulation duration |
| `dt` | 0.02 | ms | Integration time step |
| `transient` | 500 | ms | Initial transient to discard |

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Borjkhani2025,
  title={Magnesium neuroprotection in retinal ganglion cells: A computational 
         study of frequency-dependent therapeutic windows and intervention timing},
  author={Borjkhani, Mehdi and Borjkhani, Hadi and Sharif, Morteza A.},
  journal={PLOS ONE},
  year={2025},
  note={Under review}
}
```

---

## Related Publications

- Borjkhani M, Bahrami F, Janahmadi M. (2018). Computational modeling of opioid-induced synaptic plasticity in hippocampus. *PLoS ONE*, 13(3):e0193410.
- Borjkhani M, Bahrami F, Janahmadi M. (2018). Assessing the effects of opioids on pathological memory by a computational model. *Basic and Clinical Neuroscience*, 9(4):275-288.
- Borjkhani M, Borjkhani H, Sharif MA. (2022). Investigating the cocaine-induced reduction of potassium current on the generation of action potentials using a computational model. *Basic and Clinical Neuroscience*, 13(6):765-774.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Mehdi Borjkhani, PhD**  
International Centre for Translational Eye Research (ICTER)  
Institute of Physical Chemistry, Polish Academy of Sciences  
Warsaw, Poland

📧 Email: mborjkhani@ichf.edu.pl

---

## Acknowledgments

We gratefully acknowledge Professor Maciej Wojtkowski, Head of ICTER, for his support and for providing an inspiring research environment.

---

*Last updated: November 2025*
