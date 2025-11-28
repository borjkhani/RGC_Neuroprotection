# RGC_Neuroprotection

# RGC_Neuroprotection

## Magnesium Neuroprotection in Retinal Ganglion Cells: A Computational Study

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)](https://www.mathworks.com/products/matlab.html)

This repository contains the MATLAB implementation of a conductance-based retinal ganglion cell (RGC) model for investigating magnesium (Mg²⁺)-mediated neuroprotection during glutamatergic excitotoxicity.

## Overview

Retinal ganglion cells are vulnerable to excitotoxic damage in conditions such as glaucoma and retinal ischemia. This computational study systematically investigates how extracellular Mg²⁺ concentration and intervention timing affect neuroprotection through voltage-dependent NMDA receptor block.

### Key Features

- **Biophysically detailed RGC model** with Hodgkin-Huxley dynamics
- **AMPA and NMDA receptor-mediated synaptic transmission**
- **Voltage-dependent Mg²⁺ block** following the Jahr-Stevens formulation
- **Intracellular calcium dynamics** for assessing excitotoxicity
- **Systematic parameter sweeps** for therapeutic window identification

## Model Components

| Component | Description |
|-----------|-------------|
| **Intrinsic conductances** | Na⁺, K⁺ (delayed rectifier), K⁺ (A-type), Ca²⁺ (L-type), K⁺ (Ca²⁺-activated), Leak |
| **Synaptic receptors** | AMPA (τ = 3 ms), NMDA (τ = 80 ms) |
| **Mg²⁺ block** | Voltage-dependent, Jahr-Stevens formulation (η = 0.28 mM⁻¹, γ = 0.062 mV⁻¹) |
| **Calcium dynamics** | Single-compartment model with buffering and extrusion |

## Repository Structure

```
RGC_Neuroprotection/
│
├── README.md                    # This file
├── LICENSE                      # MIT License
│
├── model/
│   ├── RGC_model.m              # Main RGC model implementation
│   ├── HH_kinetics.m            # Hodgkin-Huxley gating kinetics
│   ├── synaptic_currents.m      # AMPA/NMDA receptor currents
│   └── calcium_dynamics.m       # Intracellular Ca²⁺ model
│
├── simulations/
│   ├── dose_response.m          # Mg²⁺ concentration sweeps
│   ├── frequency_analysis.m     # Stimulation frequency effects
│   ├── therapeutic_window.m     # Optimal range identification
│   └── intervention_timing.m    # Temporal constraint analysis
│
├── analysis/
│   ├── spike_detection.m        # Action potential detection
│   ├── calcium_metrics.m        # Peak Ca²⁺ and reduction calculations
│   └── protection_efficacy.m    # Intervention efficacy quantification
│
├── figures/
│   ├── Fig1_ModelOverview.m     # Model schematic and basic dynamics
│   ├── Fig2_DoseResponse.m      # Dose-response curves
│   ├── Fig3_TherapeuticWindows.m # Frequency-dependent windows
│   ├── Fig4_InterventionTiming.m # Timing analysis
│   ├── Fig5_Mechanism.m         # Mechanistic demonstration
│   └── FigS1_HighResolution.m   # Supplementary high-resolution analysis
│
└── results/
    └── (generated simulation outputs)
```

## Quick Start

### Requirements

- MATLAB R2020b or later (tested on R2023b)
- No additional toolboxes required

### Running the Simulations

1. Clone the repository:
   ```bash
   git clone https://github.com/borjkhani/RGC_Neuroprotection.git
   cd RGC_Neuroprotection
   ```

2. Open MATLAB and navigate to the repository folder

3. Run the main simulation:
   ```matlab
   % Basic simulation with default parameters
   run('simulations/dose_response.m')
   
   % Generate all figures
   run('figures/Fig1_ModelOverview.m')
   run('figures/Fig2_DoseResponse.m')
   % ... etc.
   ```

### Example Usage

```matlab
% Set simulation parameters
Mg_conc = 1.8;        % mM - extracellular Mg²⁺
freq = 80;            % Hz - stimulation frequency
duration = 3000;      % ms - simulation duration

% Run single simulation
[V, Ca, spikes] = RGC_model(Mg_conc, freq, duration);

% Analyze results
peak_Ca = max(Ca);
spike_count = length(spikes);
```

## Key Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `Mg_concentrations` | 0.2 – 2.5 | mM | Extracellular Mg²⁺ range |
| `frequencies` | 10, 30, 60, 80 | Hz | Stimulation frequencies |
| `Ca_threshold` | 1.0 | μM | Toxicity threshold |
| `spike_loss_limit` | 20 | % | Maximum acceptable spike reduction |
| `dt` | 0.02 | ms | Integration time step |

## Main Findings

1. **Therapeutic Window at 80 Hz**: 1.6 – 2.0 mM Mg²⁺
2. **Spike Loss Plateau**: 20% reduction maintained across optimal range
3. **Calcium Reduction**: Up to 82% at therapeutic concentrations
4. **Critical Timing**: Intervention must occur within 0.5 s of stress onset

## Reproducing Manuscript Figures

Each figure from the manuscript can be regenerated:

```matlab
% Figure 1: Model overview and Mg²⁺-dependent dynamics
run('figures/Fig1_ModelOverview.m')

% Figure 2: Dose-response analysis
run('figures/Fig2_DoseResponse.m')

% Figure 3: Therapeutic windows
run('figures/Fig3_TherapeuticWindows.m')

% Figure 4: Intervention timing
run('figures/Fig4_InterventionTiming.m')

% Figure 5: Mechanistic demonstration
run('figures/Fig5_Mechanism.m')

% Supplementary Figure S1: High-resolution analysis
run('figures/FigS1_HighResolution.m')
```

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Borjkhani2025,
  title={Magnesium neuroprotection in retinal ganglion cells: A computational study of frequency-dependent therapeutic windows and intervention timing},
  author={Borjkhani, Mehdi and Borjkhani, Hadi and Sharif, Morteza A.},
  journal={PLOS ONE},
  year={2025},
  note={Manuscript submitted}
}
```

## Related Publications

- Borjkhani M, Bahrami F, Janahmadi M. (2018). Computational modeling of opioid-induced synaptic plasticity in hippocampus. *PLoS ONE*, 13(3):e0193410.
- Borjkhani M, Bahrami F, Janahmadi M. (2018). Assessing the effects of opioids on pathological memory by a computational model. *Basic and Clinical Neuroscience*, 9(4):275-288.
- Borjkhani M, Borjkhani H, Sharif MA. (2022). Investigating the cocaine-induced reduction of potassium current on the generation of action potentials using a computational model. *Basic and Clinical Neuroscience*, 13(6):765-774.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Mehdi Borjkhani**  
International Centre for Translational Eye Research (ICTER)  
Institute of Physical Chemistry, Polish Academy of Sciences  
Warsaw, Poland

📧 Email: [mborjkhani@ichf.edu.pl](mailto:mborjkhani@ichf.edu.pl)

## Acknowledgments

We gratefully acknowledge Professor Maciej Wojtkowski, Head of ICTER, for his support and for providing an inspiring research environment.

---

*Last updated: November 2025*
