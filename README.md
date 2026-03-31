# An Uncertainty Principle for Neural Coding

MATLAB code for reproducing the simulations in:

> Grgurich, R. & Blair, H.T. (2020). An uncertainty principle for neural coding: Conjugate representations of position and velocity are mapped onto firing rates and co-firing rates of neural spike trains. *Hippocampus*, 30(4), 396-421. [doi:10.1002/hipo.23197](https://doi.org/10.1002/hipo.23197)

## Overview

This paper shows that neural populations can simultaneously embed two codes within their spike trains: a *firing rate code* (R) conveyed by within-cell spike intervals, and a *co-firing rate code* (R-hat) conveyed by between-cell spike intervals. These two codes behave as conjugates of one another, obeying an analog of the uncertainty principle from physics.

Two biologically inspired decoding methods are introduced:
- **Sigma decoding**: recovers information from firing rates (within-cell temporal integration)
- **Sigma-chi decoding**: recovers information from co-firing rates (between-cell temporal integration)

Simulations demonstrate that when firing rates encode position, co-firing rates encode velocity, and vice versa.

## Setup

Before running any simulation, add the project to your MATLAB path:

```matlab
run('setup.m')
```

Then `cd` into any simulation directory and run the script.

## Project Structure

```
uncertainty-principle-for-neural-coding/
├── setup.m                          % Run first: adds all paths
├── data/                            % Behavioral data files
├── lib/
│   ├── spike_generation/            % Spike train simulators
│   ├── conductance/                 % Firing rate & co-firing rate computation
│   ├── decoders/                    % Sigma and sigma-chi decoder training/testing
│   ├── tuning/                      % Tuning curve generators
│   ├── preprocessing/               % Data interpolation
│   └── external/CircStat/           % Circular statistics toolbox (Berens, 2009)
├── simulations/
│   ├── fig1_6_hd_cells/             % Head direction cell simulations
│   ├── fig7_grid_cells/             % Grid cell simulations
│   ├── fig8_11_theta_cells/         % Theta-modulated cell simulations
│   ├── fig12_speed_vs_grid/         % Speed cell vs grid cell comparison
│   └── parameter_sweeps/            % Parameter exploration scripts
└── paper/                           % Published PDF
```

## Code-to-Figure Mapping

### Figures 1-6: Head Direction Cells

These figures demonstrate that HD cell firing rates encode angular head position (recovered by sigma decoding), while co-firing rates encode angular velocity (recovered by sigma-chi decoding).

| Figure | Description | Simulation Script | Key Library Functions |
|--------|-------------|-------------------|----------------------|
| Fig 1 | Sigma decoder recovers HD from firing rates | `fig1_6_hd_cells/HDdata_Decoder_STD.m` (N=12) or `HDdata_Decoder.m` (N=32) | `vonMisesFiringProb` (Eq. 3), `Generate_Conductance` (Eq. 5-6), `Universal_Train` (Eq. 7-8) |
| Fig 2 | Accuracy-latency tradeoff (MSE vs. tau) | `Universal_Train` generates these curves internally (figure 101) | Sweep over time shifts -250:25:250 ms |
| Fig 3 | Chi rates for HD cell pairs | `Generate_Between` computes chi rates (Eq. 9) | Visualization of individual chi rates chi_{i,j} |
| Fig 4 | Chi rate matrix and sigma-chi tuning | `Generate_Between` (Eq. 10), `Universal_Test` lines 14-22 | Co-firing rate tuning curves by velocity |
| Fig 5 | Sigma-chi decoder recovers angular velocity | `fig1_6_hd_cells/HDdata_Decoder.m` | `Generate_Between` (Eq. 9-10), `Universal_Train` (Eq. 12) |
| Fig 6 | Conjugate representations (MSE vs. log tau) | `fig1_6_hd_cells/HDdata_Decoder.m` with tau sweep | Three decoders: sigma-only, sigma-chi-only, combined |

### Figure 7: Grid Cells

Grid cell firing rates encode periodic spatial position; sigma-chi units behave like speed cells.

| Figure | Description | Simulation Script | Key Library Functions |
|--------|-------------|-------------------|----------------------|
| Fig 7a-d | Position and speed decoding from grid cells | `fig7_grid_cells/Grid2Spd_Decoder_STD.m` | `vonMisesFiringProb`, `Universal_Train/Test` |
| Fig 7e | Speed tuning of sigma-chi units | `fig7_grid_cells/Universal_Decoder_grid_velfromHD.m` | Uses real speed data from `spd_data.mat` |
| Fig 7f | Latency of speed tuning | `Universal_Train` latency sweep | |

### Figures 8-11: Theta-Modulated Cells

These figures show the converse case: when theta cells encode position via oscillatory phase (co-firing rates), velocity is encoded in firing rates.

| Figure | Description | Simulation Script | Key Library Functions |
|--------|-------------|-------------------|----------------------|
| Fig 8 | Sigma-chi recovers position from theta co-firing | `fig8_11_theta_cells/Universal_Decoder_theta.m` | `PlacePhaseBurstGenerator` (Eq. 13-18), `Generate_Between` |
| Fig 9 | Sigma recovers velocity from theta firing rates | `fig8_11_theta_cells/Universal_Decoder_theta.m` | Theta cell firing rates are modulated by velocity (Eq. 17) |
| Fig 10 | Complementary ring pairs segregate channels | `fig8_11_theta_cells/Universal_Decoder_grid.m` | `PlacePhaseBurstGenerator` with +/- phase slopes |
| Fig 11 | Non-complementary ring oscillators | `fig8_11_theta_cells/Universal_Decoder_grid.m` (variant) | All positive phase slopes |

### Figure 12: Speed Cells vs. Grid Cells

Demonstrates the uncertainty tradeoff: adding position information to the firing rate channel reduces velocity information in the co-firing rate channel.

| Figure | Description | Simulation Script | Key Library Functions |
|--------|-------------|-------------------|----------------------|
| Fig 12 | Uncertainty tradeoff between speed and grid cells | `fig12_speed_vs_grid/CoFiring_Decoder_STD.m` | Joint position/velocity decoding from both channels |

### Parameter Sweeps

| Script | Purpose |
|--------|---------|
| `parameter_sweeps/Universal_Decoder.m` | Compares CV0 (regular) vs CV1 (Poisson) spike trains |
| `parameter_sweeps/HDdata_Decoder_grid.m` | Grid search over peak rate and inhibition parameters |
| `parameter_sweeps/Universal_Decoder_widthloop.m` | Sweeps tuning curve width (sigma) |
| `parameter_sweeps/Universal_Decoder_amploop.m` | Sweeps firing rate amplitude |

## Core Pipeline

All simulations follow the same pipeline:

```
Behavioral Data (head direction, position, or speed)
    │
    ▼
Spike Generation ──────────────────────────────────────────────────────┐
    │ vonMisesFiringProb (Eq. 3) — Poisson HD/grid cells              │
    │ CV0spikeGenerator — regular (CV=0) HD cells                     │
    │ PlacePhaseBurstGenerator (Eq. 13-18) — theta cells              │
    │                                                                  │
    ▼                                                                  ▼
Generate_Conductance (Eq. 5-6)              Generate_Between (Eq. 9-10)
    │ Within-cell firing rates (R)              │ Between-cell co-firing rates (R-hat)
    │ Exponential decay kernel                  │ Pairwise chi rates pooled by separation
    │                                           │
    ▼                                           ▼
Universal_Train / Universal_Test            Universal_Train / Universal_Test
    │ Sigma decoder (Eq. 7-8)                   │ Sigma-chi decoder (Eq. 12)
    │ Pseudoinverse weight fitting              │ Pseudoinverse weight fitting
    │                                           │
    ▼                                           ▼
Position (q) or Velocity (q-dot)          Velocity (q-dot) or Position (q)
depending on what the firing               (conjugate of what firing
rates encode                               rates encode)
```

## Key Equations (from paper)

| Equation | Description | Implementation |
|----------|-------------|----------------|
| Eq. 3 | Von Mises spike probability | `vonMisesFiringProb.m` |
| Eq. 5-6 | Exponential decay kernel (spike → firing rate) | `Generate_Conductance.m` |
| Eq. 7 | Weighted sum of firing rates (sigma decoder) | `Universal_Train.m` |
| Eq. 8 | Arctangent recovery of angle | `Universal_Train.m` (via `cart2pol`) |
| Eq. 9 | Chi rate between spike train pairs | `Generate_Between.m` |
| Eq. 10 | Co-firing rate (sum of chi rates by separation) | `Generate_Between.m` |
| Eq. 12 | Weighted sum of co-firing rates (sigma-chi decoder) | `Universal_Train.m` |
| Eq. 13-18 | Theta cell phase coding model | `PlacePhaseBurstGenerator.m` |
| Eq. 19 | Between-channel uncertainty principle | Demonstrated by simulation results |

## Data Files

| File | Description | Size |
|------|-------------|------|
| `HeadDirectionData.mat` | Real HD tracking data from a behaving rat | 104 KB |
| `SeedSet_10x120.mat` | Random seeds for reproducibility | 4 KB |
| `Spd_spikeMatrix.mat` | Pre-generated grid cell spike matrix | 21 MB |
| `spd_data.mat` | Real speed and position data from circular track | 3.4 MB |
| `noThetaSpks.mat` | Co-firing rate data without theta modulation | 21 MB |

## Dependencies

- MATLAB (developed and tested with R2018b+)
- [Circular Statistics Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics) (included in `lib/external/CircStat/`)

## Citation

```bibtex
@article{grgurich2020uncertainty,
  title={An uncertainty principle for neural coding: Conjugate representations of position and velocity are mapped onto firing rates and co-firing rates of neural spike trains},
  author={Grgurich, Ryan and Blair, Hugh T},
  journal={Hippocampus},
  volume={30},
  number={4},
  pages={396--421},
  year={2020},
  publisher={Wiley},
  doi={10.1002/hipo.23197}
}
```

## License

This code is provided for academic and research purposes.
