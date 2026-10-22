# Fluid Antenna System (FAS) under Finite Blocklength: BLER Analysis and Neural Network Port Prediction

## Overview

MATLAB implementation of a **Fluid Antenna System (FAS)** study for short-packet (URLLC / 6G) communications. A user device with a single RF chain selects the best of `N` configurable antenna ports over a spatially correlated Rayleigh fading channel.

The project covers:

1. **Statistical analysis** of the FAS channel gain using a *Spatial Block-Correlation (SBC)* model.
2. **Block error rate (BLER)** analysis in the **finite blocklength (FBL)** regime, which is more realistic than outage probability for short packets.
3. **Low-complexity approximations** (Gauss-Laguerre quadrature, Gamma moment matching, high-SNR BLER asymptotic).
4. **Neural network port prediction**, which estimates the best port from a few measured ports instead of scanning all of them.

Course project (Advanced Communication Networks and Systems, *MHTTT*), Faculty of Electronics and Telecommunications, The University of Danang - University of Science and Technology (DUT). Group 11, supervised by Assoc. Prof. Le Thi Phuong Mai.

## Key Features

- **FAS channel model:** Toeplitz spatial correlation matrix (Clarke model, `sinc`) and its block-diagonal SBC approximation (`BlockCorrelation` function in `Main_FBL_FAS.m`).
- **Derived PDF/CDF:** channel-power PDF of the best port, computed with **Gauss-Laguerre quadrature** (order `N_GL = 128`) using non-central chi-square functions (`ncx2cdf`, `ncx2pdf`).
- **Monte Carlo validation:** 2x10^5 correlated Rayleigh samples generated with the *full* Toeplitz matrix (Cholesky), compared against the SBC-based analytical PDF/CDF.
- **Gamma approximation:** shape/scale obtained by matching the mean and variance of the channel power.
- **FBL BLER analysis:** normal approximation with `U = 100` bits and blocklength `M = 400` (`Rc = 0.25`); exact expression vs. high-SNR asymptotic, averaged with exact PDF and Gamma PDF (four combinations), SNR from -30 dB to 30 dB.
- **Conventional baseline:** `L`-antenna MRC system over Rayleigh fading (`L = 1` by default, configurable) for comparison with FAS.
- **Neural network port prediction:** feedforward network (hidden layers 256/128, `trainscg`) predicts the channel power of all 30 ports from 5 measured ports; compared against exhaustive search using best-port accuracy, retained channel gain and outage probability.

## Project Structure

```
MHTTT_project/
├── Main_FBL_FAS.m                              # FAS statistics + FBL BLER (SBC, GL quadrature, Gamma)
├── Main_Conventional_Multiple_Antenna_System.m # Baseline: L-antenna MRC, average BLER
├── Step1_GenerateData.m                        # Correlated channel data for ML (N = 30, 5 observed ports)
├── Step2_TrainModel.m                          # Train feedforward network
├── Step3_TestAndEvaluate.m                     # Compare ML port selection with exhaustive search
├── Sigma.mat                                   # Stored MATLAB data file
└── .gitignore
```

| File | Purpose |
|------|---------|
| `Main_FBL_FAS.m` | Builds the correlation matrix and SBC block sizes, runs Monte Carlo, computes the PDF/CDF (Gauss-Laguerre), Gamma approximation and average BLER vs. SNR. |
| `Main_Conventional_Multiple_Antenna_System.m` | Average BLER vs. SNR for a conventional `L`-antenna MRC receiver (`t_L ~ Gamma(L, sigma^2)`). |
| `Step1_GenerateData.m` | Generates 200,000 correlated channel samples over 30 ports (`W = 1`); input = power at ports `[1, 8, 15, 23, 30]`, target = power at all ports. Saves `FAS_ML_Data.mat`. |
| `Step2_TrainModel.m` | Trains the network on an 80/20 split and saves `FAS_Trained_Net.mat`. |
| `Step3_TestAndEvaluate.m` | Evaluates the predicted best port against exhaustive search. |
| `Sigma.mat` | Stored MATLAB data file included in the repository. |

## Default Simulation Parameters

| Parameter | `Main_FBL_FAS.m` | ML pipeline (`Step1`-`Step3`) |
|-----------|------------------|-------------------------------|
| Ports `N` | 10 | 30 (5 observed) |
| Normalized length `W` | 0.5 | 1 |
| Blocklength `M` / info bits `U` | 400 / 100 | - |
| Channel variance `sigma^2` | 2 | - |
| Block correlation `mu^2` | 0.97 | - |
| SNR range | -30 to 30 dB | - |

The parameters in `Main_FBL_FAS.m` can be edited (for example `N = 50`, `W = 5`) to reproduce the larger-aperture scenarios discussed in the project report.

## Methodology

### 1. Channel model and SBC approximation
Port correlation follows a Toeplitz matrix with `sinc` entries. For the analysis, the full matrix is approximated by a block-diagonal matrix with `B` independent blocks of strongly correlated ports (`mu^2 = 0.97`), so the CDF of the best-port gain factorizes into a product over blocks.

### 2. PDF/CDF and approximations
The exact PDF contains infinite integrals; they are evaluated by Gauss-Laguerre quadrature. A Gamma distribution (moment matching) gives a closed-form alternative.

### 3. Finite blocklength BLER
The instantaneous BLER uses the normal approximation `Q((C - Rc) / sqrt(V / M))`. Average BLER is obtained by integrating over the FAS channel PDF, using either the exact or the Gamma PDF, and either the exact BLER or its high-SNR asymptotic.

### 4. Neural network port selection
The network predicts channel power for all ports from 5 measured ports. The predicted best port is compared with the true best port. The evaluation reports:

- **Best-port selection accuracy**
- **Retained channel gain** relative to the true maximum
- **Outage probability** of the selected port

## Reported Results (from the project report)

These values come from the group's final report and may not match a default run of the scripts; they depend on the parameters above.

- Theoretical PDF/CDF match the Monte Carlo histogram and CDF.
- The four BLER evaluation methods (exact/Gamma PDF x exact/asymptotic BLER) nearly coincide at medium and high SNR.
- With `N = 50`, `W = 5`, FAS with a single RF chain reaches BLER comparable to a conventional system with 20 antennas.
- Neural network port prediction with 5 of 30 ports measured trades roughly 1% of channel gain for a large reduction in scanning cost.

## Requirements

- MATLAB
- Toolboxes providing the functions used: Deep Learning / Neural Network Toolbox (`fitnet`, `train`) and Statistics and Machine Learning Toolbox (`ncx2cdf`, `ncx2pdf`, `gampdf`, `gamcdf`, `qfunc` may require Communications Toolbox depending on version)

## How to Run

Machine learning workflow (run in order, from the project folder):

```matlab
Step1_GenerateData
Step2_TrainModel
Step3_TestAndEvaluate
```

Analytical FAS simulation:

```matlab
Main_FBL_FAS
```

Conventional baseline (set `L` inside the script):

```matlab
Main_Conventional_Multiple_Antenna_System
```

Run the scripts from the project directory so the generated `.mat` files are saved and loaded correctly.

## Technologies

MATLAB, Fluid Antenna Systems, Rayleigh fading and spatial correlation, finite blocklength communications, BLER analysis, Gauss-Laguerre quadrature, Monte Carlo simulation, neural networks (regression).
