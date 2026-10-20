# Fluid Antenna System (FAS) Simulation and Neural Network Channel Prediction

## Overview

This project explores a **Fluid Antenna System (FAS)** for wireless communications using MATLAB. It combines statistical channel modeling, Monte Carlo simulation, finite blocklength performance analysis, and neural network-based channel prediction to investigate antenna selection and communication reliability.

The project also compares a conventional single-antenna system with FAS and evaluates whether machine learning can identify a high-performing antenna port using only partial channel observations.

## Key Features

- **FAS Channel Modeling:** Construct spatial correlation matrices and simulate correlated Rayleigh fading channels.
- **Statistical Analysis:** Estimate channel power probability density functions (PDFs) and cumulative distribution functions (CDFs).
- **Monte Carlo Simulation:** Generate channel samples and examine the distribution of the strongest antenna-port gain.
- **Gamma Approximation:** Compare the estimated channel-power distribution with a Gamma-distribution approximation.
- **Neural Network Prediction:** Train a feedforward neural network to predict channel power across antenna ports from partial observations.
- **Antenna Port Selection:** Compare the ML-selected port with exhaustive search over all ports.
- **Communication Performance Evaluation:** Analyze block error rate (BLER) against signal-to-noise ratio (SNR), including exact and high-SNR asymptotic approximations.
- **Performance Metrics:** Evaluate best-port selection accuracy, retained channel gain, and outage probability.

## Project Structure

```text
MHTTT/
├── Main_Conventional_Multiple_Antenna_System.m
├── Main_FBL_FAS.m
├── Step1_GenerateData.m
├── Step2_TrainModel.m
├── Step3_TestAndEvaluate.m
└── Sigma.mat
```

| File | Purpose |
|---|---|
| `Main_Conventional_Multiple_Antenna_System.m` | Simulates average BLER for a conventional antenna system over Rayleigh fading. |
| `Main_FBL_FAS.m` | Analyzes FAS channel statistics and evaluates BLER performance using different distribution and approximation methods. |
| `Step1_GenerateData.m` | Generates correlated channel samples and partial-observation training data. |
| `Step2_TrainModel.m` | Trains a feedforward neural network for channel-power prediction. |
| `Step3_TestAndEvaluate.m` | Evaluates predicted antenna-port selection against exhaustive search. |
| `Sigma.mat` | Stored MATLAB data file included in the repository. |

## Methodology

### 1. Channel Data Generation

The data-generation script simulates correlated Rayleigh fading across **30 antenna ports**. It selects five ports as partial observations and uses the channel power across all ports as the prediction target.

The generated dataset is saved as `FAS_ML_Data.mat`.

### 2. Neural Network Training

A feedforward neural network is trained using MATLAB's Neural Network Toolbox. The dataset is divided into 80% training data and 20% testing data.

The model uses two hidden layers with 256 and 128 neurons and the scaled conjugate gradient (`trainscg`) training algorithm.

The trained model and test data are saved to `FAS_Trained_Net.mat`.

### 3. Antenna Port Selection

The trained network predicts channel power across all ports based on observations from five selected ports. The predicted best port is then compared with the true best port found through exhaustive search.

The evaluation reports:

- **Best-port selection accuracy:** Percentage of test samples where the predicted port matches the true best port.
- **Retained channel gain:** Average channel gain achieved by the selected port relative to the true maximum.
- **Outage probability:** Probability that the selected channel gain falls below a predefined threshold.

### 4. BLER Performance Analysis

The FAS simulation evaluates BLER over an SNR range from -30 dB to 30 dB. It compares exact and high-SNR asymptotic BLER expressions, as well as channel statistics obtained from the derived PDF and Gamma approximation.

A conventional antenna system is also simulated to provide a reference for performance comparison.

## Requirements

- MATLAB
- MATLAB toolboxes providing the functions used by the scripts, including Neural Network Toolbox and relevant statistical functions

Functions used in the simulations include `fitnet`, `train`, `qfunc`, `ncx2cdf`, `ncx2pdf`, `gampdf`, and `gamcdf`. Availability depends on the installed MATLAB toolboxes and version.

## How to Run

Run the machine learning workflow in the following order from the MATLAB Command Window:

```matlab
Step1_GenerateData
Step2_TrainModel
Step3_TestAndEvaluate
```

Run the analytical and statistical FAS simulation separately:

```matlab
Main_FBL_FAS
```

Run the conventional antenna-system simulation separately:

```matlab
Main_Conventional_Multiple_Antenna_System
```

Run the scripts from the project directory so that generated `.mat` files are saved and loaded from the expected location.

## Technologies

- MATLAB
- Wireless Communications
- Fluid Antenna Systems (FAS)
- Rayleigh Fading and Spatial Correlation
- Neural Networks and Regression
- Monte Carlo Simulation
- Finite Blocklength Communications
- BLER and Outage Probability Analysis

## Project Goal

The goal is to study how statistical modeling and neural network-based channel prediction can support antenna-port selection in FAS, potentially reducing the need to measure every port while maintaining strong channel performance.

*Note: The scripts implement simulation and evaluation workflows. Performance values should be obtained by running the code in a compatible MATLAB environment; no numerical results are claimed here.*
