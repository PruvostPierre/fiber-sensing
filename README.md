# Fiber Sensing

This project provides a MATLAB-based simulation framework for Distributed Acoustic Sensing (DAS) using optical fibers. The objective is to model the interaction of light with fiber disturbances via Rayleigh backscattering, enabling the study of dynamic events such as vibrations and other external influences on optical fibers.  

---

## Project Overview  

Distributed Acoustic Sensing (DAS) is an advanced optical sensing technique that utilizes Rayleigh backscattering within optical fibers to detect and monitor environmental changes. This simulation captures key phenomena, including signal propagation, Rayleigh scattering, phase changes, and polarization effects, to create a realistic model of DAS behavior.  

---

## Repository Contents  

The repository contains the following core files:  

- **`DAS_simulation.m`**  
  This is the main script for running the DAS simulation. It integrates all modules to simulate Rayleigh scattering, signal propagation, and the detection of external disturbances.  

- **`genRayleighScattering.m`**  
  Generates Rayleigh backscatter coefficients along the fiber, accurately modeling the random scattering inherent to optical fibers.  

- **`initialize.m`**  
  Contains the parameter definitions and constants required for the simulation, such as fiber length, wavelength, and sampling rate.  

- **`RayleighModel.m`**  
  Defines the mathematical model used to calculate Rayleigh scattering coefficients for the optical fiber.  

- **`fcorrAndGetJones.m`**  
  Computes Jones matrices and accounts for polarization and phase changes due to fiber dynamics.  

---

## Features  

- **Dynamic Event Simulation**: Models the effects of vibrations and disturbances on the optical fiber.  
- **Rayleigh Backscatter Modeling**: Simulates realistic scattering effects based on fiber characteristics.  
- **Polarization and Phase Analysis**: Tracks changes in polarization and phase along the fiber length.  
- **Customizable Parameters**: Provides flexibility to adjust fiber properties and simulation parameters.  

---

## Requirements  

- MATLAB R2020b or later.  

---

## Usage  

### Setting Up  
1. Clone the repository:  
   ```bash  
   git clone <repository-url>  
   cd <repository-folder>  
