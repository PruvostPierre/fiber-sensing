# Fiber Sensing

This project provides a MATLAB-based simulation framework for Distributed Acoustic Sensing (DAS) using optical fibers to model the interaction of light with fiber disturbances via Rayleigh backscattering, enabling the study of dynamic events such as vibrations and other external influences on optical fibers.


---

## Repository Contents

The repository contains the following core files:

- **`DAS_simulation.m`**  
  This is the main script for running the DAS simulation. It integrates all modules to simulate Rayleigh scattering, signal propagation, and the detection of external disturbances.

- **`initialize.m`**  
  Contains the parameter definitions and constants required for the simulation, such as fiber length, wavelength, and sampling rate.

- **`RayleighModel.m`**  
  Defines the mathematical model used to calculate Rayleigh scattering coefficients for the optical fiber.


- **`genGolayCode.m`**  
  Generates Golay complementary codes used for probing sequences in the DAS simulation. This function produces mutually orthogonal pairs of sequences that are essential for effective signal detection.


- **`functions/`**  
  This folder contains additional subfunctions, used to simulate the DAS system's behavior (fiber dynamics, computing polarization effects etc.)


---

## Requirements

- MATLAB R2020b or later.

---

## Usage

### Setting Up

