%% DAS Simulation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of a continuous-wave DAS system + Model for Rayleigh backscattering in a fiber
% - Interrogation modes:  Golay, CAZAC or sweeps
% - Detection modes: SISO, SIMO, MISO and MIMO.
% - Frequency diversity possibility through OFDM
% - Tested only for p.rx.oversampling = 1
%
% Authors:
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

% Clear workspace and command window
clear all;
clear global;
close all;
clc;

% Subfunction Path addition
addpath('D:\fiber-sensing\functions'); %CHANGEME
%addpath('D:\Notebook\DAS\DAS model 230924 - Copy');

% Initialize environment and simulation parameters
p.env = 'exp';              % 'model' for simulation, 'exp' for experimental data
p.fibre.L = 1000;           % Fiber length in meters [m]
p.tx.fSymb =50e6;          % Symbol rate [Baud]
p.tx.ovsFactor = 2;         % Oversampling factor
p = initialize(p);          % Initialize default parameters

%birefringence event
p.pola.betaEvent = 1;  % Birefringence event (0: disabled, 1: enabled)
%% Override default parameters

p.fibre.polar=1;         % Polarization effects (0: disabled, 1: enabled)
% Transmission and reception configuration
p.ProbingMode = 'MIMO';     % 'SISO', 'SIMO', 'MISO', or 'MIMO'
p.tx.Xpol = 1;              % Use cross-polarization at TX
p.rx.Xpol = 1;              % Use cross-polarization at RX

% Statistical simulations
p.stat_NB = 1;              % Number of generated fibers for statistical analysis

% Fiber parameters
p.fibre.ScatDensity = 50;   % Backscattering points per spatial segment []
p.fibre.polCorrL = 20;      % Polarization beat length [m]

% Dynamic fiber excitation (external vibration)
p.fibre.ExcitedSegmentFlag = 0;   % Enable dynamic excitation
p.fibre.ExcitedSegmentIdx = 50;  % Segment index for vibration
p.fibre.ExcitedStrainMax = 300e-9;  % Maximum strain induced by vibration [m/m]
p.fibre.ExcitedF_event = 50;      % Vibration frequency [Hz]
p.fibre.ExcitedDynEvolution = 0;  % Linear amplitude evolution (0: none, +1/-1: increase/decrease)
p.fibre.artificialFading = 1;     % Add artificial fading (1: enabled)

% Transmitter parameters
p.tx.ProbingMethod = 'Golay';      % Probing sequence type
p.tx.seqOrderCst = 11;             % Sequence order for Golay codes
p.tx.nbCodes = 300;                % Number of transmitted codewords
p.tx.dead_zone = 0;
p.tx.dfLaser = 0;                  % Laser linewidth [Hz]
p.tx.lasernoise_on = false;        % Add laser noise at the LO
p.tx.ampli_on = false;             % Enable EDFA amplifier

% Receiver parameters
p.rx.awgnRX_on = false;            % Add AWGN at the receiver
p.rx.offset_ratio = 1;             % Offset ratio after correlation (minimum: 1)
% Display parameters
p.displ.fIdx = 10;                 % Figure index for display
p.displ.powerIndication = 1;       % Display power information (1: enabled)


%% Seed Generation (Model Only)
% Generate random seeds for reproducibility
rng('shuffle');
p.seeds.seed_scatMag = rng('shuffle');
p.seeds.seed_scatDist = rng('shuffle');
p.seeds.seed_theta = rng('shuffle');
p.seeds.seed_rotPol = rng('shuffle');
p.seeds.seed_beta = rng('shuffle');
p.seeds.seed_gamma = rng('shuffle');

%% Pre-allocate Variables for Analysis
% Initialize variables to store simulation results
mseTab = [];                      % Placeholder for Mean Squared Error results
p.stockStd = [];                  % Placeholder for standard deviation of differential phase
p.stockabsDet = [];               % Placeholder for determinant values of Jones matrices
r.threshold_comb_stock = [];      % Placeholder for detection thresholds
r.nbRemovedLastSegments = 1;      % Segments removed from analysis (default: 1)

%% Main Simulation Loop
for n = 1:p.stat_NB
    % Simulate transmission, fiber propagation, reception, and Jones matrix estimation
    [p, r] = RayleighModel(p, r);

    % Post-process data: differential phase, filtering, SOP calculation
    [p, r] = getPostProcessingfromJones(p, r);

    % Store results for analysis
    p.stockStd = cat(2, p.stockStd, p.stdDiffPhiTabSelect(2:end-1));
    p.stockabsDet = cat(2, p.stockabsDet, r.selectAbsDetTab_nonorm(2:end-1));
    r.threshold_comb_stock = cat(2, r.threshold_comb_stock, r.threshold_comb);
end

%% Display Results
% Visualize Rayleigh detection analysis
p = displayRayleighDetection(p, r);
