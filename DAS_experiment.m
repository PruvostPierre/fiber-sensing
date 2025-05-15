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
addpath('E:\03-Code\functions'); %CHANGEME
addpath('E:\03-Code');

% Initialize environment and simulation parameters
p.env = 'exp';              % 'model' for simulation, 'exp' for experimental data
p.tx.fSymb = 100e6;          % Symbol rate [Baud]
p.tx.ovsFactor = 1;
p.fibre.L = 5200; % fiber length in meters [m] 
p = initialize(p);          % Initialize default parameters

%% Override default parameters

% Transmission and reception configuration
p.ProbingMode = 'SIMO';     % 'SISO', 'SIMO', 'MISO', or 'MIMO'
p.tx.Xpol = 1;              % Use cross-polarization at TX
p.rx.Xpol = 1;              % Use cross-polarization at RX

% Transmitter parameters
p.tx.ProbingMethod = 'Golay';      % Probing sequence type
p.tx.seqOrderCst = 12;             % Sequence order for Golay codes

% Stored data files
p.rx.data_directory = 'E:\02-DATA\01-Test\';
p.rx.data_filename = 'Acquire5';

% Display parameters
p.displ.fIdx = 10;                 % Figure index for display
p.displ.powerIndication = 1;       % Display power information (1: enabled)
p.displ.lowResolFactor =  50;

%% Pre-allocate Variables for Analysis
% Initialize variables to store simulation results
mseTab = [];                      % Placeholder for Mean Squared Error results
p.stockStd = [];                  % Placeholder for standard deviation of differential phase
p.stockabsDet = [];               % Placeholder for determinant values of Jones matrices
r.threshold_comb_stock = [];      % Placeholder for detection thresholds
r.nbRemovedLastSegments = 1;      % Segments removed from analysis (default: 1)

%% Golay Code 
    [gCode, p] = genProbingSequence(p);
    %gCode(2, :) = gCode(1, :);
    p.tx.gCode = gCode;

%% Main Processing for Acquired Sequence
    % Correlation and Jones matrix estimation
    [p, r] = Fileread_and_correlation(p, r);

    % Post-process data: differential phase, filtering, SOP calculation
    [p, r] = getPostProcessingfromJones(p, r);

    % Store results for analysis
    p.stockStd = cat(2, p.stockStd, p.stdDiffPhiTabSelect(2:end-1));

    if (p.ProbingMode =='MIMO')
        p.stockabsDet = cat(2, p.stockabsDet, r.selectAbsDetTab_nonorm(2:end-1));
        r.threshold_comb_stock = cat(2, r.threshold_comb_stock, r.threshold_comb);
    end

%% Display Results
% Visualize Rayleigh detection analysis
p = displayRayleighDetection(p, r);
