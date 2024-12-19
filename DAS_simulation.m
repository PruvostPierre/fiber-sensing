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

clear global; %clear global variables p and r;
clear all;
close all;

% p: structure containing simulation parmaters and initialized in
% initialize.m
% r: structure containing generated results

% Path addition for subfunctions
addpath('D:\Notebook\DAS\DAS model 230924 - Copy\functions');
addpath('D:\Notebook\DAS\DAS model 230924 - Copy');

p.env = 'exp'; % Environnement: 'model' for simulation or 'exp' for experimentally acquired data
p.fibre.L  =        1000;  % Length of sensed fiber in meters [m]
p.tx.fSymb =        50e6;  % Symbol rate [Baud]
p.tx.ovsFactor =       1; %1 in simulation model
p = initialize(p);

%% Override default parameters

%Single or dual-polarization at TX? RX?
p.ProbingMode =  'MIMO'; % 'SISO', 'SIMO', 'MISO';
p.tx.Xpol =  1; % if SISO or SIMO, is Xpol used at TX?
p.tx.Xpol =  1; % if SISO or MISO, is Xpol used at RX?

p.stat_NB = 1; % Number of generated fibres over which stats are done

%Fibre
p.fibre.ScatDensity =   50; % nb of backscattering points per spatial segment []
p.fibre.polCorrL =     20; % Polarization beat length between 0.05 and 100m. [m]

%Parameters for Dynamic case (external vibration): a single fiber segment is excited if ExcitedSegmentFlag=1
p.fibre.ExcitedSegmentFlag = 0;     %A flag to change the fiber response as fct of time at one segment position
p.fibre.ExcitedSegmentIdx =  200;   % index of segment over which the vibration is applied
p.fibre.ExcitedStrainMax = 300e-9;   % Maximal fiber extension per meter induced by mechanical event [m/m]
p.fibre.ExcitedF_event = 50;        % Frequency of the pure sine wave perturbation [Hz]
p.fibre.ExcitedDynEvolution = 0;    % 0 (default), +1 or -1:   +1(resp.-1) linear amplitude increase over time (resp. decrease)
p.fibre.artificialFading = 1;       %add extra fading at perturbation location (1= no extra fading)

%Transmitter
p.tx.ProbingMethod =  'Golay';
p.tx.seqOrderCst =  11;%Number of basis iterations to get the final sequences
p.tx.nbCodes =   300;% Number of transmitted codewords
p.tx.dfLaser =  0;% Hz, FWHM laser linewidth

%Receiver
p.tx.lasernoise_on = false; %Add laser noise or not?
p.tx.ampli_on = false; % EDFA amplifier on or off?
p.tx.lasernoise_on = false; %Add laser noise only at LO or not?
p.tx.awgnRX_on = false; %apply AGN at RX?
p.tx.offset_ratio = 1; % ratio with respect to one code, for initial offset after correlation, minimum value: 1

% Processing and display parameters
p.displ.fIdx =  10;
p.displ.powerIndication =  1; %(model) prints in and out power for each block in cmd window
%p.displ.lowResolFactor =  10; %1;%Coarse spatial resolution factor used during initial diff phase calculation (100 =>1/100 of the rayleigh backscatters are selected
%p.displ.highestIntensSelectionRatio =     0.1;%0:none 1:all, ratio for selection of highest intensity reflectors


%% Seed generation (Model only)
rng('shuffle');
p.seeds.seed_scatMag = rng('shuffle');
p.seeds.seed_scatDist = rng('shuffle');

p.seeds.seed_theta = rng('shuffle');
p.seeds.seed_rotPol = rng('shuffle');
p.seeds.seed_beta = rng('shuffle');
p.seeds.seed_gamma = rng('shuffle');

%% CHECK these parameters
mseTab = []; %FIXME
%p.dfLTab = [0 0.005 0.1 0.5 0.75 1 5 7 10 15 30 50 75 100 150 200 400];
%p.seqOrderCsttab = [6 8 10 12 14 16];
p.stockStd = []; %FIXME
p.stockabsDet =[]; %FIXME
r.threshold_comb_stock = []; %FIXME
r.nbRemovedLastSegments =                 1;

%% Main loop
for n=1:p.stat_NB
    
    %Change the seed between runs
    %p.seeds.seed_scatMag = rng(6);
    %p.seeds.seed_scatDist = rng(7);
    %p.seeds.seed_theta = rng(8);
    %p.seeds.seed_rotPol = rng(12);%6'shuffle'
    %p.seeds.seed_beta = rng(9);%3'shuffle'
    %p.seeds.seed_gamma = rng(123845);
    
    % TX + fiber + RX + correlation and estimation of Jones matrices
    [p,r]=RayleighModel(p,r);
    
    % Post-processing: differential phase, windowing, filtering, SOP computation
    [p,r] = getPostProcessingfromJones(p,r);
    
    %Reliability trajectories  %FIXME not tested
    %p = displayRelPhi2D(p,r) ;
  
    % Store Postproc parameters
    p.stockStd = cat(2, p.stockStd,p.stdDiffPhiTabSelect(2:end-1));
    p.stockabsDet = cat(2, p.stockabsDet,r.selectAbsDetTab_nonorm(2:end-1));%/r.threshold_comb);
    r.threshold_comb_stock = cat(2,r.threshold_comb_stock, r.threshold_comb);
    
    %     p.displ.selectedIdxTabSB(n,:) =                           p.displ.selectedIdxTab;
end


%% Graphics to display results
p = displayRayleighDetection(p,r);%display function

%mseTab = cat(1, mseTab, [p.mseDetOut p.mseNormOUT p.mseDetOUTbs p.mseNormOUTbs p.tx.dfLaser]);
