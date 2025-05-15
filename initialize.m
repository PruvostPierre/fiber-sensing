function p = initialize(p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization function for DAS_simulation.m
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definition of constants
p.C =      299792458;       % celerity of light in vacuum [m/s]
p.Kg =     1.0129;          % Constant connecting the refractive index and the group refractive index of an SSMF for a wavelength around 1550nm []
p.N =      1.444;           % refractive index in the fibre core of an SSMF []
p.K=       1.3806e-23;      % Boltzmann constant [J/K]
p.T=       273.15+20;       % Temperature [K]
p.e=       1.602176565e-19; % elementary charge of an electron [C]
p.h =      6.626e-34;       % Planck's constant [J.s]

%% General parameters
p.tx.Lambda = 1536.6e-9;    % optical wavelength of the laser source [m]
p.stat_NB = 1;              % Number of generated fibres over which stats are done

p.ProbingMode =  'MIMO';    % 'SISO', 'SIMO', 'MISO';
p.tx.Xpol =  1;             % if SISO or SIMO, is Xpol used at TX?
p.rx.Xpol =  1;             % if SISO or MISO, is Xpol used at RX?

%% Fiber characteristics
p.ng =     p.Kg*p.N;        % Group velocity refractive index of the fibre
p.f_0 =    p.C/p.tx.Lambda; % central frequency of laser source [Hz]
p.fibre.LossdB =   -0.2;    % fibre loss coefficient [dB/km]
p.fibre.polCorrL =    1;   % Polarization correlation length between 0.05 and 100m. Common value for SSMF is 20cm. [m]
p.fibre.cFiber =  p.C/p.ng; % light velocity in the fibre [m/s]
p.fibre.spatialRes =     p.fibre.cFiber/(2*p.tx.fSymb*p.tx.ovsFactor); % spatial resolution induced by symbol rate [m]
                         %(factor 2: backscatter roundtrip in each fiber segment)
p.fibre.nbSegments =     floor(p.fibre.L/p.fibre.spatialRes); % number of fiber segments of length p.fibre.spatialRes 
p.fibre.TxSpatialRes =   p.fibre.spatialRes; % FIX ME

% Rayleigh scatterers
p.fibre.ScatDensity =   100;% nb of backscattering points per spatial segment []
p.fibre.MuScatMag =     0.25*sqrt(1.e-7);% magnitude portion of backscattered light per meter of fibre [1/m] % FIXME Justify values, check Sterenn's thesis
p.fibre.StdScatMag =    sqrt(1.e-9);    % standard deviation of the magnitude backscattered per distance unit for the distribution of scatterers [1/m] % FIXME Justify values, check Sterenn's thesis
p.fibre.fixed =  0;         % positions and amplitude of scaterrers can be changed or freezed

% Polarization effects
p.fibre.polar =  1;         % Should be kept at 1. For specific experiments, this removes polarization rotations % FIXME
p.fibre.AlphaPolRay =   0.00;%part of the polarized ligth that emerges in an orthogonal polarization state at the reflector []

%Seeds (Model only)
rng('shuffle');

%random scattering levels and positions
p.seeds.seed_scatMag = rng('shuffle');
p.seeds.seed_scatDist = rng('shuffle');

% random polarization rotation in first segment
p.seeds.seed_theta = rng('shuffle');
p.seeds.seed_rotPol = rng('shuffle');
p.seeds.seed_beta = rng('shuffle');
p.seeds.seed_gamma = rng('shuffle');

%% Parameters to emulate the addition of a vibration event: simple tone 
p.fibre.ExcitedSegmentFlag = 0;     %0 : no vibration is added, 1: a vibration is added
p.fibre.ExcitedSegmentIdx =  200;   % index of segment over which the vibration is applied
p.fibre.ExcitedStrainMax = 30e-9;   % Maximal fiber extension per meter induced by mechanical event [m/m]
p.fibre.ExcitedF_event = 50;        % Frequency of the pure sine wave perturbation [Hz]
p.fibre.ExcitedDynEvolution = 0;    % 0 (default), +1 or -1:   +1(resp.-1) linear amplitude increase over time (resp. decrease)
p.fibre.artificialFading = 1;       %add extra fading at perturbation location (1= no extra fading)

%% TX and RX parameters
% Transmitter 

% Probing sequence
p.tx.ProbingMethod = 'Golay'; % 'cazac', 'sweep'
p.tx.nbCodes = 50;      % Nb of the transmitted codes. For good standard deviation (std) estimates, use at least 15 codes
p.tx.CAZAC=5;           % FIX ME

p.tx.modulation =  'BPSK'; % for Golay codes or CAZAC probing
p.tx.seqBasis =     0;  % Golay seq basis: 1,10,20,26 (or 0: basis to be predefined)
p.tx.seqOrderCst =  10; % Number of iterations to get the final complementary Golay sequence
p.tx.dead_zone =    0;  % factor that defines the separation zone length when transmitting 2 consecutive complementary sequences

% for OFDM multi-carrier interrogation
p.tx.OFDM = 0; % OFDM on: 1, OFDM off: 0
p.tx.subcarriers =  1; %Number of OFDM subcarriers for probing. By default 1.
p.tx.code_comb =                    0; %0,1,..,6. For same code choose 0, for othogonal code choose any other number up to 6
% 0(codes1&1) 1(codes1&2) 2(codes1&3) 3(codes1&4)
% 4(codes2&3) 5(codes2&4) 6(codes3&4)
p.tx.syncOffset =                   0; %OFDM synchronization
p.rx.OFDMreconstruction =                 1;  %if OFDM trace, choose combination of OFDM subcarriers or full OFDM reconstruction 


% TX impairments
p.tx.lasernoise_on = true; % include or not laser phase noise
p.tx.dfLaser =  100;      % laser linewidth (full width at half maximum) [Hz]
p.tx.LaserLevel =  11;  % intensity level of the laser source [dBm]
p.tx.RinLevel =   -140; % Level of laser Rin (for Rx SNR calc.) [dBm]
p.tx.couplerFactor = 0.7;%coupling ratio between Signal and local oscillator 
p.tx.mod_loss =  8;     % losses of MZM modulator [dB]
p.tx.CirculPdldB = 0;   %PDL of the circulator. PDL emulator bypassed if CirculPdldB=0. [dB]
p.tx.circ_loss =  6;    % circulator insertion loss [dB]
p.tx.ampli_on =  false; % EDFA on or off
p.tx.nsp_edfa =  1.58;  % Spontaneous emission factor of the EDFA 
p.tx.gaindB_edfa =   9; % EDFA gain [dB]

p.tx.detectionindex = [0 0]; %  to be integrated to code % FIXME

% RX: coherent mixer, photodiodes, LO, rx noise
p.rx.PDs_Order = [1 2 3 4];
p.rx.PD_Polarity = sign( p.rx.PDs_Order ); % not used for now
p.rx.Real_Photodiodes = [1 3];%Odd positions are real
p.rx.Imag_Photodiodes = [2 4];%Even positions are imaginary

p.rx.Rsensi =  0.93;    % Rx photodiod sensitivity, ThorLabs PDB480C-AC spec [A/W]
p.rx.alphaCmrrdB = -0.1;% CMRR power imbalance level at coherent balanced mixer [dB}
p.rx.measThNoise500MHz =  1.1818e-05; % Thermal noise variance measured at Kylia mixer + Thorlabs PDB480C-AC output [V^2]
p.rx.mixLossFactLo = 10^(-10.3/10); % LO coherent mixer loss for each of the 8 outputs: -10.3dB loss measured with Kylia versus -9dB theory
p.rx.mixLossFactRx = 10^(-7/10); %Overall Rx coherent mixer loss per polar for each of the 8 outputs: -7dB according to Kylia spec, -6dB theory
p.rx.TiaGain = 13e3;% RF transimpedence gain (50ohm load) [V/A]
p.rx.RxDcBlocker = 0; %Blocking DC component (to emulate high-pass of each of the 4 ThorLabs PD outputs).
p.rx.R_ch = 50;%in Ohm, resistance of charge at Rx coherent mixer: 50 for photodiod, 500 with TIA

p.rx.coherentdetection_on = true; %emulate coherent detection or directly convol/correl.
p.rx.awgnRX_on = true; %apply awgn at the reception
p.rx.lasernoise_on = false; % apply phase noise at RX side
p.rx.dfLaser = p.tx.dfLaser; %laser linewidth (full width at half maximum), self-homodyne, same laser at TX and RX [Hz] 

% DSP parameters
p.rx.apply_total_normalization =  true;
p.rx.fSamp = 200e6; % [Hz]
p.rx.ovsFactor = p.rx.fSamp./p.tx.fSymb; %No oversampling in model; 2; %usual value at RX;
p.rx.resamplingFilter = true;

%Correlation parameters
p.rx.offset_ratio = 1.9; % ratio with respect to one code for initial offset after correlation
p.rx.corrPerBlockTwin = ceil(2*(2^20)/p.tx.fSymb);%50.e-3;% block time window to be defined if correlation is processed per block (s)
p.rx.corrPerBlockGpuFlag =  0; %if 1, the correlation process per block is computed through GPU

%% Post-processing parameters

% interpolation % FIXME NOT TESTED
p.rx.interpolation =  0; %Use a filter before downsampling at rx
p.rx.interpolationRolloff =               1;%0.8; %Spatial filtering, between 0.5 and 1 in practice here
p.displ.spatialfilter =                   0; %in case of oversampling

% decimation % FIXME NOT TESTED
p.rx.decimation = 0; %0 no decimation, 1 for uniform, 2 for lowres mode

% Segment selection criteria
p.rx.crit_detJones =  1; %if 0: criteria is intensity, 1: criteria is determinant of Jones matrix
p.rx.detTheshold =  0.01; %if |det| is comprised between 0 and 1, threshold below which |det| is considered unreliable

% Spatial low resolution parameters
p.displ.lowResolFactor =  10;%Coarse spatial resolution factor used during initial differential phase calculation 
                             %   1: all segments are selected,
                             % X>1: 1 out of X of the Rayleigh backscatters
                             % is selected, the one with the highest
                             % backscarrering
p.rx.averagingLowres = 0; %0 is low res by decimation, 1 is lowres by averaging

p.displ.getSVD =                          0; % Singular value decomposition of estimated matrices

p.displ.edgeRatio =                       0.1; % for differential phase filtering with raised cosine window
p.rx.f_cutoff =                           0;%100; %frequency below which differential phase traces are filtered out to remove DC and low freq components [Hz]
p.rx.f_cutoff_end=                        0;%1/(2*(2^p.tx.seqOrderCst/p.tx.fSymb)); % frequency above which differential phase traces are filtered out to remove high frequency components [Hz]

p.displ.polar =                     0; % compute polarization parameters

p.displ.freqvsTimeDist =                  0; %contour plot with dominant frequencies
p.displ.location =                        0;%990; % beginning of interesting events on the trace [m]
p.displ.maxloc =                          p.fibre.L; %(0.5*p.fibre.cFiber/p.rx.fSamp)*p.rx.nbOvsRayleighReflectors*0.7 [m]

p.displ.averagingvsTime =                 1; %Display averaged intensity & 2D differential phase vs time of selected backscatters (coarse resolution)
p.displ.scaled_map =                      0; % adapt time/distance phase display to enhance visibilty

%% Display options 
p.displ.fIdx =   10;% Random figure index
p.displ.powerIndication =  0; % prints input and output power for each block in command window
p.displ.detection =  1; % in getRayleighSpread, display detection of fiber start and fiber end

p.displ.sphere =                    0; %display poincare sphere
p.displ.polarPSD =                  0; %with perturbation only
p.displ.polarparam =                0;

%FIXME CHECK the parameters below

p.displ.highestIntensSelectionRatio =     0.1; %0:none 1:all, ratio for selection of highest intensity reflectors
%p.displ.cutdist = 1; %CD take value <1 if you do not wish to display the overall fiber length
% if p.tx.f_cutoff = -1 No filtering, = 0 DC filter only, >0 cut_off frequency value
p.displ.minFreq =                         p.rx.f_cutoff;
p.displ.maxFreq =                         1e3; %[Hz]

p.displ.stdPhiAlarmThres =                1;%in rad, thres of phase standard deviation beyond which an alarm is set

p.displ.StDv =                      0; % Display optical phase standard deviation
p.displ.StDvtoStrain =              1;
p.displ.IntensityPerReflector =     1; %Display backscattered intensity per reflector
p.displ.PhiIntensity =              0; %Display Phi as fct of intensity
p.displ.intensityDistance =         1;  p.displ.intensityDistdB = 0; %Display average RBS intensity as fct of fiber distance

p.displ.dispLasernoiseMat =         0; %display psd of laser noise matrix

p.displ.audio =                     0; %detect audio perturbation
p.displ.tabAudio =                  [30];%round(2000/p.displ.lowResolFactor/(p.rx.ovsFactor*logical(p.rx.decimation)))-1];%[375 480];%[13 18];%[10 500];%[202 203];%[9 500];%[86 4986]; %Indices of the segments where a perturbation is introduced
%p.displ.tabAudio =                  ginput(2); p.displ.tabAudio = round(p.displ.tabAudio(:,2).'/(p.displ.lowResolFactor*single(2*logical(p.rx.decimation))*0.5*p.fibre.cFiber/p.rx.fSamp));
p.displ.psd_audio =                 0; %display psd at points where perturbation is applied

p.displ.reliability =           0;  %Display reliability indicators of phase traces
p.displ.reliability_dist =      0;  %Display reliability indicators of phase traces
p.displ.softvalues =            0;  %Display soft bit values for selected segments

p.displ.relDistr   =            0; %distribution (proba and occurences)
p.displ.SNRdistrib =            0; %distribution of SNRdet and SNR phase. Must set a breakpoint.

p.displ.sound =                     0; %listen to sound
p.displ.averaging =                 0; %if statistical approach with N fiber models
p.file.audio1=                      'audioB1.wav';
p.file.audio2=                      'audioB2.wav';

p.displ.pp_mean=                    0;
p.displ.alarms =                    0; %alarms on traces

p.displ.averaging =                 0; %if statistical approach with N fiber models

p.displ.errCalc =                   0;
p.displ.errCalcDispl =              0;
p.displ.errDet =                    0; %|det|/AiPi

p.displ.ploterrors =                0; %plot absolute and relative errors

