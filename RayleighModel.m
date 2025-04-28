function [p,r]=RayleighModel(p,r)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates the probing sequence, models the sensed fiber, 
% the transmitter and the receiver frontends, and performs the correlation
% process and the Jones matrices extraction at the RX side.
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Erx; % global variable containing the propagated field

%% Start timing for execution
tStart = tic; % Start timer for the entire process
tCumul = 0;   % Initialize cumulative time counter

%% Probing

% Generate the probing sequence (codes for transmission)
[gCode, p] = genProbingSequence(p);
p.tx.gCode = gCode; % Store the generated probing sequence in the parameter structure

% Estimate initial reception length and set sampling frequency
p.rx.ErxLen = length(p.tx.gCode)* p.tx.nbCodes;
p.rx.fSamp = p.tx.ovsFactor * p.tx.fSymb;

% Handle polarization modes based on user configuration
if strcmpi(p.ProbingMode, 'simo') || strcmpi(p.ProbingMode, 'siso')
    if p.tx.Xpol
        gCode(2, :) = zeros(1, size(gCode, 2)); % Enable X polarization only
    else
        gCode(1, :) = zeros(1, size(gCode, 2)); % Enable Y polarization only
    end
end

% Handle OFDM or single-carrier probing sequence
if p.tx.OFDM
    % Convert codes to serial OFDM symbols for subcarriers
    tx = reshape([gCode(1, :); zeros(p.tx.subcarriers - 1, size(gCode, 2))], 1, []);
    ty = reshape([gCode(2, :); zeros(p.tx.subcarriers - 1, size(gCode, 2))], 1, []);
    gCodeSingle = cat(1, tx, ty);
else
    gCodeSingle = gCode; % Use single carrier probing sequence
end

% Normalize probing signal to 1 mW per polarization
P_gcode = mean(mean(abs(gCodeSingle.^2), 2));
gCodeSingle = gCodeSingle * sqrt(1e-3 / P_gcode);


% Generate the transmitted optical signal by replicating the sequence
Etx = repmat(gCodeSingle, 1, p.tx.nbCodes); %repeat the sequence number of codes times

% Display power levels if enabled
if p.displ.powerIndication
    fprintf('gCode level at generation: %.2f dBm per polarization\n', 10 * log10(P_gcode * 1e3));
end

%% Generate Rayleigh scattering model and fiber properties
RayStart_time = tic; % Start timing Rayleigh model generation

% Generate scatterers in the fiber 
%p = genRayleighScattering(p);
[p, Hi] = genRayleighScattering(p, gCodeSingle);

%Car passage data 
data_car = data_loading('car passage data\car_filt_high_reduced.mat');
data_car = repelem(data_car, 3, 1)/3; %repeat in space to match approx 1m per segment

%add dynamic birefringence event 
p.pola.segIdx = 1500:(1500 + size(data_car, 1) - 1); % list of segments where the birefringence event occurs
duration = p.tx.nbCodes; %time duration of the event in number of codes : if we want to apply 2 events of different duration, do zero-padding to have same size
p.pola.codeLength = size(gCodeSingle, 2);
if p.pola.betaEvent ==1
    %t = linspace(0, (p.tx.nbCodes-1)*p.tx.Tcode, p.tx.nbCodes); % Time vector for sine wave
    %data = sin(2 * pi * 50 * t)*0.2; % Generate smooth sinusoidal signal with 200 Hz frequency
    for i=1:duration
        data_car_i = squeeze(data_car(:, 800-1+i)); % Extract the data for all segments, for the i-th code (from 800 cause highest disturbance), change the index to get different data (in time)
        if i==1
            HiEvent = birefringenceEvent(p, p.pola.segIdx+1, data_car_i); %add a change in beta at segment segIdx %to put sine wave instead of data_car_i, replace data_car_i by data(i)
            p.pola.HiGen = HiEvent;
            HiRep = [HiEvent, zeros(2, 2 * (size(gCodeSingle, 2) - p.fibre.nbSegments))]; %add zeros to match the size of gCodeSingle
        else
            HiEvent_code = birefringenceEvent(p, p.pola.segIdx+1, data_car_i);
            p.pola.HiGen = [p.pola.HiGen, HiEvent_code]; %store for analysis : comparison between estimation and real data
            HiEvent_code = [HiEvent_code, zeros(2, 2 * (size(gCodeSingle, 2) - p.fibre.nbSegments))]; %add zeros to match the size of gCodeSingle
            HiRep = [HiRep, HiEvent_code]; %add a change in beta at segment segIdx
        end
    end
end

r.HiGen = HiRep;

% Display the time taken for Rayleigh model generation
fprintf('\n* Time to generate Rayleigh scattering: %.2f seconds', toc(RayStart_time));
tCumul = tCumul + toc(RayStart_time);
clear RayStart_time;

%% Generation of a laser phase noise based on a Lorentzian model
lasernoise_start = tic;

if p.tx.lasernoise_on
    [r,p,Etx] = lasernoise_model(Etx,p,r);
end

if strcmpi(p.ProbingMode,'simo') || strcmpi(p.ProbingMode,'siso')
    Etx = 2^(1/2)*Etx; % double the power at the input
end

fprintf('\n * time to generate laser noise: %.2f seconds\n', toc(lasernoise_start));tCumul = tCumul + toc(lasernoise_start);
clear lasernoise_start;

%% Amplification before the fiber
amplinoise_start = tic;
Etx = Etx*10^(-p.tx.mod_loss/20 + p.tx.LaserLevel/20)*sqrt(p.tx.couplerFactor); %laser power + modulator losses
if p.displ.powerIndication
    P_gcode = mean(mean(abs((Etx).^2),2));
    Gcode_level = 10*log10(P_gcode*1e3);%dBm
    fprintf('gCode level at ampli input : %.2f dBm per polarization\n', Gcode_level);
end

if p.tx.ampli_on
    N_ase = (10^(p.tx.gaindB_edfa/10)-1)*p.h*p.f_0*p.tx.nsp_edfa; %edfa noise per polarization
    EDFA_noise = sqrt(N_ase/2)*randn(size(Etx)) + 1j*sqrt(N_ase/2)*randn(size(Etx));
    Etx = Etx*10^(p.tx.gaindB_edfa/20)+EDFA_noise;
end

if p.displ.powerIndication
    P_gcode1 = mean(mean(abs((Etx).^2),2));
    Gcode_level1 = 10*log10(P_gcode1*1e3);% *1e3);%dBm
    fprintf('Signal level at ampli output : %.2f dBm\n', Gcode_level1);
end

Etx = Etx*10^(-p.tx.circ_loss/20); %circulator losses

if p.displ.powerIndication
    P_gcode = mean(mean(abs(Etx.^2),2));
    Gcode_level = 10*log10(P_gcode*1e3);%dBm
    fprintf('Signal level at fiber input : %.2f dBm\n', Gcode_level);
end

fprintf(' \n * time to generate amplification noise %.2f seconds\n', toc(amplinoise_start));
clear P_gcode P_gcode1 Gcode_level Gcode_level1 EDFA_noise N_ase amplinoise_start;

%% Propagation: Apply convolution of Golay codes ('Etx') by fibre Impulse Response ('HiRep')
% Input 1 (Emitted signal): the modulated code Ex & Ey per polar
% Input 2 (sensor array IR): one 2x2 complex matrix per sensor, 
% Output: Rx = hxx*Ex + hxy*Ey & Ry = hyx*Ex + hyy*Ey
% If transmission onto one Polar only: Rx = hxx*Ex, Ry = hyx*Ex

conv_start = tic;

Etx = [zeros(2, size(gCodeSingle, 2)), Etx, zeros(2, size(gCodeSingle, 2))]; % Add zeros to the left and right of the signal to do convolution of each code : the code itself and the neighboring codes (one on the left and one on the right) to account for the mixing between the codes
% Convolution of the transmitted signal with the fiber impulse response
for i=1:p.tx.nbCodes
    Hi_code =  [p.pola.HiGen(:, (i-1)*2*p.fibre.nbSegments+1:i*2*p.fibre.nbSegments), zeros(2, 2 * (size(gCodeSingle, 2) - p.fibre.nbSegments))]; % Jones matrix for the current code
    if i==1
        conv_mix = conv_fft(Etx(:, 1:3*length(gCode)), Hi_code(1,1:2:end), Hi_code(1,2:2:end), Hi_code(2,1:2:end), Hi_code(2,2:2:end)); % Convolve with the first first code and second code (and zeros on the left)
        Erx = conv_mix(:, length(gCode)+1:2*length(gCode)); % keep only part concerning the current code
    elseif i==p.tx.nbCodes
        conv_mix = conv_fft(Etx(:, (i-3)*length(gCode)+1:i*length(gCode)), Hi_code(1,1:2:end), Hi_code(1,2:2:end), Hi_code(2,1:2:end), Hi_code(2,2:2:end)); % Convolve with the last code and second to last code (and zeros on the right)
        Erx = [Erx, conv_mix(:, length(gCode)+1:2*length(gCode))]; % keep only part concerning the current code
    else
        % Convolve with the current code and the two neighboring codes
        conv_mix = conv_fft(Etx(:, (i-2)*length(gCode)+1:(i+1)*length(gCode)), Hi_code(1,1:2:end), Hi_code(1,2:2:end), Hi_code(2,1:2:end), Hi_code(2,2:2:end));
        Erx = [Erx, conv_mix(:, length(gCode)+1:2*length(gCode))]; % keep only part concerning the current code
        
    end
end

fprintf(' \n * time to apply convolution - transmission %.2f seconds',toc(conv_start)); tCumul = tCumul + toc(conv_start);

%% Rx model

rx_start = tic;
%Power levels at coherent mixer input (Local oscillator & Rx signal)
LoLevel = single(p.tx.LaserLevel + 10*log10(1-p.tx.couplerFactor));%dBm, Coherent mixer Local oscillator (LO) intensity level. 3dB: 50% coupler at laser source output attenuation
P_erx_moy = mean(mean(abs(Erx.^2), 1));
RxLevel = single(10*log10(P_erx_moy*1e3));%dBm, Signal intensity level at coherent Rx input. (after convolution)
if p.displ.powerIndication
    fprintf(' Signal level at fiber output : %.2f dBm\n', RxLevel);
end

%Photodiode&TIA : Thermal, shot noise, LO-RIN beat noise
B = p.tx.fSymb;%Bandwidth of the receiver assumed to be equal to the baud rate
alphaCmrr = sqrt(10^(p.rx.alphaCmrrdB/10));%power imbalance lin value for CMRR calculation
CMRR = single(((1-alphaCmrr^2)/(1+alphaCmrr^2))^2);%Common Mode Rejection Ratio
RIN = single(10^(p.tx.RinLevel/10));%in W, RIN of LO
P_Lo = single(p.rx.mixLossFactLo * 1e-3 * (10.^(LoLevel/10)));%in W, power of local oscillator (LO) at each of the 8 coherent mixer outputs
P_Rx = single(p.rx.mixLossFactRx * 1e-3 * (10.^(RxLevel/10)));%in W, power of Rx signal at each of the 8 coherent mixer outputs
ThermalNoiseVar = single(p.rx.measThNoise500MHz*(B/500.e6));%in V^2. Instead of 4*K_*(T_/R_ch)*B_), value actually measured without signal with B_=500MHz at TIA output

% In the following noise variance computations, (P_Rx/2) is considered assuming equal split among polarization states
ShotNoiseVar = single(2*(p.rx.TiaGain^2) * (2*p.e*p.rx.Rsensi*(P_Lo+(P_Rx/2)+2*sqrt(P_Lo*(P_Rx/2)))*B));%in V^2, shot noise viewed after balanced PD+TIA: *sqrt(2) in intensity (=>*2 in variance) due to summation of 2 independent shot noises of same statistics
RinCmrrVar = single(4*(p.rx.TiaGain^2) * (CMRR*(p.rx.Rsensi^2)*(P_Lo^2)*(RIN*2*B)));%in V^2, Rin noise viewed after balanced PD+TIA: *2 in intensity (=>*4 in variance) due to summation of 2 correlated Rin noises (assumption here: full correlation)
SNR_elect = single((p.rx.TiaGain^2)*(p.rx.Rsensi^2)*0.5*(P_Rx/2).*P_Lo./(ThermalNoiseVar+ShotNoiseVar+RinCmrrVar));%SNR estimated at PD+TIA output
p.snr_dB = single(10*log10(SNR_elect));%dB expression of SNR per polar as viewed at coherent mixer output

if p.rx.coherentdetection_on
    %Coherent receiver emulation (mixer + balanced photo-detection + TIA)
    % Creating LO - Homodyne detection: no detuning / Parameters: LO power and phase noise
    E_lo = sqrt(P_Lo).*ones(2,p.rx.ErxLen); 
    if p.rx.lasernoise_on
        % CHECK IF LASER PHASE NOISE IS WELL APPLIED AT THE RX SIDE
        laserNoiseMat = r.laserNoiseMat;
        E_lo = E_lo.*laserNoiseMat;%(:,1:hi_size); % Local oscillator at mixer output
    end
    Erx = sqrt(p.rx.mixLossFactRx) * Erx; % Er at mixer output (*2 since input is after 1st PBS (polar sep.))
    
    if p.displ.powerIndication
        P_erx_moy = mean(mean(abs(Erx.^2),2)); %Watt
        fprintf(' Signal level at coherent mixer output : %.2f dBm\n', 10*log10(P_erx_moy*1e3));
    end
    
    %Computing the 8 optical fields at mixer output
    E_cr = [Erx(1,:)*1i + E_lo(1,:)*1i; -Erx(1,:) + E_lo(1,:); Erx(1,:) + E_lo(1,:)*1i; Erx(1,:)*1i + E_lo(1,:);...   %in W. 1/2*Er(k,:): half of the optical field per polar
        Erx(2,:)*1i + E_lo(2,:)*1i; -Erx(2,:) + E_lo(2,:); Erx(2,:) + E_lo(2,:)*1i; Erx(2,:)*1i + E_lo(2,:)];
    
    % Computing the 8 photo-currents
    I_cr =p.rx.Rsensi*real(E_cr.*conj(E_cr)); %in A, square of module of each of the 8 fields multiplied by photodiod sensivity % p.rx.Rsensi*
    
    % Applying balanced photo-detection
    I_cr = [I_cr(1,:)-I_cr(2,:); I_cr(3,:)-I_cr(4,:); I_cr(5,:)-I_cr(6,:); I_cr(7,:)-I_cr(8,:)];%in A
    
    % Computing voltages after TIA
    Erx = p.rx.TiaGain*[I_cr(1,:)+1i*I_cr(2,:); I_cr(3,:)+1i*I_cr(4,:)];%in V 
    
   
    % Blocking DC component (to emulate DC block of each of the 4 balanced PD outputs)
    if p.rx.RxDcBlocker~=0
        Erx(1,:) = Erx(1,:) - mean(Erx(1,:),2);
        Erx(2,:) = Erx(2,:) - mean(Erx(2,:),2);
    end
    
    if p.displ.powerIndication
        P_erx_moy = mean(mean(abs(Erx.^2), 2));
        fprintf(' Signal level at TIA and photodiode output before awgn : %.2f dBV2\n', 10*log10(P_erx_moy));
    end
    
end

% Adding AWGN noise
if p.rx.awgnRX_on
    var_awgn = (ThermalNoiseVar+ShotNoiseVar+RinCmrrVar);
    awgn = sqrt(var_awgn/2)*complex(randn(size(Erx)),randn(size(Erx)));%Complex AWGN generation (one per polar): /2 since per Real/Imag dimension
    Erx = Erx + awgn;
    if p.displ.powerIndication
        P_erx_moy = sum(mean(abs(Erx.^2), 1)/2); %Watt
        P_awgn_moy = sum(mean(abs(awgn.^2), 1)/2); %Watt
        fprintf(' Signal level at TIA and photodiode output : %.2f dBm, Noise level %.2f dBm \n', 10*log10(P_erx_moy), 10*log10(P_awgn_moy));%*1e3));
    end
end

fprintf(' \n * time to go through Rx model %.2f seconds \n',toc(rx_start)); tCumul = tCumul + toc(rx_start);
p.rx.ErxLen = size(Erx,1);

clear RxLevel P_erx_moy awgn var_awgn I_cr LoLevel B_ alphaCmrr P_awgn_moy CMRR RIN_ ...
    P_Lo P_Rx SNR_elect RinCmrrVar ShotNoiseVar ThermalNoiseVar E_lo E_cr rx_start;

%% Multicarrier case %FIXME
corr_start = tic;

%% Correlation process at the reception
p.rx.ErxLen = size(Erx,2);
Erx=Erx.';
p.tx.detectionindex = [0 0]; %For synchronization and detection of fibre start and fibre end


for k = 1:p.tx.subcarriers
    p = rxSensingCorrelation(p); % Extract estimated Jones matrices
    p.HiTabSB(:,:,k) = p.HiTab;%Store Jones matrices for the k-th subband
    p.rx.nbOvsRayleighReflectorsSB(k) = p.rx.nbOvsReflectors;%Store nb of detected segments
    p.rx.nbDetectedCodes = floor(size(p.HiTab, 2)/p.rx.nbOvsReflectors);%Nb detected codes
    p.rx.nbDetectedCodesSB(k) = p.rx.nbDetectedCodes;
    
    fprintf('\n * time spent for parameters extraction (%d subbands case): %.2f seconds',p.tx.subcarriers, toc(corr_start)); tCumul = tCumul + toc(corr_start);
    fprintf('\n*** Rx PROC. OUTPUTS (SUBBAND %d/%d) Nb detected segments:%d/%d  Nb detected codes (%d symbols spacing):%d  Overall elapsed time:%.2f seconds ***\n Cumulated time %.2f seconds\n', ...
        k, p.tx.subcarriers, p.rx.nbOvsReflectors, floor(p.fibre.L/p.fibre.spatialRes), p.tx.Ncode, p.rx.nbDetectedCodes, toc(tStart), tCumul);
end

clear global Erx; %clear global variable here since it is not needed any longer from now


p.HiGen = [Hi(1,1:2:end); Hi(2,1:2:end); Hi(1,2:2:end); Hi(2,2:2:end)]; %Simulated matrices (before TR) %FIXME Should we reduce the number of segments after synchro?
