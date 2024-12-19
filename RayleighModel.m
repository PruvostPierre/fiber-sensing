function [p,r]=RayleighModel(p,r)
%CHANGED
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

%%
tStart = tic;
tCumul = 0;

%% Probing

% Generate probing sequence
[gCode,p] = genProbingSequence(p); %Probing codes/sweep generation
p.tx.gCode = gCode; % store coded sequence in 'p'

p.rx.ErxLen = p.tx.ovsFactor*8*2^(p.tx.seqOrderCst)*p.tx.nbCodes; %First estimate of ErxLen (recomputed later after reception), used for display purposes
p.rx.fSamp = p.tx.ovsFactor*p.tx.fSymb;

% if SIMO or SISO, choose X or Y code for probing
if strcmpi(p.ProbingMode,'simo') || strcmpi(p.ProbingMode,'siso')
    if p.tx.Xpol % 1
        gCode(2,:) = zeros(1,size(gCode,2));
    else % 0
        gCode(1,:) = zeros(1,size(gCode,2));
    end
end

if p.tx.OFDM 
    tx = reshape(([gCode(1,:); zeros(p.tx.subcarriers-1,size(gCode,2))]),1,[]); % parallel to serial OFDM symbols for both subcarriers across polar X
    ty = reshape(([gCode(2,:); zeros(p.tx.subcarriers-1,size(gCode,2))]),1,[]);
    %gCodeRep = repmat(cat(1,tx,ty),1,p.tx.nbCodes); % repeated sequence for subband1 according to number of codes
    gCodeSingle = cat(1,tx,ty);
    
    
    
else % single carrier
    gCodeSingle = gCode;
end

%normalization to 1mW per polarization
P_gcode = mean(mean(abs(gCodeSingle.^2),2));
gCodeSingle = gCodeSingle*sqrt(1e-3/P_gcode);

p.rx.ErxLen = size(gCodeSingle,2)*p.tx.nbCodes*p.rx.ovsFactor;
vec_size = p.rx.ErxLen; % number of probing symbols
Etx = repmat(gCodeSingle,1,p.tx.nbCodes);

% Oversample Etx and spply rectangular pulse shape
Etx = rectpulse(Etx.',p.tx.ovsFactor).';%CHANGED
%Etx = repelem(Etx, 1, p.tx.ovsFactor); % Oversample each symbol by replicating

if p.displ.powerIndication
    P_gcode = mean(mean(abs(gCodeSingle.^2),2));
    Gcode_level = 10*log10(P_gcode*1e3);%*1e3);%dBm
    fprintf(' gCode level at generation : %.2f dBm per polarization\n', Gcode_level);
end

%% Rayleigh model generation for the fibre of length p.fibre.L as a series of segments whose size is determined by the symbol rate p.tx.fSymb
RayStart_time = tic;

p = genRayleighScattering(p); % Generate fibre scatterers
[Hi,p] = genJonesMatrices(p); % Define random rotations and generate Jones matrices
r.HiGen = Hi;


% Insert Null Jones matrices after Hi to reach a size equal to previous gCodeRep. Necessary to get same convolution results as before
% Repeat the Hi matrices instead of repeating the transmitted codes (equivalence regarding convolution)
% Required to emulate dynamic cases with applied vibrations
%% CHANGED
HiRep = repmat( ([Hi,zeros(2,2*(size(gCodeSingle,2)-p.fibre.nbSegments))]) ,1,p.tx.nbCodes);
%HiRep = [Hi,zeros(2,2*(size(gCodeSingle,2)-p.fibre.nbSegments))];
%HiRep = [HiRep, zeros(2,length(HiRep)*(p.tx.nbCodes-1))]; %WARNING: we should not repeat both gCode and Hi before convolution
% Ensure gCodeSingle has a compatible size relative to Hi
%repFactor = floor(size(gCodeSingle, 2)*0.001 / p.fibre.nbSegments);

%HiRep = repmat(Hi, 1, repFactor);



%Process dynamic model case: recalculates the Jones matrix of the dynamic segment for each transmitted code and stores it at the right place in HiRep
if p.fibre.ExcitedSegmentFlag == 1 %FIXME needs to be done per Tcode
     p.tx.gCodesize = size(gCodeSingle,2);
    [p,HiRep] = genDynamicSegt(p,HiRep,p.fibre.ExcitedSegmentIdx(1),p.fibre.ExcitedStrainMax(1),p.fibre.ExcitedF_event(1),p.fibre.ExcitedDynEvolution(1));
end %excitedSegmentFlag

fprintf('\n* time to generate Rayleigh scatteringv2: %.2f seconds', toc(RayStart_time)); tCumul = tCumul + toc(RayStart_time);
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

Erx = Etx;

conv_fft(HiRep(1,1:2:end),HiRep(1,2:2:end),HiRep(2,1:2:end),HiRep(2,2:2:end));

Erx = (Erx(:,1:size(gCodeSingle,2)*p.tx.nbCodes*p.rx.ovsFactor));%CHANGED

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
    E_lo = sqrt(P_Lo).*ones(2,vec_size); %SG01/21 ovsFactor to use here too
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
    E_cr = [Erx(1,:)*1i + E_lo(1,:)*1i; -Erx(1,:) + E_lo(1,:); Erx(1,:) + E_lo(1,:)*1i; Erx(1,:)*1i + E_lo(1,:);... %in W. 1/2*Er(k,:): half of the optical field per polar
        Erx(2,:)*1i + E_lo(2,:)*1i; -Erx(2,:) + E_lo(2,:); Erx(2,:) + E_lo(2,:)*1i; Erx(2,:)*1i + E_lo(2,:)];
    
    % Computing the 8 photo-currents
    I_cr = p.rx.Rsensi*real(E_cr.*conj(E_cr)); %in A, square of module of each of the 8 fields multiplied by photodiod sensivity
    
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

% if p.tx.subcarriers > 1% Additional actions for multiple subbands or OFDM demodulation: NOT TESTED
%     
%     %     if p.tx.FiltBank %Each subband is oversampled at fSamp
%     %         sbSumRxSP = fft(Erx);
%     %         if p.tx.RrcShapingChoice == 1 %RRC pulse shaping at Tx side
%     %             MaskSb1RxSP = getSpectralRaisedCos(p.tx.ovsFreqDiv*p.rx.ovsFactor, p.tx.rollOffFactor, p.rx.ErxLen).^0.5;%To get a spectral root-raised cosine
%     %         else
%     %             MaskSb1RxSP = ([ones(round(p.rx.ErxLen/(2*p.tx.subcarriers)/p.rx.ovsFactor),1);zeros(p.rx.ErxLen-2*round(p.rx.ErxLen/(2*p.tx.subcarriers)/p.rx.ovsFactor),1);ones(round(p.rx.ErxLen/(2*p.tx.subcarriers)/p.rx.ovsFactor),1)]);%To get a rectangular spectral mask that selects the 1st subband only
%     %         end
%     %
%     %         if p.tx.FiltBank, rxDwnSamplFact = p.tx.ovsFreqDiv; elseif p.rx.FiltBankOvs, rxDwnSamplFact = 1; end %Define the downsampling factor
%     %
%     %         %rec_symb_SB = ones(p.tx.subcarriers*2, floor(p.rx.ErxLen/rxDwnSamplFact)); %(2=polars) will be used to store received symbols per subband (after FFT)
%     %         rec_symb_SBX = ones(floor(p.rx.ErxLen/rxDwnSamplFact),p.tx.subcarriers); %Store received symbols per subband (after FFT) for polar X
%     %         rec_symb_SBY = ones(floor(p.rx.ErxLen/rxDwnSamplFact),p.tx.subcarriers); %Store received symbols per subband (after FFT) for polar Y
%     %         for n=1:p.tx.subcarriers
%     %             sbnRx = ifft(circshift(sbSumRxSP,-round(p.rx.ErxLen*(n-1)/p.tx.subcarriers/p.rx.ovsFactor),1).*repmat(MaskSb1RxSP,1,2));%received subband n signal including base band transposition, sampled at the same rate as the Tx side
%     %             %rec_symb_SB(2*n-1:2*n,:) = [sbnRx(1:rxDwnSamplFact:end,1).'; sbnRx(1:rxDwnSamplFact:end,2).'];%CD 09/2020: WARNING See potential synchro issue here if rxDwnSamplFact>1. '1:p.tx.ovsFreqDiv:end' could be replaced by 'n:p.tx.ovsFreqDiv:end', with 1<=n<=p.tx.ovsFreqDiv
%     %             rec_symb_SBX(:,n) = sbnRx(1:rxDwnSamplFact:end,1);%CD 10/20 code à checker en mode filterbank %CD 9/20: WARNING See potential synchro issue here if rxDwnSamplFact>1. '1:p.tx.ovsFreqDiv:end' could be replaced by 'n:p.tx.ovsFreqDiv:end', with 1<=n<=p.tx.ovsFreqDiv
%     %             rec_symb_SBY(:,n) = sbnRx(1:rxDwnSamplFact:end,2);%CD 10/20 code à checker en mode filterbank
%     %         end
%     %         Erx = ([rec_symb_SBX(:,1), rec_symb_SBY(:,1)]);%Directly load Erx with 1st subband
%     %         p.tx.fSymb = p.tx.fSymb/rxDwnSamplFact; %For displays afterwards
%     %     elseif p.tx.OFDM
%     if p.tx.syncOffset ~= 0
%         Erx = circshift(Erx,p.tx.syncOffset,1);
%     end
%     
%     tic;%reshape according to nb subcarriers and to ovsFactor - WARNING (CD 10/20: valid only with p.rx.ovsFactor=2!!!)
%     size_lim = p.rx.ErxLen - mod(p.rx.ErxLen,2*p.tx.subcarriers);
%     %%% Reorganize Erx to separate polars and recombine subbands while keeping the oversampling factor (here ovsfactor assumed = 2):
%     if p.rx.ovsFactor > 1 %assumed == 2
%         rec_symb_SBX = fft(reshape(([reshape(Erx(1:p.rx.ovsFactor:size_lim,1), p.tx.subcarriers,[]);reshape(Erx(2:p.rx.ovsFactor:size_lim,1), p.tx.subcarriers,[])]),p.tx.subcarriers,[])).';%Size: nbSamplesPerCarrier, (2xnbCarriers)
%         rec_symb_SBY = fft(reshape(([reshape(Erx(1:p.rx.ovsFactor:size_lim,2), p.tx.subcarriers,[]);reshape(Erx(2:p.rx.ovsFactor:size_lim,2), p.tx.subcarriers,[])]),p.tx.subcarriers,[])).';%Size: nbSamplesPerCarrier, (2xnbCarriers)
%     else %model : ovsFact = 1
%         rec_symb_SBX = fft(reshape(Erx(1:size_lim,1),p.tx.subcarriers,[])).';
%         rec_symb_SBY = fft(reshape(Erx(1:size_lim,2),p.tx.subcarriers,[])).';
%     end
%     %%% Reorganize fft signal to form the subbands:
%     Erx = ([rec_symb_SBX(:,1), rec_symb_SBY(:,1)]);%Directly load Erx with 1st carrier
%     toc;
%     p.tx.fSymb = p.tx.fSymb/p.tx.subcarriers;
%     %    end%OFDM
%     
%     p.rx.fSamp = p.rx.ovsFactor*p.tx.fSymb;
%     p.fibre.spatialRes = p.fibre.cFiber/(2*p.tx.fSymb);
%     p.fibre.nbSegments = floor(p.fibre.L_/p.fibre.spatialRes);
% end %subcarriers > 1

%% Correlation process at the reception
p.rx.ErxLen = size(Erx,2);
Erx=Erx.';
p.tx.detectionindex = [0 0]; %For synchronization and detection of fibre start and fibre end

%In case p.rx.FiltBankOvs, correlation applies with the oversampled and pulse shaped version of gCode, and it yields a result sampled at fSamp
% if p.tx.FiltBank % the code to be used for correlation is an oversampled (if p.rx.FiltBankOvs) & pulse shaped version of initial gCode (already oversampled by p.rx.ovsFactor)
%     if p.tx.FiltBank, rxUpSamplFact = 1; elseif p.rx.FiltBankOvs, rxUpSamplFact = p.tx.ovsFreqDiv; end %Define the upsampling factor for the correlation code
%     gCodeRx = ([reshape(repmat(gCode(1,:),rxUpSamplFact,1),1,rxUpSamplFact*size(gCode,2)) ; reshape(repmat(gCode(2,:),rxUpSamplFact,1),1,rxUpSamplFact*size(gCode,2))]).';%Duplicate the codes: codes are twice longer, the now spread over half the bandwidth,in base band: this feeds the first subband
%     if p.tx.RrcShapingChoice == 1
%         rcSP = getSpectralRaisedCos(p.rx.ovsFactor*rxUpSamplFact, p.tx.rollOffFactor, size(gCodeRx,1));%To get a spectral raised cosine
%     else
%         rcSP = ([ones(round(size(gCodeRx,1)/(2*p.tx.subcarriers)/p.rx.ovsFactor),1);zeros(size(gCodeRx,1)-2*round(size(gCodeRx,1)/(2*p.tx.subcarriers)/p.rx.ovsFactor),1);ones(round(size(gCodeRx,1)/(2*p.tx.subcarriers)/p.rx.ovsFactor),1)]);%CD WARNING: bug on upper line corrected
%     end
%
%     gCodeRxRc = ifft(fft(gCodeRx) .* repmat(rcSP,1,2)).';
%     gCode = gCodeRxRc; %The new code to be used for correlation is now an oversampled & pulse shaped version of initial gCode
% end %FiltBank(ovs)

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