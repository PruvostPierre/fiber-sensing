function [gCode,p]  = genGolayCode(p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Complementary pair sequences generation 
% Sequence length, inter-seq length and sampling rate to be jointly chosen according to the length of the impulse response: N.T>4.Tir with N=2.(NG+Ndz)
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the Golay sequence basis is set to 0
if p.tx.seqBasis == 0 
    % Define the complementary pair mates for Golay sequences
    % These are the two pairs of mutually orthogonal codes for the Golay sequence
    ga1 = [1, -1, -1, -1]; 
    gb1 = [-1, 1, -1, -1]; 
    ga2 = [-1, -1, 1, -1]; 
    gb2 = [1, 1, 1, -1]; 
end

% Determine the sequence order based on the modulation type
if strcmpi(p.tx.modulation, 'qpsk')
    % For QPSK modulation, the sequence order is adjusted to ensure same code length as BPSK
    seqOrder = p.tx.seqOrderCst + 1; 
else
    % For other modulation schemes, use the provided sequence order
    seqOrder = p.tx.seqOrderCst;
end

% Generate the two complementary Golay sequences based on the selected sequence basis and order
[gxat, gxbt] = GenGolay(p.tx.seqBasis, seqOrder, ga1, gb1); % Generate the first Golay pair
[gyat, gybt] = GenGolay(p.tx.seqBasis, seqOrder, ga2, gb2); % Generate the second Golay pair

%% Mapping and modulation of the codes onto the 2 polarization axes
if strcmpi(p.tx.modulation, 'bpsk') 
    % BPSK Modulation (Embodiment 1 & 2):
    % One code per polarization, with complementary sequences sent successively
    
    % Check code time length compatibility (BPSK case)
    % If nbSymbFbgSep is not a multiple of 4, Tcode must be longer than 4*Tir
    NG = length(gxat); % Length of the Golay sequence
    Ndz = round(NG * p.tx.dead_zone); % Sequence and inter-sequence lengths
    p.tx.Ncode = 2 * (NG + Ndz); % Total code length for BPSK modulation
    
    % Map complementary Golay sequences onto two polarizations
    gCode = zeros(2, p.tx.Ncode); 
    gCode(:, 1:NG) = [gxat; gyat]; 
    gCode(:, NG + Ndz(1) + 1:2*NG + Ndz(1)) = [gxbt; gybt]; % Successively send the complementary sequences
    
elseif strcmpi(p.tx.modulation, 'qpsk') 
    % QPSK Modulation (Embodiment 3):
    % One code per polarization, with both complementary sequences mapped together in QPSK
    
    % Check symbol rate compatibility (QPSK case)
    NG = length(gxat); % Length of the Golay sequence
    Ndz = round(NG * p.tx.dead_zone); % Sequence and inter-sequence lengths
    p.tx.Ncode = NG + Ndz; % Total code length for QPSK modulation
    
    % Map the Golay pairs to QPSK, combining both sequences for each polarization
    gCode = zeros(2, p.tx.Ncode); 
    gCode(:, 1:NG) = (sqrt(2) / 2) * ([gxat + 1i * gxbt; gyat + 1i * gybt]); % QPSK mapping
    
elseif strcmpi(p.tx.modulation, 'emb4') 
    % EMB4 Modulation (Embodiment 4):
    % QPSK with factor 2 oversampling and sign inversion for odd bits
    
    % Check symbol rate compatibility (EMB4 case)
    NG = 2 * length(gxat); % Length of the Golay sequence with factor 2 oversampling
    Ndz = 2 * round(NG * p.tx.dead_zone); % Sequence and inter-sequence lengths including oversampling
    p.tx.Ncode = NG + Ndz; % Total code length for EMB4, taking oversampling into account
    
    % Perform oversampling and sign inversion for the sequences
    gxat = reshape([gxat; gxat], 1, []); gxbt = reshape([gxbt; -gxbt], 1, []); % Repeat and invert Gb symbols
    gyat = reshape([gyat; gyat], 1, []); gybt = reshape([gybt; -gybt], 1, []); % Repeat and invert Gb symbols
    
    % Map the oversampled sequences to QPSK modulation
    gCode = zeros(2, p.tx.Ncode); 
    gCode(:, 1:NG) = (sqrt(2) / 2) * ([gxat + 1i * gxbt; gyat + 1i * gybt]); % QPSK mapping of oversampled signal
end

% Calculate the time duration for one code, based on code length and symbol rate
p.tx.Tcode = p.tx.Ncode / p.tx.fSymb; % Time sounding duration for one code (in seconds)

