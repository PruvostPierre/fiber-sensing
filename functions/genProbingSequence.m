function [gCode,p] = genProbingSequence(p)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%'genProbingSequence' generates probing code (coding of two mutually 
% orthogonal pairs of complementary mates) or sweep probing codes/sweep generation 
%
% FIXME For OFDM case, the symbol rate fSymb and the nb of taps Ncode in each
% subcarrier are divided by 2^n with (n>0) in order to come up, 
% after OFDM modulation, with a data flow at the % same rate as the single 
% carrier reference case. Therefore, in the multicarrier case, fSymb 
% stands for the symbol rate of the OFDM flow.
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Select probing method based on user input
switch (lower(p.tx.ProbingMethod))
    
    % Golay Code: Generates two mutually orthogonal pairs of complementary codes
    case 'golay'
        [gCode, p] = genGolayCode(p); % Probing codes generation using Golay code method
        
        maxLengthFactor = (1/4); % Golay codes are more compact, so maximum length factor is 1/4
        
    % Cazac Code: Uses two translated PSK sequences for probing
    case 'cazac'
        % Calculate the probing sequence length and the time duration for one sequence
        p.tx.Ncode = 8 * 2^p.tx.seqOrderCst; % Sequence probing length based on sequence order
        p.tx.Tcode = p.tx.Ncode / p.tx.fSymb; % Time duration for one sequence (in seconds)
        
        % Generate the Cazac code for probing, with complex exponential components
        gCode(1, :) = exp(1i * 2 * pi / sqrt(p.tx.Ncode) .* (mod((1:p.tx.Ncode) - 1, sqrt(p.tx.Ncode)) + 1) .* (floor(((1:p.tx.Ncode) - 1) / sqrt(p.tx.Ncode)) + 1)) * exp(1i * pi / 4);
        gCode(2, :) = circshift(gCode(1, :).', p.tx.Ncode / 2).'; % Shift the second sequence by half the length
        
        maxLengthFactor = (1/2); % Cazac codes are less compact, so maximum length factor is 1/2
       
        
    % Default case: If an unknown probing method is provided, an error is raised
    otherwise
        error('ERROR: Unknown probing method.')
        
end

% Oversample the probing code to match the sampling rate
gCode = ([reshape(repmat(gCode(1, :), p.rx.ovsFactor, 1), 1, []); reshape(repmat(gCode(2, :), p.rx.ovsFactor, 1), 1, [])]); 

% Display probing sequence and fiber details

% fprintf('\n*** Rx PROC. INPUTS  Probing:%s  fSymb:%.0fMHz  Ncode:%d symbols (sb%d)  Tcode:%.2fus  BW:%.2fkHz  MaxProbingDist:%.3fkm  EnteredFiberLength:%.1fkm  Signal duration:%.4fs ***\n', ...
%     p.tx.ProbingMethod, p.tx.fSymb * 1.e-6, p.tx.Ncode, p.tx.seqOrderCst, p.tx.Tcode * 1.e6, 0.5e-3 / p.tx.Tcode, 1.e-3 * maxLengthFactor * (0.5 * p.tx.Tcode * p.fibre.cFiber), 1.e-3 * p.fibre.L, 8 * 2^(p.tx.seqOrderCst) * p.tx.nbCodes / (p.rx.ovsFactor * p.tx.fSymb));
fprintf('\n*** Rx PROC. INPUTS  Probing:%s  fSymb:%.0fMHz  Ncode:%d symbols (sb%d)  Tcode:%.2fus  BW:%.2fkHz  MaxProbingDist:%.3fkm  EnteredFiberLength:%.1fkm ***\n', ...
    p.tx.ProbingMethod, p.tx.fSymb * 1.e-6, p.tx.Ncode, p.tx.seqOrderCst, p.tx.Tcode * 1.e6, 0.5e-3 / p.tx.Tcode, 1.e-3 * maxLengthFactor * (0.5 * p.tx.Tcode * p.fibre.cFiber), 1.e-3 * p.fibre.L);


% Warning: Check if probing sequence length is compatible with fiber length
if p.tx.Tcode < (1 / maxLengthFactor) * 2 * p.fibre.L / p.fibre.cFiber || p.tx.Tcode > 1 / (pi * p.tx.dfLaser)
    fprintf('\nWARNING: Probing sequence (type: %d, order %d) length (%.1fkm) can process fiber length up to %.3fkm (current one is %.1fkm) to avoid spatial aliasing AND shorter than laser coherence length (%.0fkm). Potential inconsistent phase detection!\n\n', ...
        p.tx.ProbingMethod, p.tx.seqOrderCst, 1.e-3 * p.fibre.cFiber * p.tx.Tcode, 1.e-3 * maxLengthFactor * (0.5 * p.tx.Tcode * p.fibre.cFiber), 1.e-3 * p.fibre.L, 1.e-3 * p.fibre.cFiber / (pi * p.tx.dfLaser));
end

