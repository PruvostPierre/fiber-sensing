function [gSweep,p]  = genSweep(p)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generation of two sweep signals for each of the two polarization channels
%Input parameters:
%p.tx.F0sweep: sweep start frequency 
%p.tx.F1sweep: sweep stop frequency 
%p.tx.MUsweep: Sweep curvature param (0<MU<0.5)
%p.tx.DFsweep: normalized frequency between 2 parallel linear sweeps (0<DF<0.5)
%
%Output: gSweep: the 2 sweep signals for probing the two polar axes
%(sweep length is fixed to match the same length as orth. code)
%
% Authors: 
% Original code by C. Dorize - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




p.tx.Ncode = 8*2^p.tx.seqOrderCst; %sweep length, matches length of the codes
p.tx.Tcode = p.tx.Ncode/p.tx.fSymb;%in s, time sounding duration for 1 code

tVect = (1/p.tx.fSymb)*(0:p.tx.Ncode-1); 

gSweep(1,:) = sin(2*pi * ( (((p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + p.tx.F0sweep.*tVect) );%Linear sweep, freq increases with time
gSweep(2,:) = ([gSweep(1,p.tx.Ncode/2+1:end) gSweep(1,1:p.tx.Ncode/2)]);%Translation (circshift) by half the sweep length

% switch(p.tx.ProbingMethod)
%     case(p.tx.SWEEP_CROSSED_PROBING)  %2 linear sweep signals with opposite slope
%         gSweep(1,:) = sin(2*pi * ( (((p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + p.tx.F0sweep.*tVect) );%Linear sweep, freq increases with time
%         gSweep(2,:) = sin(2*pi * ( ((-(p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + p.tx.F1sweep.*tVect) + p.PhiSweep2);%Linear sweep, freq decreases with time
%     case(p.tx.SWEEP_CURVED_PROBING)   %2 curved (MUsweep param dependent) sweep signals
%         gSweep(1,:) = sin(2*pi * ( (((p.tx.F1sweep-p.tx.F0sweep)/(p.tx.Tcode^(1+p.tx.MUsweep)))/(1+(1+p.tx.MUsweep)))*(tVect.^(1+(1+p.tx.MUsweep))) + p.tx.F0sweep.*tVect) );%Concave sweep, order (1+mu)
%         gSweep(2,:) = sin(2*pi * ( (((p.tx.F1sweep-p.tx.F0sweep)/(p.tx.Tcode^(1/(1+p.tx.MUsweep))))/(1+(1/(1+p.tx.MUsweep))))*(tVect.^(1+(1/(1+p.tx.MUsweep)))) + p.tx.F0sweep.*tVect) );%Convex sweep, order 1/(1+mu)
%     case(p.tx.SWEEP_PARALLEL_PROBING) 
%         df = (p.tx.F1sweep-p.tx.F0sweep)*p.tx.DFsweep;
%         gSweep(1,:) = sin(2*pi * ( (((1-p.tx.DFsweep)*(p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + p.tx.F0sweep.*tVect) );%Linear sweep, freq increases with time
%         gSweep(2,:) = sin(2*pi * ( (((1-p.tx.DFsweep)*(p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + (p.tx.F0sweep+df).*tVect) );%Linear sweep, freq increases with time        
%     case(p.tx.SWEEP_TRANSLATED_PROBING) %2 linear sweep signals with half-duration translation
%         %df = (p.tx.F1sweep-p.tx.F0sweep)*p.tx.DFsweep;
%         gSweep(1,:) = sin(2*pi * ( (((p.tx.F1sweep-p.tx.F0sweep)/p.tx.Tcode)/2)*(tVect.^2) + p.tx.F0sweep.*tVect) );%Linear sweep, freq increases with time
%         gSweep(2,:) = ([gSweep(1,p.tx.Ncode/2+1:end) gSweep(1,1:p.tx.Ncode/2)]);%Translation (circshift) by half the sweep length
% end

end
