%p = genRayleighScattering(p);
function [p, Hi] = genRayleighScattering(p, gCodeSingle)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates the Rayleigh scatterers in the fiber segments 
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1)Draw scatterer positions and attenuations along the fibre
avgScatStep = p.fibre.TxSpatialRes/p.fibre.ScatDensity; %/!\ TxSpatialRes -- in m, average distance between 2 adjacent scatterers
nbScats = floor(p.fibre.L/avgScatStep);%nb backscatters over the fibre

rng(p.seeds.seed_scatDist);
scatDistTab = single(cumsum(rand(nbScats,1)));%
p.fibre.scatDistTab = p.fibre.L*(scatDistTab./max(scatDistTab));%in m, table of scatterer positions ('zm') over the whole fibre length (uniform random distribution)

rng(p.seeds.seed_scatMag);
p.fibre.scatMagTab = single((p.fibre.MuScatMag+p.fibre.StdScatMag*randn(nbScats,1)));%Table of power attenuation ('am') per scatterer, linear. No fibre attenuation added. 

scatDistTab = p.fibre.scatDistTab; %z_m
scatMagTab = p.fibre.scatMagTab; %a_m

%% 2)Distribution of the scatterers between the fibre segments
xMaxTab = single(p.fibre.TxSpatialRes*(1:p.fibre.nbSegments)); %Table of starting positions of the fibre segments

%find zm and am of the scatterers in each segment
for n=1: p.fibre.nbSegments
   tmpIdxTab = find((scatDistTab>xMaxTab(n)-p.fibre.TxSpatialRes) & scatDistTab<=xMaxTab(n));
   nbScatsPerSegment = length(tmpIdxTab);
   r_(1:nbScatsPerSegment,n) = xMaxTab(n)-scatDistTab(tmpIdxTab); % Distances of scatterers from the current segment start
   a_(1:nbScatsPerSegment,n) = scatMagTab(tmpIdxTab); % Magnitudes of scatterers in the current segment 
end

% Handle case where scatterer positions for one segment change in time (local dynamic model). 
% It will allow for local recalculation of optical field as fct of time
if p.fibre.ExcitedSegmentFlag==1 %If a segment is excited, store the magnitude & initial distance (static case) of its scatterers
     p.fibre.ExcitedSegment_xMax = xMaxTab(p.fibre.ExcitedSegmentIdx); %End of segment coordinate
     p.fibre.ExcitedSegmentNbScats = length(find((scatDistTab>p.fibre.ExcitedSegment_xMax-p.fibre.TxSpatialRes) & scatDistTab<=p.fibre.ExcitedSegment_xMax));
     p.fibre.ExcitedSegmentScatDistInit = r_(1:p.fibre.ExcitedSegmentNbScats,p.fibre.ExcitedSegmentIdx);%Initial distances of scatterers for the excited segment
     a_(1:p.fibre.ExcitedSegmentNbScats,p.fibre.ExcitedSegmentIdx) = p.fibre.artificialFading*a_(1:p.fibre.ExcitedSegmentNbScats,p.fibre.ExcitedSegmentIdx);  % add fading to perturbation location
     p.fibre.ExcitedSegmentScatMag = a_(1:p.fibre.ExcitedSegmentNbScats,p.fibre.ExcitedSegmentIdx);%Magnitudes of scatterers for the excited segment
 end

% Calculation of the backscattered optical field on two polarizations
p.fibre.Ei = single(sum(a_.*exp(1j*4*pi*p.N/p.tx.Lambda*r_)) .* exp(1j*4*pi*p.N/p.tx.Lambda*xMaxTab));%Backscattered optical field per segment
p.fibre.Ai = single(10.^(1e-3*p.fibre.LossdB*xMaxTab/10));% fibre loss (round-trip) - x0.5 factor since optical field and x2 since round-trip
p.fibre.Ei(1) = 0.01; %put big reflection at beginning to make sure we detect the beginning of the fiber
p.fibre.Ei(end) = 0.01; %put big reflection at end to make sure we detect the end of the fiber
p.fibre.Ai(1) = 3; %put big reflection at beginning to make sure we detect the beginning of the fiber
p.fibre.Ai(end) = 3; %put big reflection at end to make sure we detect the end of the fiber

%% ADDED Generate Jones Matrices before dynamic update
[Hi, p] = genJonesMatrices(p); % Compute Jones Matrices

if ~exist('gCodeSingle', 'var') || isempty(gCodeSingle)
    error('gCodeSingle is undefined. Ensure it is passed to genRayleighScattering.');
end

HiRep = repmat([Hi, zeros(2, 2 * (size(gCodeSingle, 2) - p.fibre.nbSegments))], 1, p.tx.nbCodes);

%% ADDED Handle dynamic fiber model update
if p.fibre.ExcitedSegmentFlag == 1
    [p, HiRep] = genDynamicSegt(p, HiRep, p.fibre.ExcitedSegmentIdx(1), ...
        p.fibre.ExcitedStrainMax(1), p.fibre.ExcitedF_event(1), ...
        p.fibre.ExcitedDynEvolution(1));
end

end
