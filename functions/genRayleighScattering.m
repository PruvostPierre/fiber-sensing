function p = genRayleighScattering(p)

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

scatDistTab = p.fibre.scatDistTab;
scatMagTab = p.fibre.scatMagTab;

%% 2)Distribution of the scatterers between the fibre segments
xMaxTab = single(p.fibre.TxSpatialRes*(1:p.fibre.nbSegments)); %Table of starting positions of the fibre segments

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

end