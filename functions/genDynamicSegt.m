function [p,HiRep] = genDynamicSegt(p,HiRep,ExcitedSegmentIdx,ExcitedStrainMax,ExcitedF_event,ExcitedDynEvolution)
%CHANGED
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that transforms the selected static fiber segment in a dynamic segment by locally applying a predefined perturbation
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%hMax = (p.fibre.spatialRes) * 0.5 * sqrt((1+ExcitedStrainMax)^2-1^2);%Maximal transversal deformation in the
% middle of the segment: EA 190224 not required

hMax = p.fibre.spatialRes*ExcitedStrainMax;%Maximal axial deformation in the % middle of the segment
%hInstantTab = hMax*sin(2*pi*ExcitedF_event*(1 : p.tx.nbCodes-1)*p.tx.Tcode);%Table of Instantaneous transversal deformations (sine wave) over the acquisition period
hInstantTab = hMax*sin(2*pi*ExcitedF_event*(1 : p.tx.nbCodes-1)*p.tx.Tcode);% + 0.3*hMax*cos(4*pi*ExcitedF_event*(1 : p.tx.nbCodes-1)*p.tx.Tcode);


if ExcitedDynEvolution==1 %Linear increase of perturbation
        hInstantTab = hInstantTab.*((1/p.tx.nbCodes)*(1 : p.tx.nbCodes-1));%Instantaneous transversal deformation
elseif ExcitedDynEvolution==-1 %Linear increase of perturbation
        hInstantTab = hInstantTab.*((1/p.tx.nbCodes)*(p.tx.nbCodes-1:-1:1));%Instantaneous transversal deformation
elseif ExcitedDynEvolution~=0
    fprintf('Error in genDynamicSegt(). Wrong Value entered for ExcitedDynEvolution. Exit'); exit();
end

scatDistExtensionTab = (p.fibre.spatialRes + hInstantTab)./p.fibre.spatialRes; %Table of instantaneous fiber extension factors induced by deformation at excited segment position: axial strain

figure(110); hold off;
plot(scatDistExtensionTab); hold on; xlabel('Time (Nb codes)'); title('scatterer distance extension factor at dynamic segment position');

EiDynCorrection= NaN(1,p.tx.nbCodes-1); %Table to store (for display purpose) the correction factor applied to the excited fiber segment along the time axis
expSumReflectorsRef = sum(double(p.fibre.ExcitedSegmentScatMag).*exp(1j*4*pi*p.N*double(p.fibre.ExcitedSegmentScatDistInit)/p.tx.Lambda));%Sum of phasors for static case (ref) %p.fibre.Ei(p.fibre.ExcitedSegmentIdx); %
%CHANGED
for nn=1 : p.tx.nbCodes-1 %Loop that recalculates the Jones matrix of the dynamic segment for each transmitted code (apart from the first ref one: static case) and stores it at the right place in HiRep
    %expSumReflectorsInstant = sum(double(p.fibre.ExcitedSegmentScatMag).*exp(1j*4*pi*p.N*double(p.fibre.ExcitedSegmentScatDistInit).*scatDistExtensionTab(nn)./p.tx.Lambda ));%Sum of phasors for dyn case at current time code index
    expSumReflectorsInstant = sum(double(p.fibre.ExcitedSegmentScatMag).*exp(1j*4*pi*p.N*double(p.fibre.ExcitedSegmentScatDistInit+hInstantTab(nn))./p.tx.Lambda ));%Sum of phasors for dyn case at current time code index
    
    %EiDynCorrection(nn) = expSumReflectorsInstant/expSumReflectorsRef; %The correction factor to apply at the current time index to the excited segment Jones matrix Hi
    EiDynCorrection(nn) = exp(1j*4*pi*p.N*double(hInstantTab(nn))./p.tx.Lambda);
    
    offsetExcitation = 2*(ExcitedSegmentIdx + nn*p.tx.Ncode*p.tx.ovsFactor) -1; %Current dynamic segment matrix position in HiRep    
    HiRep(:,offsetExcitation:offsetExcitation+1) = EiDynCorrection(nn)*HiRep(:,offsetExcitation:offsetExcitation+1) ;%Update phasor and intensity only (polar params assumed unchanged for now)
       
    mmTab = 2*nn*p.tx.Ncode + (2*(p.fibre.ExcitedSegmentIdx+1) - 1 : 2*p.fibre.nbSegments);%Table of all Hi indices to update (all fiber segments that follow the perturbation)
    HiRep(:,mmTab) = EiDynCorrection(nn) + HiRep(:,mmTab);
    %HiRep(:,mmTab) = (EiDynCorrection(nn)/abs(EiDynCorrection(nn)))*HiRep(:,mmTab); % CHECK: update phase changes in following segments
end

%% 
figure(131);
plot(unwrap(angle(HiRep(1,:))));
title('Phase evolution of HiRep');

figure(120); hold off;
plot(unwrap(angle(EiDynCorrection)));
hold on; xlabel('Time (Nb codes)'); 
title('angle of correction factor at dynamic segment position');

figure(121); hold off;
plot(abs(EiDynCorrection)); 
hold on; xlabel('Time (Nb codes)'); 
title('module of correction factor at dynamic segment position');



end

