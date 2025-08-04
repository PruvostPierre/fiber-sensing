%% MecaFunctions.m - Functions to show the curves

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MecaFunctions.m - Classe contenant les fonctions de méca
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef MecaFunctions

   methods(Static)
        
    function microstrain = StrainfromPHI(DPhi, LamdaLaser,indRef,CphotoElas,Lchannel)
        
        microstrain = (DPhi * LamdaLaser*10^-9)/(indRef*CphotoElas*Lchannel*2*pi)*10^6;
    end

          
    function StressScalOutput = StressScal(E, strain, rowStudied)
        StressScalOutput = E * strain(rowStudied, :);
    end 

   end   
end

 