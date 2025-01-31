function out = getAnglePlusMinusPiOver2( Phi )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: an angle in radian expressed within interval[-inf;+inf]
%Output: the same angle in radian but expressed within interval[-pi/2;+pi/2]
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


out = sign(Phi).*mod(abs(Phi),pi);%Constraint Phi elts between -pi/2 and +pi/2
out = (out>pi/2)*(-pi) + out;
out = (out<-pi/2)*(pi) + out;


