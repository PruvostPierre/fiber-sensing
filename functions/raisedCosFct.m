function sigOut = raisedCosFct(sigIn, edgeRatio)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smooth the start & stop edges of vector 'sigIn' with a raised cosine .
%The smoothing width is fixed by param 'edgeRatio'
%(0.1=>start edge is 10% of the signal length and the same for stop edge)
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

incr = pi/round(edgeRatio*size(sigIn,2));
raisedStart = 0.5*(1+cos(-pi+incr:incr:-incr)); raisedStop = raisedStart(end:-1:1);
winRaised = ([raisedStart  ones(1, size(sigIn,2)-2*length(raisedStart)) raisedStop]);

sigOut=sigIn.*winRaised;
