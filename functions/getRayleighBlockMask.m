function [RayleighBlockMask] = getRayleighBlockMask(blockLen, blockIndex, xStartRay, xStopRay, codeLen)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returning a mask of the useful part of the Rayleigh backscattering
% indices for the current block.
% We select the useful part of the impulse response corresponding to the
% actual fiber by detecting first and last reflectors.
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nStart = ceil( ((blockIndex-1)*blockLen - xStartRay)/codeLen );%nStart: an
%integer
nStart = ceil( ((blockIndex-1)*blockLen - (xStartRay-1))/codeLen );
xStartRayCurrentBlock = xStartRay + nStart*codeLen - (blockIndex-1)*blockLen; %Take the first Rayleigh start position in the current block

%nStop = ceil( ((blockIndex-1)*blockLen - xStopRay)/codeLen );
nStop = ceil( ((blockIndex-1)*blockLen - (xStopRay-1))/codeLen );
xStopRayCurrentBlock = xStopRay + nStop*codeLen - (blockIndex-1)*blockLen; %Take the first Rayleigh stop position in the current block

codeRayleighMask = int8(zeros(1,codeLen));%A mask having the length of the code to store the Rayleigh backscatter zone 
if xStartRayCurrentBlock<=xStopRayCurrentBlock
    codeRayleighMask(xStartRayCurrentBlock:xStopRayCurrentBlock) = 1;
else
    codeRayleighMask(1:xStopRayCurrentBlock) = 1; codeRayleighMask(xStartRayCurrentBlock:end) = 1;
end

RayleighBlockMask = repmat(codeRayleighMask, 1, floor(blockLen/codeLen)); %Duplicate the mask per code along the block
RayleighBlockMask = ([RayleighBlockMask codeRayleighMask(1:mod(blockLen,codeLen))]);
