function p=fcorrAndGetJones(gCode, brutBlockLen, p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation with codes and Jones matrices extraction. 
% Subblock partitioning for correlation and Jones matrix extractions 
% The subblock must contain a minimum number of time codes (say 10) 
% to ensure a reliable Rayleigh detection from the first subblock 
% and then the extraction of the Jones matrices inside this 1st subblock.
% Selected Jones matrices HiTab are now returned in struct p
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Erx;
disp(size(Erx))
p.rx.ErxLen = brutBlockLen; 
Erx = Erx(1:brutBlockLen,:); 

hLen  = length(gCode)*p.rx.ovsFactor/p.tx.ovsFactor;

% Process 1st block apart to include time synchro (Rayleigh start & stop positions) prior to Jones matrix extraction
offsetInStart = 1; blockIdx = 1;

%gCode = repmat(gCode, 1, p.tx.nbCodes);
% Code spectrum is zero-padded to reach length 'brutBlockLen'
gCode_fft(:,1) = fft(gCode(1,:).',brutBlockLen);
gCode_fft(:,2) = fft(gCode(2,:).',brutBlockLen);

tmpCorr = fctCorrBlock(Erx(offsetInStart:1:offsetInStart+brutBlockLen-1,:), gCode_fft, hLen,p.rx.offset_ratio);

%for i=1:p.tx.nbCodes
%    if i==1
%        tmpCorr = fctCorrBlock(Erx(offsetInStart:1:length(gCode),:), gCode_fft, hLen,p.rx.offset_ratio);
%    else
%        tmpCorr = [tmpCorr, fctCorrBlock(Erx(offsetInStart+(i-1)*length(gCode):1:i*length(gCode),:), gCode_fft, hLen,p.rx.offset_ratio)];
%    end
%end
netBlockLen = length(tmpCorr);
nbFullBlocks = floor((p.rx.ErxLen-(hLen-1))/netBlockLen);
nbRemainder = mod(p.rx.ErxLen-(hLen-1),netBlockLen);

%% Time synchronization: search for Rayleigh Start & Stop from the first processed subblock
ovsCodeLength = p.rx.ovsFactor*p.tx.Ncode;%nb of samples between 2 consecutive codes in the current trace

nbExpectedCodes = floor(length(tmpCorr)/ovsCodeLength);%nbExpectedCodes in the 1st subblock

%sumAbsTmpCorr = tmpCorr(1,:).*conj(tmpCorr(1,:))+tmpCorr(2,:).*conj(tmpCorr(2,:))+tmpCorr(3,:).*conj(tmpCorr(3,:))+tmpCorr(4,:).*conj(tmpCorr(4,:));%Frobenius norm 
sumAbsTmpCorr = abs(tmpCorr(1,:).*tmpCorr(4,:)+tmpCorr(2,:).*tmpCorr(3,:));% |det| criterion

if p.tx.detectionindex(1) == 0 && p.tx.detectionindex(2) == 0 %no detection conducted yet
    [xStartRay,xStopRay,estimatedFiberLength,p] = getRayleighSpread(sumAbsTmpCorr, ovsCodeLength, nbExpectedCodes, p); %Get Start & Stop positions of Rayleigh backscatter, processed here from the 1st block of data
    p.tx.detectionindex = [xStartRay,xStopRay] ;
    
else
    estimatedFiberLength = p.fibre.L;
    xStartRay = p.tx.detectionindex(1); xStopRay = p.tx.detectionindex(2);
end

p.fibre.L = estimatedFiberLength;

if nbRemainder<hLen %Useless to carry out the remaining truncated block if its size < hLen %FIXME not tested
    p.HiTab = single(NaN(4,nbFullBlocks*p.rx.nbOvsReflectors*(nbExpectedCodes+1)));%CD: prealloc for speed issue, with an overevaluated table size
else
    p.HiTab = single(NaN(4,(nbFullBlocks+1)*p.rx.nbOvsReflectors*(nbExpectedCodes+1)));%CD: prealloc for speed issue, with an overevaluated table size
end
rayleighBlockMask = getRayleighBlockMask(netBlockLen, blockIdx, xStartRay, xStopRay, ovsCodeLength);%Get the mask that localizes the series of Rayleigh Jones matrices contained in the current subblock 'tmpCorr'

idxStart = 1; nbMaskElts = sum(rayleighBlockMask==1);
p.HiTab(:,idxStart:idxStart+nbMaskElts-1) = tmpCorr(:,rayleighBlockMask==1);%Extract all the Rayleigh Jones matrices from the current block and store them in p.HiTab

%fprintf('\nSpatial resolution (ovs): %.2fm   Estimated fibre length:%.1fm   Laser coherence length:%.0fkm  code length (SB%d):%.0fkm\n', 0.5*p.fibre.cFiber/p.rx.fSamp, estimatedFiberLength, 1.e-3*p.fibre.cFiber/(pi*p.tx.dfLaser), p.tx.seqOrderCst, 1.e-3*p.fibre.cFiber*p.tx.Tcode);
fprintf('\nSpatial resolution (ovs): %.2fm   Estimated fibre length:%.1fm   Laser coherence length:%.0fkm  code length (SB%d):%.0fkm\n', 0.5*p.fibre.cFiber/(2*p.tx.fSymb), estimatedFiberLength, 1.e-3*p.fibre.cFiber/(pi*p.tx.dfLaser), p.tx.seqOrderCst, 1.e-3*p.fibre.cFiber*p.tx.Tcode);

%% Process the next blocks to iteratively extract their Jones matrices
for blockIdx=2 : nbFullBlocks %FIXME not tested yet
    offsetInStart = offsetInStart + netBlockLen; %input block translation is netBlockLen
% B2B deconvolution prior to correlation
%     if p.file.DECONV_B2B==1 ||  p.file.DECONV_B2B==2 %the 2 received channels are deconvolved with the same filter hB2B
%         Erx(offsetInStart:offsetInStart+brutBlockLen-1,1) = ifft( fft(Erx(offsetInStart:offsetInStart+brutBlockLen-1,1))./sp_hB2B);
%         Erx(offsetInStart:offsetInStart+brutBlockLen-1,2) = ifft( fft(Erx(offsetInStart:offsetInStart+brutBlockLen-1,2))./sp_hB2B);
%     elseif p.file.DECONV_B2B==3 %the 2 received polar channels are deconvolved with  dedicated filters hB2B_1 and hB2B_2 resp.
%         Erx(offsetInStart:offsetInStart+brutBlockLen-1,1) = ifft( fft(Erx(offsetInStart:offsetInStart+brutBlockLen-1,1))./sp_hB2B_1 );
%         Erx(offsetInStart:offsetInStart+brutBlockLen-1,2) = ifft( fft(Erx(offsetInStart:offsetInStart+brutBlockLen-1,2))./sp_hB2B_2 );
%     end
       
    %%
    tmpCorr = fctCorrBlock(Erx(offsetInStart+brutBlockLen-1:-1:offsetInStart,:), gCode_fft, hLen,p.rx.offset_ratio);
    rayleighBlockMask = getRayleighBlockMask(netBlockLen, blockIdx, xStartRay, xStopRay, ovsCodeLength);%Get the mask that localizes the series of Rayleigh Jones matrices contained in the current subblock 'tmpCorr'
    idxStart = idxStart+nbMaskElts; nbMaskElts = sum(rayleighBlockMask==1);
    p.HiTab(:,idxStart:idxStart+nbMaskElts-1) = tmpCorr(:,rayleighBlockMask==1);%Extract all the Rayleigh Jones matrices from the current block and store them in p.HiTab
end

%Process the remaining block (size<netBlockLen)
if nbRemainder>=hLen %Useless to carry out the remaining truncated block if its size < hLen %FIXME not tested

    gCode_fft(:,1) = fft(gCode(1,:).',nbRemainder+hLen-1);%code spectrum zero-padded to reach length of processed block
    gCode_fft(:,2) = fft(gCode(2,:).',nbRemainder+hLen-1);
    offsetInStart = offsetInStart + netBlockLen; %input block translation is netBlockLen
    clear tmpCorr; tmpCorr = single(zeros(4,nbRemainder));
    
%%%%
% % for B2B deconvolution prior to correlation on the remaining received data block
% if p.file.DECONV_B2B==1 ||  p.file.DECONV_B2B==2 %the 2 received channels are deconvolved with the same filter hB2B
%         sp_hB2B = fft(hB2B, nbRemainder);%Spectral response after zero padding to match size of the received data block onto which the deconvolution will be applied
%        Erx(offsetInStart:offsetInStart+nbRemainder-1,1) = ifft( fft(Erx(offsetInStart:offsetInStart+nbRemainder-1,1))./sp_hB2B );
%        Erx(offsetInStart:offsetInStart+nbRemainder-1,2) = ifft( fft(Erx(offsetInStart:offsetInStart+nbRemainder-1,2))./sp_hB2B );
% elseif p.file.DECONV_B2B==3 %the 2 received polar channels are deconvolved with  dedicated filters hB2B_1 and hB2B_2 resp.
%         sp_hB2B_1 = fft(hB2B_R1, nbRemainder);%Spectral response after zero padding to match size of the received data block onto which the deconvolution will be applied
%         sp_hB2B_2 = fft(hB2B_R2, nbRemainder);%Spectral response after zero padding to match size of the received data block onto which the deconvolution will be applied
%        Erx(offsetInStart:offsetInStart+nbRemainder-1,1) = ifft( fft(Erx(offsetInStart:offsetInStart+nbRemainder-1,1))./sp_hB2B_1 );
%        Erx(offsetInStart:offsetInStart+nbRemainder-1,2) = ifft( fft(Erx(offsetInStart:offsetInStart+nbRemainder-1,2))./sp_hB2B_2 );
% end  

    tmpCorr = fctCorrBlockRev(Erx(offsetInStart+nbRemainder+hLen-1-1:-1:offsetInStart,:), gCode_fft, hLen);
    blockIdx = nbFullBlocks + 1;
    rayleighBlockMask = getRayleighBlockMask(netBlockLen, blockIdx, xStartRay, xStopRay, ovsCodeLength);%Get the mask that localizes the series of Rayleigh Jones matrices contained in the current subblock 'tmpCorr'
    idxStart = idxStart+nbMaskElts; nbMaskElts = sum(rayleighBlockMask(1:nbRemainder)==1);
    p.HiTab(:,idxStart:idxStart+nbMaskElts-1) = tmpCorr(:,rayleighBlockMask(1:nbRemainder)==1);%Extract all the Rayleigh Jones matrices from the current block and store them in p.HiTab
    p.HiTab = p.HiTab(:,1:idxStart+nbMaskElts-1); %remove the zeros
else
    p.HiTab = p.HiTab(:,1:idxStart+nbMaskElts-1);
end

if xStartRay > xStopRay   
    p.HiTab = circshift(p.HiTab, -xStopRay, 2); 
end
