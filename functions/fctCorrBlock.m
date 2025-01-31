function out = fctCorrBlock(Y, gCode_fft, Hlen,offset_ratio)
%CHANGED
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation process
% Computes per block the 4 cross-correlations between x(1,:), x(2,:) and h(1,:), h(2,:) 
% The input block length is already a power of 2. 
% The input block Y is a time-domain signal
% gCode_fft is the precalculated spectrum of the code, same size as the input block. 
% 1)xSeg is reversed prior to the fct call to save time,
% 2)output 2D table is preallocated in the upper fct 
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1=conj(fft(Y(:,1))); % conjugate of Fast Fourier transform for x(1,:)
X2=conj(fft(Y(:,2))); % conjugate of Fast Fourier transform for x(2,:)

% ifft of spectral products 
y11=ifft(X1.*gCode_fft(:,1)); 
y21=ifft(X2.*gCode_fft(:,1));
y12=ifft(X1.*gCode_fft(:,2));
y22=ifft(X2.*gCode_fft(:,2));     

% We need to flip the correlation result and shift by one sample to get the
% correct result
y11 = circshift(flipud(y11),1);
y12 = circshift(flipud(y12),1);
y21 = circshift(flipud(y21),1);
y22 = circshift(flipud(y22),1);
%CHANGED
normalize = true;
if normalize
    scaleFactor = Hlen;
else
    scaleFactor = 1;                                                                              
end

% Output with adjusted scaling
out(1,:) = y11(round(offset_ratio*Hlen)+1:end) ./ scaleFactor;
out(2,:) = y21(round(offset_ratio*Hlen)+1:end) ./ scaleFactor;
out(3,:) = y12(round(offset_ratio*Hlen)+1:end) ./ scaleFactor;
out(4,:) = y22(round(offset_ratio*Hlen)+1:end) ./ scaleFactor;

