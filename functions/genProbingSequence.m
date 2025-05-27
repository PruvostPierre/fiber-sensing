function out = fctCorrBlock(Y, gCode_fft, codeLen,offset_ratio)
%CHANGED
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation process
% Computes per block the 4 cross-correlations between Y(1,:), Y(2,:) and gCode_fft(1,:), gCode_fft(2,:) 
% The input block length is already a power of 2. 
% The input block Y is a time-domain signal
% gCode_fft is the precalculated spectrum of the code, same size as the input block. 
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the convolution theorem to implement correlation in Fourier domain
% TF ( x(t) * y'(-t)) = TF(x).TF(y)' where (.)' stands for complex
% conjugate and * stands for convolution operator

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

normalize = true; % FIX ME: replace it by a parameter in p.rx? - EA 27/05/25
if normalize
    scaleFactor = codeLen;
else
    scaleFactor = 1;                                                                              
end

% Output with adjusted scaling and first codeLen-long part discarded due to
% possible time aliasing
out(1,:) = y11(round(offset_ratio*codeLen):end) ./ scaleFactor;
out(2,:) = y21(round(offset_ratio*codeLen):end) ./ scaleFactor;
out(3,:) = y12(round(offset_ratio*codeLen):end) ./ scaleFactor;
out(4,:) = y22(round(offset_ratio*codeLen):end) ./ scaleFactor;



% FIXME: test and compare this alternative
% % CD Fast complex correlation derived from FCONV Fast Convolution (spectral domain implementation), called by fcorrCD_PerBlock_4X()
% % Computes per block the 4 cross-correlations between x(1,:), x(2,:) and h(1,:), h(2,:) 
% %The input block length is already power of 2. spHseg is the precalculated spectrum of h, same size as the input block. 
% 
% X1=fft(conj(xSeg(end:-1:1,1))); % Fast Fourier transform for x(1,:)
% X2=fft(conj(xSeg(end:-1:1,2))); % Fast Fourier transform for x(2,:)
% 
% y11=ifft(X1.*spHseg(:,1)); % Inverse fast Fourier transform of spectral products, no zero padding needed here since already a power of 2 
% y21=ifft(X2.*spHseg(:,1));
% y12=ifft(X1.*spHseg(:,2));
% y22=ifft(X2.*spHseg(:,2));     
% 
% out(:,1) = conj(y11(end:-1:Hlen))./Hlen; % Take just the last (length(xSeg)-Hlen+1) elements, in reverse order, and then normalize
% out(:,2) = conj(y21(end:-1:Hlen))./Hlen;
% out(:,3) = conj(y12(end:-1:Hlen))./Hlen;
% out(:,4) = conj(y22(end:-1:Hlen))./Hlen;
