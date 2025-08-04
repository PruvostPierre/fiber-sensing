function conv_fft(h11,h12,h21,h22)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution (spectral domain implementation)
% Computes the convolutions between Erx(1,:), Erx(2,:) and the channel
% impulse response h11, h12, h21, h22
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Erx;
Ly=length(Erx)+length(h11)-1;    % Output size for each of the 4 convolutions
Ly2=2^(nextpow2(Ly));            % Find smallest power of 2 that is > Ly

X=fft(Erx(1,:), Ly2);                 % Fast Fourier transformfor for Erx(1,:)
Y11=ifft(fft(h11, Ly2).*X, Ly2);      % Inverse fast Fourier transform of spectral products
Y21=ifft(fft(h21, Ly2).*X, Ly2);

X=fft(Erx(2,:), Ly2);                 % Fast Fourier transformfor for Erx(2,:)
Y12=ifft(fft(h12, Ly2).*X, Ly2);      % Inverse fast Fourier transform of spectral products
Y22=ifft(fft(h22, Ly2).*X, Ly2);

Erx = ([Y11(1:Ly)+Y12(1:Ly); Y21(1:Ly)+Y22(1:Ly)]);
end