function Erx = conv_fft(Erx, h11,h12,h21,h22)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution (spectral domain implementation)
% Computes the convolutions between Erx(1,:), Erx(2,:) and the channel
% impulse response h11, h12, h21, h22
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Erx = double(Erx);
%disp(size(Erx));
%disp(size(h11));
Ly=length(Erx)+length(h11)-1;    % Output size for each of the 4 convolutions
Ly2=2^(nextpow2(Ly));            % Find smallest power of 2 that is > Ly
% Apply windowing to the impulse responses to reduce temporal leakage
%h11 = h11 .* hamming(length(h11))'.^2;
%h12 = h12 .* hamming(length(h12))'.^2;
%h21 = h21 .* hamming(length(h21))'.^2; 
%h22 = h22 .* hamming(length(h22))'.^2;
%h11 = h11 .* hann(length(h11))'.^2;
%h12 = h12 .* hann(length(h12))'.^2;
%h21 = h21 .* hann(length(h21))'.^2; 
%h22 = h22 .* hann(length(h22))'.^2;

%disp(size(Erx(1,:)));
X=fft(Erx(1,:), Ly2);                 % Fast Fourier transformfor for Erx(1,:)
Y11=ifft(fft(h11, Ly2).*X, Ly2);      % Inverse fast Fourier transform of spectral products
Y21=ifft(fft(h21, Ly2).*X, Ly2);

X=fft(Erx(2,:), Ly2);                 % Fast Fourier transformfor for Erx(2,:)
Y12=ifft(fft(h12, Ly2).*X, Ly2);      % Inverse fast Fourier transform of spectral products
Y22=ifft(fft(h22, Ly2).*X, Ly2);

%disp('First 10 samples of the 4 convolutions:');
%disp([Y11(1:10); Y12(1:10); Y21(1:10); Y22(1:10)]);
Erx = ([Y11(1:Ly)+Y12(1:Ly); Y21(1:Ly)+Y22(1:Ly)]);
%disp(size(Erx));
end