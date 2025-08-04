function [r,p,gCodeRep] = lasernoise_model(gCodeRep,p,r)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASERNOISE_MODEL adds phase noise to simulated laser input. 
% White frequency noise model (Lorentzian lineshape)
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

codelength = size(gCodeRep,2); 

pn_x = sqrt(4*pi*(p.tx.dfLaser/2)/p.tx.fSymb).*cumsum(randn(1,codelength));%(dfLaser/2) is the std dev of the laser freq. bias
noiseSeq = exp(1i*pn_x);%The phase noise vector (Lorentzian model)
r.laserNoiseMat = repmat(noiseSeq,2,1);%Duplicate phase noise vector for the 2 polars

 gCodeRep = gCodeRep.*r.laserNoiseMat(:,1:codelength ); %Apply laser phase noise to the transmitted signal (same noise for both polar)
figure;
subplot(2,1,1);
plot(angle(noiseSeq)); title('Phase Noise Over Time'); xlabel('Samples'); ylabel('Phase (rad)');
subplot(2,1,2);
plot(abs(fft(noiseSeq))); title('Power Spectrum of Phase Noise'); xlabel('Frequency (Hz)'); ylabel('Magnitude');



 %% To test: betaline and 1/f model
%  f_betaline = 3.8e3; %Hz, see Di domenico 2010 : is frequency up to which contribution on spectrum is 1/f not white
%  
%  t_probing = p.tx.nbCodes * p.tx.Tcode; % probing duration
%  f_obs = 2/t_probing; %minimum frequency we can capture
%  f_code = 1/p.tx.Tcode; %max frequency captured
%  
%  if nu ==1
%      oneoverf_lw = f_betaline*8*log(2)*sqrt(log(f_betaline*t_probing))/pi^2; %sqrt(8*log(2)*sum(Sv_f,2));% Hz % see PhD Omar Sahni p.20
%  elseif nu > 1
%      oneoverf_lw = f_betaline*8*log(2)*sqrt(((f_betaline*t_probing)^(nu-1)-1)/(nu-1)); %DiDomenico2010
%  end
%  oneoverf_level = sqrt(pi*oneoverf_lw/p.tx.fSymb); % same as for noise_level (shouldnt be)
%  
%  oneoverf_noise = oneoverf_level.*cumsum(randn(1,length(gCodeRep))); %1/f noise contribution to linewidth
%  
%  fprintf('Low frequencies min : %1.2d Hz, fq max %2.2f Hz \n', f_obs, p.tx.fSymb);
