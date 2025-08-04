%% main.m - main of the meca test

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing of experimental measurements for a coded DAS system
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Addpath to the code
addpath('D:\01-PHD_TelecomParis\01-TelecomParis_1ère_année_2025-2026\13-Git\fiber-sensing\functions')
addpath('D:\01-PHD_TelecomParis\01-TelecomParis_1ère_année_2025-2026\13-Git\fiber-sensing\.git')
addpath('D:\01-PHD_TelecomParis\01-TelecomParis_1ère_année_2025-2026\13-Git\fiber-sensing\Meca')
%%
TestInProgress = initialization();
disp(TestInProgress.Config);

Display.quickdisplay(TestInProgress)

microstrainTAB = MecaFunctions.StrainfromPHI(TestInProgress.DAS.diffPhiTabSelect, TestInProgress.Lambda, TestInProgress.n, TestInProgress.CphotoElas, TestInProgress.Lchannel);
strainTAB = microstrainTAB/10^6;
microstrainMeanCol = sum(microstrainTAB, 1);

Display.DisplayComp(TestInProgress.IMC.LoadCell50KN_time,TestInProgress.IMC.LoadCell50KN,microstrainMeanCol,'LoadCell','microstrain cumulée sur la fibre')


TestInProgress.C_SEG=Coherency.CohSeg(strainTAB,0.1,TestInProgress);
smooth_C_SEG=smooth(TestInProgress.C_SEG);
smooth_C_SEG50=smooth(TestInProgress.C_SEG,50);

Display.PLOT_Lin(TestInProgress.C_SEG,'Cohenrency raw','on','channels', 'cohenrency')
Display.PLOT_Lin(smooth_C_SEG,'Cohenrency with smooth 1','on','channels', 'cohenrency')
Display.PLOT_Lin(smooth_C_SEG50,'Cohenrency with smooth 50','on','channels', 'cohenrency')

StressTAB = MecaFunctions.StressScal(60000,strainTAB,TestInProgress.rowStudied);
Display.PLOT_Lin(StressTAB,'Stress Optical Fiber','on','sample', 'stress amplitude Pa')

% Compute cross-correlation

rowStudied = TestInProgress.rowStudied;
selectedRow = TestInProgress.DAS.diffPhiTabSelect(rowStudied, :);
LenDAS = length(selectedRow);
LenMeca = length(TestInProgress.IMC.LoadCell50KN);

% Interpolation
x_original = linspace(1, LenDAS, LenDAS);
x_target = linspace(1, LenDAS, LenMeca);
selectedRow_resized = interp1(x_original, selectedRow, x_target, 'linear');



A = xcorr(selectedRow_resized, TestInProgress.IMC.Acc_Mono_Z,'normalized');
[B, C] = Coherency.CircCorr(selectedRow_resized, TestInProgress.IMC.Acc_Mono_Z);

AA = xcorr(selectedRow_resized, selectedRow_resized,'normalized');


Display.PLOT_Lin(A,'Corrélation Das Row studied and Acc Mono Z','on','samples','Correlation coeff')
Display.PLOT_Lin(B,'Corrélation Das Row with and Acc Mono Z circulary shifted','on','samples','Correlation coeff')
Display.PLOT_Lin(C,'Max xcorr coef','on','samples','Correlation coeff')


Display.PLOT_Lin(AA,'Autocorrélation Das studied and Acc Mono Z','on','samples','Correlation coeff')





t1 = 1;
t2 = 1.01;

s1=t1*length(selectedRow_resized)/TestInProgress.timeRecord
s2=t2*length(selectedRow_resized)/TestInProgress.timeRecord

% Assuming FFTDAS is your input signal
segment = FFTDAS(s1:s2); % Extract the segment from sample 20000 to 42000

% Compute the FFT of the segment
Y = fft(segment);

% Optional: Plot the magnitude of the FFT
Fs = 200; % Example sampling frequency (adjust as necessary)
L = length(segment);
f = Fs*(0:(L/2))/L; % Frequency axis

% Compute the two-sided spectrum
P2 = abs(Y/L);

% Compute the single-sided spectrum
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Plot the single-sided amplitude spectrum
figure;
loglog(f, P1);
title('Single-Sided Amplitude Spectrum');
xlabel('Frequency (f)');
ylabel('|P1(f)|');


DSP_Freq.FFT_Wiwdow(selectedRow_resized,1,1.1,200,TestInProgress.timeRecord,'LogLog');