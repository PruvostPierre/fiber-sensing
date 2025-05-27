function [p,r]=Fileread_and_correlation(p,r)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates the probing sequence, models the sensed fiber, 
% the transmitter and the receiver frontends, and performs the correlation
% process and the Jones matrices extraction at the RX side.
%
% Authors: 
% Original code by E. Awwad, P. Pruvost, D. Prato - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Erx; % global variable containing the propagated field
tStart = tic;
%% Start timing for execution
rx_start = tic; % Start timer for the entire process
tCumul = 0;   % Initialize cumulative time counter


%% Reading acquired files
if (p.rx.data_filename(end-2:end)=='dat') 
    % GageStream2Disk saves all channels in one file in a time division
    % multiplexed way, for instance if we call the four channels
    % a1,a2,a3,a4, the samples in the single-column file are a1(0), a2(0),
    % a3(0),a4(0),a1(1), a2(1),a3(1),a4(1),etc.
    fileID = fopen([p.rx.data_directory,p.rx.data_filename]);
    B = fread(fileID,p.rx.samples,'int16');

    % rescaling all channels to values in Volts
    B = B.*p.rx.acq_card_range./2^16;

    A(:,1)=B(1:4:end);
    A(:,2)=B(2:4:end);
    A(:,3)=B(3:4:end);
    A(:,4)=B(4:4:end);
    fclose(fileID);
    clear B;
else
  for i = 1:4
        A(:,i) = readmatrix([p.rx.data_directory,p.rx.data_filename,'_CH',num2str(i),'.txt'],'NumHeaderLines',12);
  end
end
Erx(1,:)= A(:,1)+1i*A(:,2); % Associate channels to PolX I and Q
Erx(2,:)= A(:,3)+1i*A(:,4); % Associate channels to PolY I and Q
clear A;

% Estimated number of acquired codes
p.tx.nbCodes = round(length(Erx)./length(p.rx.gCode));

    % Blocking DC component (to emulate DC block of each of the 4 balanced PD outputs)
    if p.rx.RxDcBlocker~=0
        Erx(1,:) = Erx(1,:) - mean(Erx(1,:),2);
        Erx(2,:) = Erx(2,:) - mean(Erx(2,:),2);
    end
    
    if p.displ.powerIndication
        P_erx_moy = mean(mean(abs(Erx.^2), 2));
        fprintf(' Signal level at TIA and photodiode output before awgn : %.2f dBV2\n', 10*log10(P_erx_moy));
    end
    
fprintf(' \n * time to go through Rx model %.2f seconds \n',toc(rx_start)); tCumul = tCumul + toc(rx_start);

corr_start = tic;

%% Correlation process at the reception
p.rx.ErxLen = size(Erx,2);
Erx=Erx.';
p.tx.detectionindex = [0 0]; %For synchronization and detection of fibre start and fibre end


for k = 1:p.tx.subcarriers
    p = rxSensingCorrelation(p); % Extract estimated Jones matrices
    p.HiTabSB(:,:,k) = p.HiTab;%Store Jones matrices for the k-th subband
    p.rx.nbOvsRayleighReflectorsSB(k) = p.rx.nbOvsReflectors;%Store nb of detected segments
    p.rx.nbDetectedCodes = floor(size(p.HiTab, 2)/p.rx.nbOvsReflectors);%Nb detected codes
    p.rx.nbDetectedCodesSB(k) = p.rx.nbDetectedCodes;
    fprintf('\n * time spent for parameters extraction (%d subbands case): %.2f seconds',p.tx.subcarriers, toc(corr_start)); tCumul = tCumul + toc(corr_start);
    fprintf('\n*** Rx PROC. OUTPUTS (SUBBAND %d/%d) Nb detected segments:%d/%d  Nb detected codes (%d symbols spacing):%d  Overall elapsed time:%.2f seconds ***\n Cumulated time %.2f seconds\n', ...
        k, p.tx.subcarriers, p.rx.nbOvsReflectors, floor(p.fibre.L/p.fibre.spatialRes), p.tx.Ncode, p.rx.nbDetectedCodes, toc(tStart), tCumul);
end

clear global Erx; %clear global variable here since it is not needed any longer from now


%p.HiGen = [Hi(1,1:2:end); Hi(2,1:2:end); Hi(1,2:2:end); Hi(2,2:2:end)]; %Simulated matrices (before TR) %FIXME Should we reduce the number of segments after synchro?