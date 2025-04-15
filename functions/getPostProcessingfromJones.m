function [p,r] = getPostProcessingfromJones(p,r)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that fills the structure p and r with estimation results
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tStart = tic;

%% Calculation of intensity, phase & soft bit tables for each Jones matrix
nbNetUpJonesMat = p.rx.nbDetectedCodes*p.rx.nbOvsReflectors;

%Computing Norm and Determinant of all estimated matrices
intJonesTab = 0.5*sum(p.HiTab(:,1:nbNetUpJonesMat).*conj(p.HiTab(:,1:nbNetUpJonesMat)),1); 
intJonesTab = reshape(intJonesTab,p.rx.nbOvsReflectors,p.rx.nbDetectedCodes);
detJonesTab = p.HiTab(1,1:nbNetUpJonesMat).*p.HiTab(4,1:nbNetUpJonesMat)-p.HiTab(2,1:nbNetUpJonesMat).*p.HiTab(3,1:nbNetUpJonesMat);
detJonesTab = reshape(detJonesTab,p.rx.nbOvsReflectors,p.rx.nbDetectedCodes);
listJonesMatrices = [p.HiTab(1,1:nbNetUpJonesMat); p.HiTab(2,1:nbNetUpJonesMat); p.HiTab(3,1:nbNetUpJonesMat); p.HiTab(4,1:nbNetUpJonesMat)];

%disp(size(listJonesMatrices));
listJonesMatrices = reshape(listJonesMatrices, 4, p.rx.nbOvsReflectors,p.rx.nbDetectedCodes);
%disp(size(listJonesMatrices));
%calculate product of two consecutive jones matrices Hdagger*H for each code
% Initialize the product matrix
productJonesMatrices = zeros(4, p.rx.nbOvsReflectors-1, p.rx.nbDetectedCodes);
p.pola.estimate_pm = zeros(1, p.rx.nbOvsReflectors, p.rx.nbDetectedCodes);

%estimate of phasor
for codeIdx = 1:p.rx.nbDetectedCodes
    % Loop through each segment
    for segIdx = 1:p.rx.nbOvsReflectors
        JonesMatrix = reshape(listJonesMatrices(:, segIdx, codeIdx), 2, 2);
        p.pola.estimate_pm(:, segIdx, codeIdx) = sqrt(det(JonesMatrix));
    end
end 

%estimated phasor for generated matrices
p.pola.estimate_pmGen = zeros(1, p.fibre.nbSegments);
listJonesMatricesGen = reshape(p.HiGen, 2, 2, p.fibre.nbSegments);

% Loop through each segment
for segIdx = 1:p.fibre.nbSegments
    JonesMatrix = listJonesMatricesGen(:, :, segIdx);
    p.pola.estimate_pmGen(:, segIdx) = sqrt(det(JonesMatrix));
end


%plot estimated phasor (after tranmission and reception), phasor of HiGen, and true phasor for code 100

figure;
plot(1:p.rx.nbOvsReflectors, abs(p.pola.estimate_pm(1,:,100)), 'r--');
hold on;
plot(1:p.fibre.nbSegments, abs(p.pola.estimate_pmGen(1,:)), 'g--');
hold on;
plot(1:p.fibre.nbSegments, abs(p.fibre.Ei(1,1:p.fibre.nbSegments).*p.fibre.Ai(1, 1:p.fibre.nbSegments)), 'b+');
title('Estimated and True Phasor for Code 100');
legend('Estimated Phasor', 'Estimated Phasor for HiGen', 'True Phasor');




%plotn
% Loop through each code
for codeIdx = 1:p.rx.nbDetectedCodes
    % Loop through each segment
    for segIdx = 1:p.rx.nbOvsReflectors-1
        % Extract the current and next Jones matrices
        currentJonesMatrix = reshape(listJonesMatrices(:, segIdx, codeIdx), 2, 2);
        nextJonesMatrix = reshape(listJonesMatrices(:, segIdx+1, codeIdx), 2, 2);
        %disp(size(currentJonesMatrix));
        % Calculate the determinants
        detCurrent = currentJonesMatrix(1,1)*currentJonesMatrix(2,2) - currentJonesMatrix(1,2)*currentJonesMatrix(2,1);
        detNext = nextJonesMatrix(1,1)*nextJonesMatrix(2,2) - nextJonesMatrix(1,2)*nextJonesMatrix(2,1);


        
        % Normalize the matrices by the square root of their determinants
        currentJonesMatrix = currentJonesMatrix/sqrt(detCurrent);%(sqrt(abs(detCurrent))*exp(0.5*1j*angle(detCurrent)));
        nextJonesMatrix = nextJonesMatrix/sqrt(detNext);%(sqrt(abs(detNext))*exp(0.5*1j*angle(detNext)));
        %disp('current');
        %disp(currentJonesMatrix);
        %disp('next');
        %disp(nextJonesMatrix);
        %[Uc, Sc, Vc] = svd(currentJonesMatrix);
        %denoised_currentJonesMatrix = Uc*Vc';
        %[Un, Sn, Vn] = svd(nextJonesMatrix);
        %denoised_nextJonesMatrix = Un*Vn';
        %disp(det(currentJonesMatrix));
        % Calculate the product of the dagger of the current matrix and the next matrix
        productMatrix = currentJonesMatrix' * nextJonesMatrix;
        
        % Store the result
        productJonesMatrices(:, segIdx, codeIdx) = productMatrix(:);
    end
end
disp(det(reshape(listJonesMatrices(:,1,1), 2, 2)));
% Initialize the eigenvalues matrix
p.pola.eigenvaluesJonesMatrices = zeros(2, p.rx.nbOvsReflectors-1, p.rx.nbDetectedCodes);
p.pola.birefringenceJonesMatrices = zeros(1, p.rx.nbOvsReflectors-1, p.rx.nbDetectedCodes);
p.pola.birefringenceJonesMatrices2 = zeros(1, p.rx.nbOvsReflectors-1, p.rx.nbDetectedCodes);
% Loop through each code
for codeIdx = 1:p.rx.nbDetectedCodes
    % Loop through each segment
    for segIdx = 1:p.rx.nbOvsReflectors-1
        % Extract the product Jones matrix
        productMatrix = reshape(productJonesMatrices(:, segIdx, codeIdx), 2, 2);
        %productMatrix = productMatrix/sqrt(det(productMatrix));%(sqrt(abs(det(productMatrix)))*exp(0.5*1j*angle(det(productMatrix)));
        % Compute the eigenvalues
        eigenvalues = eig(productMatrix);
        

        % Store the eigenvalues
        p.pola.eigenvaluesJonesMatrices(:, segIdx, codeIdx) = eigenvalues;

        %sin2beta = (1/(2i))*(eigenvalues(1) - eigenvalues(2));
        %p.pola.birefringenceJonesMatrices2(:, segIdx, codeIdx) = asin(sin2beta);
        %disp("eigenvalues");
        %disp(eigenvalues);
        p.pola.birefringenceJonesMatrices(:, segIdx, codeIdx) =  (1/2)*getAnglePlusMinusPiOver2(angle(eigenvalues(1))) - (1/2)*getAnglePlusMinusPiOver2(angle(eigenvalues(2)));%(1/2)*((angle(eigenvalues(1)))) - (angle(eigenvalues(2))); 
    end
end

% Process the generated matrices stored in p.HiGen
% Reshape p.HiGen into 2x2 matrices for each segment
listJonesMatricesGen = reshape(p.HiGen, 2, 2, p.fibre.nbSegments);

% Initialize the product matrix for generated matrices
productJonesMatricesGen = zeros(2, 2, p.fibre.nbSegments - 1);

% Loop through each segment to calculate the product of consecutive matrices
for segIdx = 1:p.fibre.nbSegments - 1
    % Extract the current and next Jones matrices
    currentJonesMatrixGen = listJonesMatricesGen(:, :, segIdx);
    nextJonesMatrixGen = listJonesMatricesGen(:, :, segIdx + 1);

    % Calculate the determinants
    detCurrentGen = det(currentJonesMatrixGen);
    detNextGen = det(nextJonesMatrixGen);

    % Normalize the matrices by the square root of their determinants
    currentJonesMatrixGen = currentJonesMatrixGen / sqrt(detCurrentGen);
    nextJonesMatrixGen = nextJonesMatrixGen / sqrt(detNextGen);

    % Calculate the product of the dagger of the current matrix and the next matrix
    productMatrixGen = currentJonesMatrixGen' * nextJonesMatrixGen;

    % Store the result
    productJonesMatricesGen(:, :, segIdx) = productMatrixGen;
end

% Initialize the eigenvalues and birefringence matrices for generated matrices
p.pola.eigenvaluesJonesMatricesGen = zeros(2, p.fibre.nbSegments - 1);
p.pola.birefringenceJonesMatricesGen = zeros(1, p.fibre.nbSegments - 1);

% Loop through each segment to compute eigenvalues and birefringence
for segIdx = 1:p.fibre.nbSegments - 1
    % Extract the product Jones matrix
    productMatrixGen = productJonesMatricesGen(:, :, segIdx);

    % Compute the eigenvalues
    eigenvaluesGen = eig(productMatrixGen);

    % Store the eigenvalues
    p.pola.eigenvaluesJonesMatricesGen(:, segIdx) = eigenvaluesGen;

    % Compute the birefringence
    p.pola.birefringenceJonesMatricesGen(:, segIdx) = ...
        (1/2) * getAnglePlusMinusPiOver2(angle(eigenvaluesGen(1))) - ...
        (1/2) * getAnglePlusMinusPiOver2(angle(eigenvaluesGen(2)));
end

% Plot the birefringence evolution through the segments for the generated matrices
figure;
plot(1:p.fibre.nbSegments-1, abs(p.pola.birefringenceJonesMatricesGen(1, :)), 'mo');
hold on;
plot(2*(abs(p.fibre.evolBetaTab(2:end))), 'b+');
legend('Birefringence estimated from Generated Matrices', 'True Birefringence of Fiber Segments');
title('Birefringence Evolution through Segments for Generated Matrices');
xlabel('Segment Index');
ylabel('Birefringence (radians)');
grid on;


figure;
for codeIdx = 50:min(50, p.rx.nbDetectedCodes)
    %subplot(1, 1, codeIdx);
    plot(1:p.rx.nbOvsReflectors-1, abs(p.pola.birefringenceJonesMatrices(1,:,codeIdx)), 'r--');
    hold on;
    plot(2*(abs(p.fibre.evolBetaTab(2:end))), 'b');
    title(['Birefringence Evolution through Segments for Code ', num2str(codeIdx)]);
    xlabel('Segment Index');
    ylabel('Birefringence (radians)');
    grid on;
end

% Display birefringence through the codes for the first five segments
figure;
for segIdx = 1:min(5, p.rx.nbOvsReflectors-1)
    subplot(5, 1, segIdx);
    plot(1:p.rx.nbDetectedCodes, squeeze(p.pola.birefringenceJonesMatrices(1,segIdx,:)));
    title(['Birefringence through Codes for Segment ', num2str(segIdx)]);
    xlabel('Code Index');
    ylabel('Birefringence (radians)');
    grid on;
end


figure;
for segIdx = p.pola.segIdx:min(p.pola.segIdx, p.rx.nbOvsReflectors-1)
    figure;
    plot(1:p.rx.nbDetectedCodes-1, squeeze(p.pola.birefringenceJonesMatrices(1,segIdx,1:end-1)));
    title(['Birefringence through Codes for Segment ', num2str(segIdx)]);
    xlabel('Code Index');
    ylabel('Birefringence (radians)');
    grid on;
end
% Display birefringence evolution through the segments for the first code (using birefringenceJonesMatrices2)
%figure;
%for codeIdx = 1:min(1, p.rx.nbDetectedCodes)
%    subplot(1, 1, codeIdx);
%    plot(1:p.rx.nbOvsReflectors-1, abs(p.pola.birefringenceJonesMatrices2(1,:,codeIdx)), 'g--');
%    hold on;
%    plot(2*abs(p.fibre.evolBetaTab(2:end)), 'b');
%    title(['Birefringence Evolution (Method 2) through Segments for Code ', num2str(codeIdx)]);
%    xlabel('Segment Index');
%    ylabel('Birefringence (radians)');
%    grid on;
%end

% Display birefringence through the codes for the first five segments (using birefringenceJonesMatrices2)
% figure;
% for segIdx = 1:min(5, p.rx.nbOvsReflectors-1)
%     subplot(5, 1, segIdx);
%     plot(1:p.rx.nbDetectedCodes, squeeze(p.pola.birefringenceJonesMatrices2(1,segIdx,:)));
%     title(['Birefringence (Method 2) through Codes for Segment ', num2str(segIdx)]);
%     xlabel('Code Index');
%     ylabel('Birefringence (radians)');
%     grid on;
% end
% 
% % Display birefringence evolution through the segments for code number 100 (using birefringenceJonesMatrices2)
% figure;
% codeIdx = 100;
% if codeIdx <= p.rx.nbDetectedCodes
%     plot(1:p.rx.nbOvsReflectors-1, abs(p.pola.birefringenceJonesMatrices2(1,:,codeIdx)), 'g--');
%     hold on;
%     plot(2*abs(p.fibre.evolBetaTab(2:end)), 'b');
%     title(['Birefringence Evolution (Method 2) through Segments for Code ', num2str(codeIdx)]);
%     xlabel('Segment Index');
%     ylabel('Birefringence (radians)');
%     grid on;
% else
%     warning('Code index exceeds the number of detected codes.');
% end

% Plot the standard deviation of the birefringence across time
%figure;
%stdBirefringence = std(p.pola.birefringenceJonesMatrices(1, :, 1:end-1), 0, 3); % Compute std across time (3rd dimension)
%plot(1:p.rx.nbOvsReflectors-1, stdBirefringence, 'b-o');
%title('Standard Deviation of Birefringence Across Time');
%xlabel('Segment Number');
%ylabel('Standard Deviation of Birefringence (radians)');
%grid on;

% % Plot the standard deviation of the birefringence across time (method 2)
% figure;
% stdBirefringence2 = std(p.pola.birefringenceJonesMatrices2(1, :, 1:end-1), 0, 3); % Compute std across time (3rd dimension)
% plot(1:p.rx.nbOvsReflectors-1, stdBirefringence2, 'b-o');
% title('Standard Deviation of Birefringence Across Time (method 2)');
% xlabel('Segment Number');
% ylabel('Standard Deviation of Birefringence (radians)');
% grid on;

r.intJonesTabOvs =          intJonesTab; 
p.rx.AvgIntensPerReflectorTabOvs = mean(intJonesTab,2); %Time averaged intensity per reflector  [1 x Nbsegt]
p.rx.AvgAbsDetPerReflectorTabOvs = mean(abs(detJonesTab),2); %Time averaged abs(det) per reflector  [1 x Nbsegt]


%% Interpolation before downsampling (of singleband of combined bands data)
% A correct downsampling implies an interpolation prior to the decimation step
detTab = detJonesTab; 
intTab = intJonesTab;

if p.rx.ovsFactor > 1 && p.rx.interpolation %FIXME not tested
    
    Hf = getSpectralRaisedCos(p.rx.ovsFactor, p.rx.interpolationRolloff, p.rx.nbOvsReflectors);%A spectral raised cos fct to filter out high spatial frequencies according to ovsFactor and rollOffFactor
    HfRep = repmat(Hf,1,p.rx.nbDetectedCodes);%Repeat Hf a nb of times equal to the nb of codes
    intJonesTabSpatialFiltered = ifft(fft(intTab).*HfRep); %Spatial filtering in the spectral domain
    detJonesTabSpatialFiltered = ifft(fft(detTab).*HfRep); %Spatial filtering in the spectral domain
    
    p.rx.AvgAbsDetPerReflectorTabOvs = mean(abs(detJonesTabSpatialFiltered),2); %Time averaged abs(det) per reflector
    p.rx.AvgIntensPerReflectorTabOvs = abs(mean(intJonesTabSpatialFiltered,2));% 'omitnan'); %Time averaged intensity per reflector
    
    if p.displ.spatialfilter
        %Display for filtering test
        fStep=1/p.rx.nbOvsReflectors;
        fAxis = (0:fStep:1-fStep).*(p.rx.ovsFactor/p.fibre.spatialRes);%Spatial frequency axis -> proportional to the inverse of the (Rx ovs) spatial resolution
        figure(220); hold off;
        %plot(fAxis,20*log10(abs(fft(detJonesTabSpatialFiltered)))); hold on;
        plot(fAxis,20*log10(mean(abs(fft(intJonesTabSpatialFiltered).'))),'b.-'); hold on;
        plot(fAxis,20*log10(mean(abs(fft(detJonesTabSpatialFiltered).'))),'g.-'); hold on;
        plot(fAxis,20*log10(Hf),'r.-'); hold on;
        axis([fAxis(1) fAxis(end) -80 0]); legend('Int', 'Det', 'Filter');
        xlabel('Spatial frequency (Hz)'); ylabel('PSD (dB)');
        title('PSD of abs(det) along the distance dimension'); grid on;
        %End display
    end
    
else %no interpolation
    detJonesTabSpatialFiltered = detTab;
    intJonesTabSpatialFiltered = intTab;
end %ovsFactor>1

clear intTab detTab;

%% Selection of best reflectors for decimation
if p.rx.ovsFactor > 1 && logical(p.rx.decimation) %FIXME not tested
    
    if p.rx.crit_detJones
        IntValOvs = p.rx.AvgAbsDetPerReflectorTabOvs;
    else %p.rx.crit = intensity
        IntValOvs = p.rx.AvgIntensPerReflectorTabOvs;
    end
    
    if p.rx.decimation == 1 %regular downsampling
        downSampMask = true(1, p.rx.nbOvsReflectors); downSampMask(2:p.rx.ovsFactor:end) = false;
    elseif p.rx.decimation % == 2 %downsampling=2: det criteria decimation
        nbBlocks = floor(p.rx.nbOvsReflectors/p.rx.ovsFactor);%Nb selected reflectors (target)
        blockLength = floor(p.rx.nbOvsReflectors/nbBlocks);%Averaged distance between 2 selected reflectors
        [~,XX]=max(reshape(IntValOvs(1: blockLength*nbBlocks), [blockLength,nbBlocks]));%XX is a tab of max values per block of consecutive reflectors
        downSampMask = sum(  ((1:p.rx.nbOvsReflectors) == single([XX + (0:nbBlocks-1)*blockLength p.rx.nbOvsReflectors]).' ), 1);%Index table for selected reflectors (1 selection per block)
        downSampMask = logical(downSampMask);% + circshift((downSampMask == 2),-1)); % patch in case (last) segment is selected twice
        clear XX;
    end%decimation mode
    if p.tx.subcarriers == 1
        detJonesTab = detJonesTabSpatialFiltered(downSampMask,:); %Dim:[SegtsxTime], Decimation to go back to the native spatial resolution
        intJonesTab = intJonesTabSpatialFiltered(downSampMask,:); %Decimation to go back to the native spatial resolution
        r.intJonesTab = intJonesTab; 
        p.rx.nbReflectors = size(detJonesTab,1); %sum(downSampMask,'all'); %
    else %multicarrier %FIXME not tested
        rotaDet_sum = detJonesTabSpatialFiltered; %will be downsampled below
        r.intJonesTabOvs = intJonesTabSpatialFiltered;
        p.rx.nbReflectors = sum(downSampMask,'all'); % size(rotaDet_sum,1); %
    end
    clear HfRep intJonesTabSpatialFiltered detJonesTabSpatialFiltered IntValOvs;
else %ovsFact == 1 or decimation==0
    p.rx.nbReflectors = p.rx.nbOvsReflectors;
    downSampMask = true(1,p.rx.nbReflectors);
end

nbNetJonesMat =  p.rx.nbReflectors*p.rx.nbDetectedCodes;
p.rx.AvgIntensPerReflectorTab = p.rx.AvgIntensPerReflectorTabOvs(downSampMask);
p.rx.AvgAbsDetPerReflectorTab = p.rx.AvgAbsDetPerReflectorTabOvs(downSampMask);
r.fullDetJonesTab =         detJonesTab;

%% Selection criterion for subset of "best" reflectors extraction
if p.rx.crit_detJones
    IntVal = p.rx.AvgAbsDetPerReflectorTab;
else %p.rx.crit = intensity
    IntVal = p.rx.AvgIntensPerReflectorTab;
end
p.rx.AvgSoftBitPerReflectorTab = IntVal./max(IntVal);%Give a soft decision value to each estimated angle. For now, this is just abs(det)) normalized to spread between 0 & 1 (1=full confidence)

%% Selection of best reflectors: averaging or choose the best among p.displ.lowResolFactor
if p.displ.lowResolFactor > 1
    if p.rx.averagingLowres %FIXME not tested
       p.rx.nbOvsSelectedReflectors =  floor(p.rx.nbReflectors/p.displ.lowResolFactor);
       residual = mod(p.rx.nbReflectors, p.displ.lowResolFactor); 
       lowresDet_full = reshape(detJonesTab(1:end-residual,:), p.displ.lowResolFactor , p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes); %dim [nbNetJonesMat, lowresfact]
       lowresDet_time = reshape(detJonesTab(1:end-residual,:), p.displ.lowResolFactor, p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); %dim [lowresfact,segts, time]

       %V1: mean value
       %detJonesTab = squeeze(mean(lowresDet,2));
       %V2: rotated vector sum
       reliablePhase = sum(lowresDet_time.*abs(lowresDet_time),3)./sum(abs(lowresDet_time),3); %sum in time of weighted det.
       rotval_tab = conj(reliablePhase)./abs(reliablePhase); % Angle of which the matrices are to be rotated. size [NbCarriers x NbSegments]
       rotated_detJonesTab = lowresDet_full .* repmat(rotval_tab,1,p.rx.nbDetectedCodes);  %Directly rotate determinant values from detJonesTabSB, so that we don't have to re-compute the determinants from rotated HiTabSum
       detJonesTab = reshape(sum(rotated_detJonesTab,1),p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes)/p.displ.lowResolFactor; %[Nb segts x Nb codes] New constrictively summed det  %Dimensions : [nbSegments x NbJonesMat]
   
       %p.displ.selectedIdxTab = (1:p.rx.nbOvsSelectedReflectors);
       p.displ.selectedIdxTab = single(1:p.rx.nbReflectors); % Ensure all reflectors are considered

       r.intJonesTab = abs(detJonesTab); 
    else %choose best segment
        nbBlocks = floor(p.rx.nbReflectors/p.displ.lowResolFactor);%Nb selected reflectors (target)
        blockLength = floor(p.rx.nbReflectors/nbBlocks);%Averaged distance between 2 selected reflectors
        [~,XX]=max(reshape(IntVal(1: blockLength*nbBlocks), [blockLength,nbBlocks]));%XX is a tab of max values per block of consecutive reflectors
        p.displ.selectedIdxTab = single([XX + (0:nbBlocks-1)*blockLength p.rx.nbReflectors]);%Index table for selected reflectors (1 selection per block)
        p.rx.nbOvsSelectedReflectors = length(p.displ.selectedIdxTab);%Nb of actually selected reflectors
        r.intJonesTab = abs(detJonesTab);
        clear XX;
    end
else
    p.displ.selectedIdxTab = single(1:p.rx.nbReflectors);
    p.rx.nbOvsSelectedReflectors = p.rx.nbReflectors;
end
clear IntVal;

if p.displ.testsimomimo == 0
    phiJonesTabSelect = 0.5*angle(detJonesTab(p.displ.selectedIdxTab,:));%Absolute phase is calculated & stored for the selected Jones matrices only
    %disp(size(phiJonesTabSelect));
else
    phiJonesTabSelect = 0.5*angle(detJonesTab);%Absolute phase is calculated & stored for all Jones matrices
    
end

p.pola.detCurrent = zeros(1, p.rx.nbOvsSelectedReflectors-1, p.rx.nbDetectedCodes);
productJonesMatrices = zeros(4, p.rx.nbOvsSelectedReflectors-1, p.rx.nbDetectedCodes);
% Test of birefringence estimation for best reflectors
JonesMatricesSelect = listJonesMatrices(:, p.displ.selectedIdxTab, :); % Extract the selected Jones matrices
p.pola.JonesMatricesSelect = reshape(JonesMatricesSelect, 4, p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); % Reshape the selected Jones matrices
p.pola.testJonesMatricesSelect = zeros(4, p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes);
disp(size(p.pola.JonesMatricesSelect));
for codeIdx = 1:p.rx.nbDetectedCodes
    % Loop through each segment
    for segIdx = 1:p.rx.nbOvsSelectedReflectors-1
        % Extract the current and next Jones matrices
        currentJonesMatrix = reshape(p.pola.JonesMatricesSelect(:, segIdx, codeIdx), 2, 2);
        nextJonesMatrix = reshape(p.pola.JonesMatricesSelect(:, segIdx+1, codeIdx), 2, 2);
        
        % Calculate the determinants
        p.pola.detCurrent(:, segIdx, codeIdx) = det(currentJonesMatrix);
        detNext = det(nextJonesMatrix);
        
        % Normalize the matrices by the square root of their determinants
        currentJonesMatrix = currentJonesMatrix / sqrt(p.pola.detCurrent(:, segIdx, codeIdx));
        nextJonesMatrix = nextJonesMatrix / sqrt(detNext);
        p.pola.testJonesMatricesSelect(:, segIdx, codeIdx) = reshape(currentJonesMatrix, 4, 1);
        % Calculate the product of the dagger of the current matrix and the next matrix
        productMatrix = currentJonesMatrix' * nextJonesMatrix;
        
        % Store the result
        productJonesMatrices(:, segIdx, codeIdx) = productMatrix(:);
    end
end

% Initialize the eigenvalues matrix
p.pola.eigenvaluesJonesMatricesSelect = zeros(2, p.rx.nbOvsSelectedReflectors-1, p.rx.nbDetectedCodes);
p.pola.birefringenceJonesMatricesSelect = zeros(1, p.rx.nbOvsSelectedReflectors-1, p.rx.nbDetectedCodes);

% Loop through each code
for codeIdx = 1:p.rx.nbDetectedCodes
    % Loop through each segment
    for segIdx = 1:p.rx.nbOvsSelectedReflectors-1
        % Extract the product Jones matrix
        productMatrix = reshape(productJonesMatrices(:, segIdx, codeIdx), 2, 2);
        
        % Compute the eigenvalues
        eigenvalues = eig(productMatrix);
        
        % Store the eigenvalues
        p.pola.eigenvaluesJonesMatricesSelect(:, segIdx, codeIdx) = eigenvalues;
        
        % Compute the birefringence
        p.pola.birefringenceJonesMatricesSelect(:, segIdx, codeIdx) = ...
            (1/2) * getAnglePlusMinusPiOver2(angle(eigenvalues(1))) - ...
            (1/2) * getAnglePlusMinusPiOver2(angle(eigenvalues(2)));

            % Compute the birefringence
        %p.pola.birefringenceJonesMatricesSelect(:, segIdx, codeIdx) = ...
        %    (1/2) * getAnglePlusMinusPiOver2(angle(eigenvalues(1))) - ...
        %    (1/2) * getAnglePlusMinusPiOver2(angle(eigenvalues(2)));
    end
end

% Calculate the briefringence evolution through the segments for the selected reflectors for the generated matrices
p.pola.detCurrentGen = zeros(1, p.rx.nbOvsSelectedReflectors-1);
productJonesMatricesGen = zeros(4, p.rx.nbOvsSelectedReflectors-1);
% Test of birefringence estimation for generated matrices
listJonesMatricesGen = reshape(p.HiGen, 4, p.fibre.nbSegments);
% Extract the selected Jones matrices for generated matrices
JonesMatricesGenSelect = listJonesMatricesGen(:, p.displ.selectedIdxTab); % Extract the selected Jones matrices
p.pola.JonesMatricesGen = reshape(JonesMatricesGenSelect, 4, p.rx.nbOvsSelectedReflectors); % Reshape the selected Jones matrices
p.pola.testJonesMatricesGenSelect = zeros(4, p.rx.nbOvsSelectedReflectors);
for segIdx = 1:p.rx.nbOvsSelectedReflectors-1
    % Extract the current and next Jones matrices
    currentJonesMatrixGen = reshape(p.pola.JonesMatricesGen(:,segIdx), 2, 2);
    nextJonesMatrixGen = reshape(p.pola.JonesMatricesGen(:,segIdx + 1), 2, 2);
    
    % Calculate the determinants
    p.pola.detCurrentGen(:, segIdx) = det(currentJonesMatrixGen);
    detNextGen = det(nextJonesMatrixGen);
    
    % Normalize the matrices by the square root of their determinants
    currentJonesMatrixGen = currentJonesMatrixGen / sqrt(p.pola.detCurrentGen(:, segIdx));
    nextJonesMatrixGen = nextJonesMatrixGen / sqrt(detNextGen);
    p.pola.testJonesMatricesGen(:, segIdx) = reshape(currentJonesMatrixGen, 4, 1);
    
    % Calculate the product of the dagger of the current matrix and the next matrix
    productMatrixGen = currentJonesMatrixGen' * nextJonesMatrixGen;
    
    % Store the result
    productJonesMatricesGen(:, segIdx) = productMatrixGen(:);
end

% Initialize the eigenvalues matrix
p.pola.eigenvaluesJonesMatricesGen = zeros(2, p.rx.nbOvsSelectedReflectors-1);
p.pola.birefringenceJonesMatricesGenSelect = zeros(1, p.rx.nbOvsSelectedReflectors-1);

% Loop through each segment
for segIdx = 1:p.rx.nbOvsSelectedReflectors-1
    % Extract the product Jones matrix
    productMatrixGen = reshape(productJonesMatricesGen(:,segIdx), 2, 2);
    
    % Compute the eigenvalues
    eigenvaluesGen = eig(productMatrixGen);
    
    % Store the eigenvalues
    p.pola.eigenvaluesJonesMatricesGen(:, segIdx) = eigenvaluesGen;
    
    % Compute the birefringence
    p.pola.birefringenceJonesMatricesGenSelect(:, segIdx) = ...
        (1/2) * getAnglePlusMinusPiOver2(angle(eigenvaluesGen(1))) - ...
        (1/2) * getAnglePlusMinusPiOver2(angle(eigenvaluesGen(2)));
end


%adapt p.evol.betaTab to the selected reflectors : need to sum the birefringence of the non selected reflectors in between consecutive selected reflectors
p.fibre.evolBetaTabSelect = zeros(1, p.rx.nbOvsSelectedReflectors);
% Compute the birefringence evolution between selected reflectors
for segIdx = 2:p.rx.nbOvsSelectedReflectors
    % Sum the birefringences of the non-selected segments between the selected reflectors
    startIdx = p.displ.selectedIdxTab(segIdx - 1);
    endIdx = p.displ.selectedIdxTab(segIdx);
    p.fibre.evolBetaTabSelect(segIdx) = sum(p.fibre.evolBetaTab(startIdx+1:endIdx));
end
% Set the birefringence of the first selected reflector
p.fibre.evolBetaTabSelect(1) = p.fibre.evolBetaTab(p.displ.selectedIdxTab(1));

% Plot the birefringence evolution through the segments for the first code
figure;
for codeIdx = 154:min(154, p.rx.nbDetectedCodes)
    plot(1:p.rx.nbOvsSelectedReflectors-1, abs(p.pola.birefringenceJonesMatricesSelect(1, :, codeIdx)), 'r--');
    hold on;
    plot(abs(getAnglePlusMinusPiOver2((((2*p.fibre.evolBetaTabSelect(2:end)))))), 'b+');
    hold on;
    plot(1:p.rx.nbOvsSelectedReflectors-1, abs(p.pola.birefringenceJonesMatricesGenSelect(1, :)), 'g');
    title(['Absolute Birefringence Evolution through Segments for Code ', num2str(codeIdx)]);
    legend('Estimated Matrices', 'sum of true birefringence','Generated Matrices');
    xlabel('Segment Index');
    ylabel('Birefringence (radians)');
    grid on;
end

% figure;
% for codeIdx = p.pola.codeIdx-1:min(p.pola.codeIdx-1, p.rx.nbDetectedCodes)
%     plot(1:p.rx.nbOvsSelectedReflectors-1, abs(p.pola.birefringenceJonesMatricesSelect(1,:,codeIdx)), 'r--');
%     hold on;
%     plot(1:p.rx.nbOvsSelectedReflectors-1, abs(p.pola.birefringenceJonesMatricesGenSelect(1, :)), 'b');
%     title(['Absolute Birefringence Evolution through Segments for Code ', num2str(codeIdx)]);
%     legend('Estimated Matrices', 'Generated Matrices');
%     xlabel('Segment Index');
%     ylabel('Birefringence (radians)');
%     grid on;
% end


%standard deviation on best reflectors
figure;
p.pola.stdBirefringenceSelect = std(p.pola.birefringenceJonesMatricesSelect(1, :, 21:end-1), 0, 3); % Compute std across time (3rd dimension), discarding the first 20 codes
plot(1:p.rx.nbOvsSelectedReflectors-1, p.pola.stdBirefringenceSelect, 'b-o');
title('Standard Deviation of Birefringence Across Time for selected reflectors only');
xlabel('Segment Index');
ylabel('Standard Deviation of Birefringence (radians)');
grid on;


%% Save softBit table(s)
r.maxdet = max(abs(detJonesTab(p.displ.selectedIdxTab,:)),[],'all');
r.selectAbsDetTab_nonorm =  abs(detJonesTab(p.displ.selectedIdxTab,:));  %store all softbits
r.selectAbsDetTab = abs(detJonesTab(p.displ.selectedIdxTab,:))/r.maxdet; %store averaged softbits
r.selectIntFrobTab = r.intJonesTab(p.displ.selectedIdxTab,:); 
r.softBitDiffSelect = [r.selectAbsDetTab(1,:); min(r.selectAbsDetTab(2:end,:),r.selectAbsDetTab(1:end-1,:))];  %r.softBitDiffSelect = r.softBitDiffSelect./max(r.softBitDiffSelect,[],'all');%sum softbit per segment (softbit for diffphase)

r.softBitSelectSTD = std(r.selectAbsDetTab.'); r.muABSDET = mean(r.selectAbsDetTab,2);
r.threshold_comb = r.maxdet/p.tx.subcarriers^2 * p.rx.detTheshold; %doesnt depend on OFDM comparison mode

%% Compute SVD (of raw jones matrices)
r.SVD.exists = 0;
if p.displ.getSVD
    select_ = repmat(p.rx.nbReflectors.*(0:1:p.rx.nbDetectedCodes-1),p.rx.nbOvsSelectedReflectors,1)+p.displ.selectedIdxTab.';
    select_ = sort(reshape(select_,1,[]));
    
    HiSelect = p.HiTab(:,select_);
    
    Ui = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    Vi = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    Li = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    for n=1:p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes %size(selectedIdxLarge,2)
        [Ui(:,:,n),Li(:,:,n),Vi(:,:,n)]=svd([HiSelect(1,n) HiSelect(2,n); HiSelect(3,n) HiSelect(4,n)]);
    end
    r.SVD.exists = 1;
    r.SVD.U = Ui;
    r.SVD.L = Li;
    r.SVD.V = Vi;
    
    Ug = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    Vg = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    Lg = zeros(2,2,p.rx.nbOvsSelectedReflectors);
    for n=1:p.fibre.nbSegments
        [Ug(:,:,n),Lg(:,:,n),Vg(:,:,n)]=svd([p.HiGen(1,n) p.HiGen(2,n); p.HiGen(3,n) p.HiGen(4,n)]);
    end
    r.SVD.Ug = Ug;
    r.SVD.Lg = Lg;
    r.SVD.Vg = Vg;
end %p.displ.getSVD

%% Test differences SIMO MIMO
if p.displ.testsimomimo %FIXME not tested
    temp = p.HiTab;
    p.HiTab = p.HiTab(:,1:nbNetJonesMat);
    if p.fulllength
        if p.rx.siso %XX ou YY - either RXX, RXY or RYY, RYX are null
            if p.rx.Xpol && p.tx.Xpol %SISO XX
                phiJonesTabSelect = angle(p.HiTab(1,:));%p
                intJonesTab =  p.HiTab(1,:).*conj(p.HiTab(1,:));
            elseif p.rx.Xpol || p.tx.Xpol %SISO XY or SISO YX
                phiJonesTabSelect = angle(p.HiTab(2,:) + p.HiTab(3,:));%
                intJonesTab = p.HiTab(2,:).*conj(p.HiTab(2,:))+p.HiTab(3,:).*conj(p.HiTab(3,:));
            else %SISO YY
                phiJonesTabSelect = angle(p.HiTab(1,:) + p.HiTab(4,:));%
                intJonesTab = p.HiTab(4,:).*conj(p.HiTab(4,:));
            end
        elseif p.rx.simo
            phiJonesTabSelect = angle(p.HiTab(1,:) + p.HiTab(2,:) +p.HiTab(3,:) + p.HiTab(4,:));%
        elseif p.rx.miso
            if p.rx.Xpol %MISO X
                phiJonesTabSelect = angle(p.HiTab(1,:) + p.HiTab(3,:) );%
                intJonesTab =  p.HiTab(1,:).*conj(p.HiTab(1,:))+p.HiTab(3,:).*conj(p.HiTab(3,:));
            else %MISO Y
                phiJonesTabSelect = angle(p.HiTab(2,:) + p.HiTab(4,:) );
                intJonesTab = p.HiTab(2,:).*conj(p.HiTab(2,:))+p.HiTab(4,:).*conj(p.HiTab(4,:));
            end
        else %mimo
        end
        phiJonesTabSelect = reshape(phiJonesTabSelect,[p.rx.nbReflectors,p.rx.nbDetectedCodes]);%Reshape determinant vector to get reflectors per row and codes per column
        p.rx.AvgIntensPerReflectorTab = mean(reshape(intJonesTab,p.rx.nbReflectors,p.rx.nbDetectedCodes)  ,2);
        
        phiJonesTabSelect = phiJonesTabSelect(p.displ.selectedIdxTab,:);
        p.HiTab = temp; %repair line 86
    else %not full length
        if p.rx.siso %either RXX, RXY or RYY, RYX are null
            diffPhiTabSelect = reshape(p.HiTab(:,2) +p.HiTab(:,3),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);%Reshape determinant vector to get reflectors per row and codes per column
            diffPhiTabSelect = angle(diffPhiTabSelect);%angle for one pola only
        elseif p.rx.simo
            diffPhiTabSelect = reshape(angle(p.HiTab(:,1) + p.HiTab(:,2) +p.HiTab(:,3) + p.HiTab(:,4)),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);%Reshape determinant vector to get reflectors per row and codes per column
        else %%MIMO
            diffPhiTabSelect = phiJonesTabSelect;%RxCorrXX.*RxCorrYY-RxCorrYX.*RxCorrXY);
        end
        %Then process phase : unwrap, windowing
        diffPhiTabSelect = getAnglePlusMinusPiOver2(([diffPhiTabSelect(1,:) ; diff(diffPhiTabSelect)]).').';%Get differential phase between the selected reflectors
        diffPhiTabSelect = (unwrap(2*diffPhiTabSelect.')/2).';%Unwrap to get rid of PI phase jumps between reflectors
        diffPhiTabSelect = raisedCosFct(diffPhiTabSelect,p.displ.edgeRatio);
        
        p.stdDiffPhiTabSelect = diffPhiTabSelect;%std(diffPhiTabSelect(:,round(0.1*p.rx.nbDetectedCodes):round(0.9*p.rx.nbDetectedCodes)).');%Differential phase std of the selected backscatters after HP filtering (Edges in time ignored)
        return
    end %fulllength
end%p.displ.testsimomimo

%% Phase unwrapping, windowing and filtering

%1) Angle +/- pi/2 and unwrap, make diff phase and then Raised Cosined edge Windowing
diffPhiTabSelect = getAnglePlusMinusPiOver2(([phiJonesTabSelect(1,:) ; diff(phiJonesTabSelect)]).').';%Get differential phase between the selected reflectors
if strcmpi(p.env,'model') 
    diffPhiTabSelect = ((diffPhiTabSelect-mean(diffPhiTabSelect,2)));
end
diffPhiTabSelect = (unwrap(2*diffPhiTabSelect.')/2).';%Unwrap to get rid of PI phase jumps between reflectors
diffPhiTabSelect = raisedCosFct(diffPhiTabSelect,p.displ.edgeRatio); %window with defined edgeratio (cf. global flags)
clear phiJonesTabSelect;

%2) Apply filtering (High-pass & low-pass in a single step) if needed
if (p.rx.f_cutoff >=0 || p.rx.f_cutoff_end >=0)
    halfNbCodes = ceil(0.5*p.rx.nbDetectedCodes);
    nbEdgeVals = 0; nbEdgeVals_end = 0;
    if p.rx.f_cutoff == 0
        nbEdgeVals = 1;%Get rid of DC component
    elseif p.rx.f_cutoff >0
        nbEdgeVals = round(halfNbCodes*p.rx.f_cutoff*2*p.tx.Tcode);
    end
    
    if (p.rx.f_cutoff_end > p.rx.f_cutoff) && (p.rx.f_cutoff_end<0.5/p.tx.Tcode)
        nbEdgeVals_end = round(halfNbCodes*(1-p.rx.f_cutoff_end*2*p.tx.Tcode));
    end
    
    p.rx.spWeightingTab = single([zeros(1,nbEdgeVals) ones(1,halfNbCodes-nbEdgeVals-nbEdgeVals_end) zeros(1,nbEdgeVals_end)]);%First part of the mask (from DC to fMax=fSamp/2)
    if mod(p.rx.nbDetectedCodes,2)%Impose a symetric spectral mask over [DC:fSamp] to get a real filtered signal after ifft
        p.rx.spWeightingTab = ([p.rx.spWeightingTab p.rx.spWeightingTab(end:-1:2)]); %if nbDetectedCodes is odd
    else
        p.rx.spWeightingTab = ([p.rx.spWeightingTab 1 p.rx.spWeightingTab(end:-1:2)]); %if nbDetectedCodes is even
    end
    
    %Filtering : apply window mask in the spectral domain, then ifft
    diffPhiTabSelect = ifft(repmat(p.rx.spWeightingTab.', [1,p.rx.nbOvsSelectedReflectors]).*fft(diffPhiTabSelect.')).';
end



%% Indicator parameters
p.stdDiffPhiTabSelect = std(diffPhiTabSelect(:,round(0.1*p.rx.nbDetectedCodes):round(0.9*p.rx.nbDetectedCodes)).');%Differential phase std of the selected backscatters after HP filtering (Edges in time ignored)
r.diffPhiTabSelect = diffPhiTabSelect; %used for displays (time/freq on selected segments)

%% Jones to Mueller transformation
if p.displ.polar
    downSampHSum = HiTabSum(:,repmat(downSampMask,1,p.rx.nbDetectedCodes));
    downSampHi = p.HiTab(:,repmat(downSampMask,1,p.rx.nbDetectedCodes));
    selectedSegments = sum( ((1:p.rx.nbReflectors) == p.displ.selectedIdxTab.'),1); %Table of indices
    % patch in case (last) segment is selected twice
    selectedSegments = logical(selectedSegments + circshift((selectedSegments == 2),-1));
    if p.tx.subcarriers == 1
        selectedMatrices = downSampHi(:,logical(repmat(selectedSegments, 1, p.rx.nbDetectedCodes)));
    else
        selectedMatrices = downSampHSum(:,logical(repmat(selectedSegments, 1, p.rx.nbDetectedCodes)));
    end
    
    p = getMuellerParam(p,reshape(selectedMatrices ,2,2,[]));
    clear downSampHi downSampHSum
end




%p = getMuellerParam(p,reshape(p.HiTab ,2,2,[]));
% Add this code at the end of your getMuellerParam function

% Plotting the standard deviation of differential polarization parameters
%figure;
%subplot(3,1,1);
% plot(p.pola.stddiffPolaS1);
% title('Standard Deviation of Differential Polarization S1');
% xlabel('Reflector Index');
% ylabel('Standard Deviation');
% 
% subplot(3,1,2);
% plot(p.pola.stddiffPolaS2);
% title('Standard Deviation of Differential Polarization S2');
% xlabel('Reflector Index');
% ylabel('Standard Deviation');
% 
% subplot(3,1,3);
% plot(p.pola.stddiffPolaS3);
% title('Standard Deviation of Differential Polarization S3');
% xlabel('Reflector Index');
% ylabel('Standard Deviation');



%% Display averaged intensity & differential phase vs time of selected backscatters (coarse resolution)

r.max_dist_idx = p.rx.nbReflectors*(p.displ.maxloc/p.fibre.L);
r.max_time_idx = p.rx.nbDetectedCodes;

if p.displ.freqvsTimeDist %FIXME not tested
    normFact_dBPerHz = (p.rx.fSamp/p.rx.ErxLen)*(1/sqrt(size(diffPhiTabSelect,1)));
    freqdiffPhiTabSelect = normFact_dBPerHz.*abs(fft(diffPhiTabSelect.').');
    freqdiffPhiTabSelect = 20*log10(freqdiffPhiTabSelect(:,1:ceil(0.5*p.rx.nbDetectedCodes)));%Store the 1st half spectrum part only
    
    fStep = 1/p.rx.nbDetectedCodes; fTab = single(0:fStep:1-fStep);
    fTab = (1/p.tx.Tcode)*fTab(1:size(freqdiffPhiTabSelect,2));
    
    p.displ.fIdx=p.displ.fIdx+1;figure(p.displ.fIdx);
    colorbar; colormap(jet);
    contour(ftab,p.displ.selectedIdxTab*(0.5*p.fibre.cFiber/p.rx.fSamp),freqdiffPhiTabSelect,[-15,-10,0,5,10,15,20,35],'fill','on');
    
    minFreqLeveldB = -30;%A min level threshold to handle the color map
    offsetColMap = 20;%An offset to transpose the color map
    freqdiffPhiTabSelect(freqdiffPhiTabSelect < minFreqLeveldB) = minFreqLeveldB-1;
    mesh(fTab,p.displ.selectedIdxTab*(0.5*p.fibre.cFiber/p.rx.fSamp),freqdiffPhiTabSelect,real(freqdiffPhiTabSelect-offsetColMap));
    xlabel('Freq (Hz)'); ylabel('distance (m)');zlabel('Power Spectral Density (dB rad^2/Hz)');
    axis([p.displ.minFreq p.displ.maxFreq p.displ.location p.displ.maxloc -30 max(freqdiffPhiTabSelect,[],'all')-1]);
    
    clear fTab freqdiffPhiTabSelect;
end

if p.displ.averagingvsTime  %FIXME not tested
    r.max_dist_idx = p.rx.nbReflectors*(p.displ.maxloc/p.fibre.L);
    r.max_time_idx = p.rx.nbDetectedCodes;
    compLevelDist = 1;% to scale the magnitude display of the reflector phase (but jointly scale the fiber distance observed)
    p.displ.fIdx=p.displ.fIdx+1;
    if p.displ.scaled_map
        prompt = {'linear scaling factor ?', 'Exponential scaling factor ?'};
        dlgtitle = 'Scaled time/dist map inputs';
        definput = {'1', '3'};
        answer = inputdlg(prompt,dlgtitle,[1 40],definput);
        
        p=displayDiffPhi2D_scale(diffPhiTabSelect, str2double(answer{1}), str2double(answer{2}), r.max_dist_idx, r.max_time_idx, p);
    else
        p=displayDiffPhi2D_cat(diffPhiTabSelect, compLevelDist, r.max_dist_idx, r.max_time_idx, p);
    end
end

clear compLevelDist ;%r.max_dist_idx r.max_time_idx;

%% End of function
fprintf('\n** Post-processing completed in %.1f seconds \n', toc(tStart));
fprintf('\n*** POSTPROCESSED (%d BANDS) Nb segments:%d, Nb detected codes (duration %.0f �s):%d, spatial resolution %.2fm (native %.2fm), mech bandwidth %.2fkHz, fibre length %.1fkm \n \n', ...
    p.tx.subcarriers, p.rx.nbOvsSelectedReflectors, p.tx.Ncode/p.tx.fSymb*1e6, p.rx.nbDetectedCodes,p.fibre.spatialRes*p.displ.lowResolFactor,  p.fibre.spatialRes, 1e-3/(2*p.tx.Tcode), p.fibre.L*1e-3);

end

