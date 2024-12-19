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

% if strcmpi(p.env,'model')
%     if p.fibre.ExcitedSegmentFlag %Add artificial attenuation to emulate fading at perturbation location
%         locationIdxTime = cat(1, (0:p.rx.nbDetectedCodes-1)*p.rx.nbOvsReflectors+2*p.fibre.ExcitedSegmentIdx-1,(0:p.rx.nbDetectedCodes-1)*p.rx.nbOvsReflectors+2*p.fibre.ExcitedSegmentIdx) ;
%         locationTabInTime = reshape( cat(1,locationIdxTime-2,locationIdxTime,locationIdxTime+2), 1, [])  ;
%         p.HiTab(:,locationTabInTime) = p.fibre.artificialFading*p.HiTab(:,locationTabInTime);
%     end
% end

%Computing Norm and Determinant of all estimated matrices
intJonesTab = 0.5*sum(p.HiTab(:,1:nbNetUpJonesMat).*conj(p.HiTab(:,1:nbNetUpJonesMat)),1); 
intJonesTab = reshape(intJonesTab,p.rx.nbOvsReflectors,p.rx.nbDetectedCodes);
detJonesTab = p.HiTab(1,1:nbNetUpJonesMat).*p.HiTab(4,1:nbNetUpJonesMat)-p.HiTab(2,1:nbNetUpJonesMat).*p.HiTab(3,1:nbNetUpJonesMat); 
detJonesTab = reshape(detJonesTab,p.rx.nbOvsReflectors,p.rx.nbDetectedCodes);

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
   
       p.displ.selectedIdxTab = (1:p.rx.nbOvsSelectedReflectors);
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
else
    phiJonesTabSelect = 0.5*angle(detJonesTab);%Absolute phase is calculated & stored for all Jones matrices
end

%% Save softBit table(s)
r.maxdet = max(abs(detJonesTab(p.displ.selectedIdxTab,:)),[],'all');
r.selectAbsDetTab_nonorm =  abs(detJonesTab(p.displ.selectedIdxTab,:));  %store all softbits
r.selectAbsDetTab = abs(detJonesTab(p.displ.selectedIdxTab,:))/r.maxdet; %store averaged softbits
r.selectIntFrobTab = r.intJonesTab(p.displ.selectedIdxTab,:); % backscattered intensity from selected segments 
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
    %contour(ftab,p.displ.selectedIdxTab*(0.5*p.fibre.cFiber/p.rx.fSamp),freqdiffPhiTabSelect,[-15,-10,0,5,10,15,20,35],'fill','on');
    
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
    compLevelDist = 1;%(p.rx.ovsFactor*logical(p.rx.decimation));% to scale the magnitude display of the reflector phase (but jointly scale the fiber distance observed)
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

