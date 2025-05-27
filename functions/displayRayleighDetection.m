function p = displayRayleighDetection(p,r)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays features from extracted Jones matrices
%
% Authors:
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by A. Sahu - 2024 adrish.sahu@ip-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Differential phase std versus distance plot
if p.rx.ovsFactor == 1 || (p.rx.decimation == 0)
    p.rx.nbReflectors = p.rx.nbOvsReflectors;
    if p.rx.averagingLowres
        dist_axis = 0.001*0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab*p.displ.lowResolFactor;
    else
        dist_axis = 0.001*0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab;
    end
    dist_axis_all = 0.001*0.5*p.fibre.cFiber/p.rx.fSamp*(1:p.rx.nbOvsSelectedReflectors);
else
    if p.rx.averagingLowres
        dist_axis = 0.001*0.5*p.fibre.cFiber/p.tx.fSymb*p.displ.selectedIdxTab*p.displ.lowResolFactor;
    else
        dist_axis = 0.001*0.5*p.fibre.cFiber/p.tx.fSymb*p.displ.selectedIdxTab;
    end
    dist_axis_all = 0.001*0.5*p.fibre.cFiber/p.tx.fSymb*(1:p.rx.nbOvsSelectedReflectors);
end

p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
plot(dist_axis(2:end-1), p.stdDiffPhiTabSelect(2:end-1),'.-'); hold on;
ylabel('Diff Phase Std (rad)'); xlabel('Distance (km)'); grid on;

%% Jones matrices plot
% SINGLE_JonesMatTab = repmat([p.HiGen, [0;0;0;0]], 1, p.rx.nbDetectedCodes);%p.HiTab(:,1:p.rx.nbDetectedCodes*p.rx.nbOvsReflectors); %
% h11err = p.HiTab(1,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors) - SINGLE_JonesMatTab(1,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors); %absolute error
% h11_relerr = reshape(abs(h11err)./abs(SINGLE_JonesMatTab(1,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors)), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); %relative err
% 
% h12err = p.HiTab(2,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors) - SINGLE_JonesMatTab(2,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors); %absolute error
% h12_relerr = reshape(abs(h12err)./abs(SINGLE_JonesMatTab(2,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors)), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); %relative err
% 
% h21err = p.HiTab(3,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors) - SINGLE_JonesMatTab(3,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors); %absolute error
% h21_relerr = reshape(abs(h21err)./abs(SINGLE_JonesMatTab(3,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors)), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); %relative err
% 
% h22err = p.HiTab(4,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors) - SINGLE_JonesMatTab(4,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors); %absolute error
% h22_relerr = reshape(abs(h22err)./abs(SINGLE_JonesMatTab(4,1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors)), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes); %relative err
% 
% p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
% plot(abs(h11err), '*'); hold on;
% plot(abs(h12err), 'o'); hold on;
% plot(abs(h21err), 'x'); hold on;
% plot(abs(h22err), '<'); hold on;
% legend('h11','h12','h21','h22');
% ylabel('Absolute error on |Jones mat terms|'); xlabel('Time index'); grid on;


% p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
% plot(dist_axis, reshape(abs(h11err),p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes), '*'); hold on;
% plot(dist_axis, reshape(abs(h12err), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes), 'o'); hold on;
% plot(dist_axis, reshape(abs(h21err), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes), 'x'); hold on;
% plot(dist_axis, reshape(abs(h22err), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes), '<'); hold on;
% legend('h11','h12','h21','h22');
% ylabel('Absolute error on |Jones mat terms|'); xlabel('Distance (km)'); grid on;
% 
% p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
% plot(dist_axis_all(1:end-1), h11_relerr(1:end-1,2:end), '*'); hold on;
% plot(dist_axis_all(1:end-1), h12_relerr(1:end-1,2:end), 'o'); hold on;
% plot(dist_axis_all(1:end-1), h21_relerr(1:end-1,2:end), 'x'); hold on;
% plot(dist_axis_all(1:end-1), h22_relerr(1:end-1,2:end), '<'); hold on;
% ylabel('Relative error on Jones mat terms'); yyaxis right;
% plot(dist_axis_all(1:end-1),p.stdDiffPhiTabSelect(1:end-1));
% legend('h11','h12','h21','h22');
% ylabel('\sigma \phi (rad)'); xlabel('Distance (km)'); grid on;

if p.rx.ovsFactor == 1 || (p.rx.decimation == 0)
    p.rx.nbReflectors = p.rx.nbOvsReflectors;
    dist_axis = 0.001*0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab;
    dist_axis_all = 0.001*0.5*p.fibre.cFiber/p.rx.fSamp*(1:p.rx.nbReflectors);
else %ovsFact = 2
    dist_axis = 0.001*0.5*p.fibre.cFiber/p.tx.fSymb*p.displ.selectedIdxTab;
    dist_axis_all = 0.001*0.5*p.fibre.cFiber/p.tx.fSymb*(1:p.rx.nbReflectors);
end
%% Display standard deviation as fct of fiber length
taxis = (0:p.rx.nbDetectedCodes-1).*p.tx.Tcode;

if p.displ.StDv && p.env.EXP
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
    plot(dist_axis(2:end-1), p.stdDiffPhiTabSelect(2:end-1),'.-'); hold on;
    ylabel('Diff Phase Std (rad)');
    if p.displ.perSubband
        hold off;
        plot(dist_axis(2:end-1), r.stdDiffPhiTabSelectSB(:,2:end-1),'-');hold on;
        plot(dist_axis(2:end-1), p.stdDiffPhiTabSelect(2:end-1),'r.-','LineWidth',2); hold on;
        legend_str = string(1:p.tx.subcarriers+1); legend_str(end) = "Combination";
        legend('Subband '+legend_str);
    else  

    end
    xlabel('Fiber distance (km)');
    title(sprintf('Optical phase standard deviation as fct of fiber length \n(BP filter:%d:%dHz) MinPhi:%.2e MaxPhi:%.2e AvgPhi:%.2e', p.rx.f_cutoff,p.rx.f_cutoff_end, min(p.stdDiffPhiTabSelect(2:end-1)), max(p.stdDiffPhiTabSelect(2:end-1)), mean(p.stdDiffPhiTabSelect(2:end-1))));
    grid on;
end

if p.displ.StDv && p.env.MODEL
    legend_str = string(1:p.tx.subcarriers+1); legend_str(end) = 'Combined';
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
    %eventuellement : for n=1:nbSimu
    if p.displ.perSubband
        plot(dist_axis(2:end-2),r.stdDiffPhiTabSelectSB(:,2:end-2),'.-');hold on;
    end
    plot(dist_axis(2:end-1), p.stdDiffPhiTabSelect(2:end-1),'r.-','LineWidth',2); hold on;
    xlabel('Fiber distance (km)'); ylabel('Diff Phase Std (rad)');
    legend('Subband '+legend_str);
    if p.displ.perSubband
        axis([0 max(dist_axis(1,2:end)) 0 (1e-7)+(1+0.01)*max(r.stdDiffPhiTabSelectSB(:,2:end-2),[],'all')]); grid on;
    else
        axis([0 max(dist_axis(1,2:end)) 0 (1e-7)+(1+0.01)*max(p.stdDiffPhiTabSelect(:,2:end-2),[],'all')]); grid on;
    end
    
    title(sprintf('Optical phase standard deviation as fct of fiber length (HP filter:%dHz).', p.rx.f_cutoff));
end

if p.displ.StDvtoStrain
    epsilon_strain = p.tx.Lambda * p.stdDiffPhiTabSelect / (4*pi*p.ng*0.78) ; %*p.fibre.spatialRes
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
    plot(dist_axis(2:end-1), 1e9*epsilon_strain(2:end-1),'.-'); hold on;
    xlabel('Fiber distance (km)'); ylabel('Strain detected (n\epsilon)');
    grid on; 
end %stdvtostrain

%% Display intensity per reflector
if p.displ.IntensityPerReflector
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
    %plot(dist_axis(2:end-1),p.rx.AvgAbsDetPerReflectorTab(p.displ.selectedIdxTab(2:end-1)),'b.-');hold on;
    plot(dist_axis(2:end-1),p.stdDiffPhiTabSelect(2:end-1),'g.-');hold on;
    xlabel('Distance (km)'); ylabel('Avg (in time) intensity per reflector');
    axis([0 dist_axis(end) 0 1]); grid on;
    %axis([1 0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end-1)) 0 1]); grid on;
    title(sprintf('Phi & Intensity as fct of selected reflector index (LP filter:%dHz).\n MinPhi:%.2e MaxPhi:%.2e AvgPhi:%.2e',...
        p.rx.f_cutoff, min(p.stdDiffPhiTabSelect(2:end-1)), max(p.stdDiffPhiTabSelect(2:end-1)), mean(p.stdDiffPhiTabSelect(2:end-1))));
end

%% Display average RBS intensity as fct of fiber distance
if p.displ.intensityDistance %FIXME not tested
    if p.displ.intensityDistdB
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
        plot(dist_axis_all,10*log10(p.rx.AvgIntensPerReflectorTab./max(p.rx.AvgIntensPerReflectorTab)), '-');
        if p.displ.perSubband
            legent_str_int = string(1:p.tx.subcarriers+2); legent_str_int(1) = "Avg Subbands";
            hold on;%ff;
            for k = 1:p.tx.subcarriers
                plot(dist_axis_all,10*log10( mean(  squeeze(r.fullIntJonesTab(k,:,:))./max(p.rx.AvgIntensPerReflectorTab), 2)), '--'); hold on;
                legent_str_int(k+1) = "Int sb "+string(k);
            end
            legent_str_int(p.tx.subcarriers+2) = "Rotated sum vector";
            plot(dist_axis_all,10*log10(mean(r.intJonesTab,2)./max(p.rx.AvgIntensPerReflectorTab)), '.-');
            legend(legent_str_int);
        end
        xlabel('distance(km)'); ylabel('RBS intensity native res (dB)');
        title('RBS intensity as function of distance, dB');
        grid on;
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
        plot(dist_axis,10*log10(mean(r.intJonesTab(p.displ.selectedIdxTab,:),2)./max(mean(r.intJonesTab(p.displ.selectedIdxTab,:),2))), '-');
        xlabel('distance(km)'); ylabel('RBS intensity (dB) ');
        title('RBS intensity as function of distance, dB');
        grid on;
        
    else %linear intensity displayed
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold on;

        axis([ 0 dist_axis_all(end) 0 1*p.tx.subcarriers]);
        xlabel('distance(km)'); ylabel('RBS intensity'); grid on;
        title('RBS intensity as function of distance, average intensity on all codes');
    end
end

%% Display Phi as fct of intensity and |det|
if p.displ.PhiIntensity
    if p.tx.subcarriers == 1
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
        plot(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab(2:end-1)),p.stdDiffPhiTabSelect(2:end-1),'b.','MarkerSize', 12);hold on;
        xlabel('Rayleigh Intensity'); ylabel('Diff Phase Std (rad)');grid on;
        title(sprintf('Phi as fct of intensity for selected reflectors.\n MinPhi:%.2e MaxPhi:%.2e AvgPhi:%.2e', min(p.stdDiffPhiTabSelect(2:end-1)), max(p.stdDiffPhiTabSelect(2:end-1)), mean(p.stdDiffPhiTabSelect(2:end-1))));
        sprintf('Phi & Intensity as fct of selected reflector index (LP filter:%dHz). MinPhi:%.2e MaxPhi:%.2e AvgPhi:%.2e', p.rx.f_cutoff, min(p.stdDiffPhiTabSelect(2:end-1)), max(p.stdDiffPhiTabSelect(2:end-1)), mean(p.stdDiffPhiTabSelect(2:end-1)));
        
        p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
        semilogx(p.rx.AvgAbsDetPerReflectorTab(p.displ.selectedIdxTab(2:end-1)),p.stdDiffPhiTabSelect(2:end-1),'o','MarkerSize', 5);hold on;
        %        semilogx(p.rx.AvgAbsDetPerReflectorTab(p.displ.selectedIdxTab(2:end-1))/r.maxdet,p.stdDiffPhiTabSelect(2:end-1),'o','MarkerSize', 5);hold on;
        xlabel('|det|'); ylabel('Diff Phase Std (rad)');grid on;
    else
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
        loglog(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab(2:end-1)),mean(r.stdDiffPhiTabSelectSB(:,2:end-1),1),'o','MarkerSize', 4);hold on;
        loglog(mean(r.intJonesTab(p.displ.selectedIdxTab(2:end-1),:),2),p.stdDiffPhiTabSelect(2:end-1),'.','MarkerSize', 12);hold on;
        legend('Mean all subcarriers', 'Combined subcarriers');
        xlabel('Rayleigh Intensity'); ylabel('Diff Phase Std (rad)');grid on;
        title('\sigma_\phi as function of intensity');
        
        p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
        loglog(p.rx.AvgAbsDetPerReflectorTab(p.displ.selectedIdxTab(2:end-1)),p.stdDiffPhiTabSelect(2:end-1),'o','MarkerSize', 5);hold on;
        %        loglog(p.rx.AvgAbsDetPerReflectorTab(p.displ.selectedIdxTab(2:end-1))/r.maxdet,p.stdDiffPhiTabSelect(2:end-1),'o','MarkerSize', 5);hold on;
        xlabel('|det|'); ylabel('Diff Phase Std (rad)');grid on;
    end
end

%% Store selected phases (audio) and softbit
if p.displ.audio || p.displ.psd_audio
    selectedReflectorIdx1=min(p.displ.tabAudio);%206;%656;%197;%982;
    if length(p.displ.tabAudio) == 1
        selectedReflectorIdx2 = selectedReflectorIdx1+1 ;
    else
        selectedReflectorIdx2=max(p.displ.tabAudio);%313;%191;%10;%206;%656;%197;%982;
    end
    audioBuffer1 = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx1,:),p.displ.edgeRatio);
    audioBuffer2 = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx2,:),p.displ.edgeRatio);
    audioBuffer1m = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx1-1,:),p.displ.edgeRatio);
    audioBuffer2m = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx2-1,:),p.displ.edgeRatio);
    audioBuffer1p = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx1+1,:),p.displ.edgeRatio);
    audioBuffer2p = raisedCosFct(r.diffPhiTabSelect(selectedReflectorIdx2+1,:),p.displ.edgeRatio);
    
    spatialRes = 0.5*p.fibre.cFiber*single(2*logical(p.rx.decimation))/(2*p.tx.fSymb); %spatialRes = spatialRes*p.displ.lowResolFactor;
    legend_str_soft = string(1:4);
    legend_str_soft(2) =  sprintf(' %0.f m',p.displ.selectedIdxTab(selectedReflectorIdx1)*spatialRes );
    legend_str_soft(3) =  sprintf(' %0.f m',p.displ.selectedIdxTab(selectedReflectorIdx2)*spatialRes );
    legend_str_soft(1) =  sprintf(' %0.f m',p.displ.selectedIdxTab(selectedReflectorIdx1-1)*spatialRes );
    legend_str_soft(4) =  sprintf(' %0.f m',p.displ.selectedIdxTab(selectedReflectorIdx2+1)*spatialRes );
    
    
    intBuffer_1 = r.intJonesTab(selectedReflectorIdx1,:);
    intBuffer_2 = r.intJonesTab(selectedReflectorIdx2,:);
    intBuffer_1m = r.intJonesTab(selectedReflectorIdx1-1,:);
    intBuffer_2p = r.intJonesTab(selectedReflectorIdx2+1,:);
end

if p.displ.softvalues
    absDetBuffer_soft1 = r.selectAbsDetTab(selectedReflectorIdx1,:);
    absDetBuffer_soft2 = r.selectAbsDetTab(selectedReflectorIdx2,:);
    absDetBuffer_soft1m = r.selectAbsDetTab(selectedReflectorIdx1-1,:);
    absDetBuffer_soft2p = r.selectAbsDetTab(selectedReflectorIdx2+1,:);
    if p.tx.subcarriers > 1
        absdet_all_soft1m = squeeze(r.fullAbsDetTable(:,selectedReflectorIdx1-1,:));
        absdet_all_soft1 = squeeze(r.fullAbsDetTable(:,selectedReflectorIdx1,:));
        absdet_all_soft2 = squeeze(r.fullAbsDetTable(:,selectedReflectorIdx2,:));
        absdet_all_soft2p = squeeze(r.fullAbsDetTable(:,selectedReflectorIdx2+1,:));
    end
    softbitBuffer_soft1 = r.softBitDiffSelect(selectedReflectorIdx1,:); %is product of successive det. Use r.selectAbsDetTab for point det.
    softbitBuffer_soft2 = r.softBitDiffSelect(selectedReflectorIdx2,:);
    softbitBuffer_soft1m = r.softBitDiffSelect(selectedReflectorIdx1-1,:);
    softbitBuffer_soft2p = r.softBitDiffSelect(selectedReflectorIdx2+1,:);
    
    
    legend_str_soft = string(1:4);
    legend_str_soft(2) =  sprintf(' %.1f m',p.displ.selectedIdxTab(selectedReflectorIdx1)*spatialRes );
    legend_str_soft(3) =  sprintf(' %.1f m',p.displ.selectedIdxTab(selectedReflectorIdx2)*spatialRes );
    legend_str_soft(1) =  sprintf(' %.1f m',p.displ.selectedIdxTab(selectedReflectorIdx1-1)*spatialRes );
    %legend_str_soft(4) =  sprintf('Soft %0.2f m',p.displ.selectedIdxTab(selectedReflectorIdx2-1)*spatialRes );
    %legend_str_soft(5) =  sprintf('Soft %0.2f m',p.displ.selectedIdxTab(selectedReflectorIdx1+1)*spatialRes );
    legend_str_soft(4) =  sprintf(' %.1f m',p.displ.selectedIdxTab(selectedReflectorIdx2+1)*spatialRes );
end%softvalues

%% Display selected phases (audio)
if p.displ.audio
    %%% Display selected phase of best subcarrier per segment according to softbit criteria
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
    plot(taxis, audioBuffer1m ,'--'); hold on;
    plot(taxis,audioBuffer1,'m-'); hold on;
    plot(taxis,audioBuffer2,'b-'); hold on;
    %plot(taxis, audioBuffer_soft2m ,'--'); hold on;
    %plot(taxis, audioBuffer_soft1p ,'--'); hold on;
    plot(taxis, audioBuffer2p ,'--'); %hold on;
    legend("Buffer"+legend_str_soft);
    xlabel('t (s)'); ylabel('audioBuffer phase');grid on;
    title('Time domain selected audio buffers');
end %audio

%% Display psd at points where perturbation is applied
if p.displ.psd_audio
    df = 1/(p.tx.Tcode*r.max_time_idx); faxis = -1/(2*p.tx.Tcode):df:1/(2*p.tx.Tcode)-df;
    buffersize = size(audioBuffer1,2);
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
    plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer1m)))),'--'); hold on;
    plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer1)))),'m-'); hold on;
    plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer2)))),'b-'); hold on; grid on;
    %plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer_soft2m)))),'--'); hold on;
    %plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer_soft1p)))),'--'); hold on;
    plot(faxis,fftshift(20*log10((1/sqrt(buffersize))*abs(fft(audioBuffer2p)))),'--'); hold on; grid on;
    legend("Phase"+legend_str_soft);
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB rad^2/Hz)');
    xlim([0 0.5/p.tx.Tcode]); ylim([-40 10]);grid on;%ylim([-20 60]);grid on;
    title('PSD of soft-selected audio buffers');
    
end

%% Display reliability metric on selected phases
if p.displ.reliability %can be done for multi- and single carrier . Single : r.softBitDiffSelect Multi : r.softBitDiffSelect
    %softBitTable=colormap(spring); softBitTable=softBitTable(end:-2:1,:); %hot %winter %parula
    softBitTable=colormap(jet); softBitTable=softBitTable(33:end,:);
    resol = size(softBitTable,1) ; %resolution : number of color levels
    stepX = (1000-100)/(resol-1);
    tabThr=single(1./([1 (100+stepX : stepX : 1000)])); %table of threshold reliability values
    
    sizeMap=reshape(repmat((5:1:12),resol/8,1),1,resol); % point size increases as reliability decreases
    
    tab_ = [selectedReflectorIdx1-1, selectedReflectorIdx1, selectedReflectorIdx2, selectedReflectorIdx2+1];
    for nn = 1:4
        selectedIndex = tab_(nn); %selectedReflectorIdx2; %reflector index inspected
        selectedsoft = r.softBitDiffSelect(selectedIndex,:); %softbitBuffer_soft2; %corresponding soft bit table
        selectedphase = raisedCosFct(r.diffPhiTabSelect(selectedIndex,:),p.displ.edgeRatio); %corresponding phase
        tabMask0 = repmat((1:length(selectedsoft)).',1, length(tabThr)).*(selectedsoft.'<tabThr); %Mask : location of unreliable values
        
        p.displ.fIdx=p.displ.fIdx+1;figure(p.displ.fIdx);hold off;
        plot(taxis,selectedphase, 'g:','lineWidth', 1); hold on; %Phase evolution, in time
        for pp=1: length(tabThr) %for all reliability levels, assess reliability for all instants
            %UnderThresh = nonzeros(tabMask0(:,pp)); %Time indices where phase is under tabThr value
            Phasepoints = selectedphase(nonzeros(tabMask0(:,pp))); %Corresponding phase values at these points
            plot((nonzeros(tabMask0(:,pp))-1)*p.tx.Tcode,Phasepoints,'.','color', softBitTable(pp,:),'MarkerSize',sizeMap(pp)); hold on;
        end
        colormap(softBitTable); colorbar('Ticks',[0,1], 'TickLabels',{'+','-'});
        %colorbar.Label.String = 'Phase estimation';
        ylabel('Phase (rad)'); xlabel(sprintf('Time (s)'));% index (%.2es res., %.2fs overall)',p.tx.Tcode,p.tx.Tcode*r.max_time_idx));
        xlim([0, p.tx.Tcode*r.max_time_idx]);
        title(sprintf('Phase vs time at d=%.1fm (index %d) - (%d:%dHz) BP filter - soft mean %.2f e-3', p.displ.selectedIdxTab(selectedIndex)*(0.5*p.fibre.cFiber/p.rx.fSamp),selectedIndex, p.rx.f_cutoff, p.rx.f_cutoff_end,1e3*mean(selectedsoft)));
        grid on;
    end
    clear selectedIndex selectedsoft selectedphase UnderThresh Phasepoints pp;
end %p.displ.reliability


%% Overall reliability indicator
if p.displ.reliability_dist
    Msb =  min(r.softBitDiffSelect,[],2); %Minimal soft bit value per segment
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
    plot(dist_axis(2:end), Msb(2:end),'.-'); hold on;
    if p.tx.subcarriers > 1
        fact_norm = max( abs(r.fullDetJonesTab) , [], 'all') / r.maxdet; %to apply same normalization to combined and individual subbands
        absdettable = r.fullAbsDetTable(1:end,p.displ.selectedIdxTab,:) * fact_norm; %selected segments for individual subcarriers %rawabsdet for no normalization
        %absdettable = absdettable/max(absdettable, [],'all'); %normalized
        absdettable = cat(2, absdettable(:,1,:), min(absdettable(:,1:end-1,:),absdettable(:,2:end,:))); %differential soft bit
        Mdt = min(absdettable,[], 3);
        plot(dist_axis(2:end), Mdt(:,2:end), '--'); hold on;
        legend_str_ovrel = string(0:p.tx.subcarriers);
        legend_str_ovrel(1) = 'Combined bands';
        legend('Min rel. for band ' + legend_str_ovrel);
    else
        legend('Min reliability');
    end
    xlabel('Fibre distance (km)'); ylabel('MIN Soft bit value per segment');
    grid on;
    
end

if p.displ.relDistr
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
    h = histogram(reshape(r.softBitDiffSelect,1,[]), 'Normalization','probability');
    legend('Combined bands'); grid on;
    ylabel('Probability'); xlabel('|det| values');
    
    if p.tx.subcarriers > 1
        fact_norm = max( abs(r.fullDetJonesTab) , [], 'all') / r.maxdet; %to apply same normalization to combined and individual subbands
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
        absdettable = r.fullAbsDetTable(1:end,p.displ.selectedIdxTab,:)*fact_norm; %selected segments for individual subcarriers %rawabsdet for no normalization
        absdettable = cat(2, absdettable(:,1,:), min(absdettable(:,1:end-1,:),absdettable(:,2:end,:))); %differential soft bit
        hist( cat(1, reshape(absdettable,p.tx.subcarriers,[]), reshape(r.softBitDiffSelect,1,[])   ).',h.NumBins); hold on;
        strr = string(1:p.tx.subcarriers+1); strr(end) = 'All';
        legend('sub-band ' +  strr); grid on;
        ylabel('Occurences'); xlabel('|det| values');
    end
    
end %Rel distrib

if p.displ.SNRdistrib
    %% %%%% snr phase
    if p.env.EXP
        p.stockStd = p.stdDiffPhiTabSelect;
        p.stockabsDet = r.selectAbsDetTab_nonorm(2:end-1);
        p.stockabsDet_t = r.selectAbsDetTab(2:end-1)*r.maxdet/r.threshold_comb;
    end
    SNRphase = 10*log10(1./p.stockStd.^2);
    SNRdet = 10*log10(p.stockabsDet);

    p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
    histogram( SNRdet,floor(200),'Normalization','pdf','BinMethod','auto'); hold on;
xlabel(sprintf('SNR det/N^2 (dB)')); ylabel('Probability density function');  grid on; 
%xlabel(sprintf('SNR det (dB), threshold %.2f',p.rx.detTheshold)); ylabel('Probability density function');  grid on; 

    p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
    histogram( p.stockabsDet,floor(200),'Normalization','pdf','BinMethod','auto'); hold on;
    xlabel(sprintf('Determinant, normalized, threshold %.2f',p.rx.detTheshold)); ylabel('Probability density function');

    p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
    histogram( SNRphase(2:end-1),floor(200),'Normalization','pdf','BinMethod','auto'); hold on;
    xlabel('SNR phase'); ylabel('Probability density function');
    
    p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
    histogram( p.stockStd,'Normalization','pdf','BinMethod','auto'); hold on;
    xlabel('\sigma_\phi (rad)'); ylabel('Probability density function');
    
     return;
    %% histogram to pdf
    single_x = h_single.BinEdges(2:end); single_y = h_single.Values;%;
    ofdm2_x = h_ofdm2.BinEdges(2:end); ofdm2_y = h_ofdm2.Values;%(3:end);
    ofdm4_x = h_ofdm4.BinEdges(2:end); ofdm4_y = h_ofdm4.Values;%(3:end);
    %ofdm8_x = h_ofdm8.BinEdges(2:end); ofdm8_y = h_ofdm8.Values;%(3:end);
    figure(); plot(single_x, single_y); legend('SINGLE'); grid on; 
    xlabel('SNR_\epsilon (strain, dB)'); ylabel('Probability density function');
    figure(); plot(ofdm2_x, ofdm2_y); legend('OFDM2');
    figure(); plot(ofdm4_x, ofdm4_y); legend('OFDM4');
    %figure(); plot(ofdm8_x, ofdm8_y); legend('OFDM8');
    xlabel('SNR |det|'); ylabel('Probability density function');
end %SNRdistrib

%% Reliable subcarriers and softbit values
if p.displ.softvalues
    %%% Most reliable subcarriers for selected Idx %%%
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;absDetBuffer_soft1
    plot(taxis, absDetBuffer_soft1m ,'k-'); hold on;
    plot(taxis, absDetBuffer_soft1 ,'m-'); hold on;
    plot(taxis, absDetBuffer_soft2 ,'b-'); hold on;
    plot(taxis, absDetBuffer_soft2p ,'g-'); hold on;
    ylabel('|det| value');
    legend("|det| at "+legend_str_soft);
    if p.tx.subcarriers > 1
        plot(taxis, absdet_all_soft1m ,'k--'); hold on;
        plot(taxis, absdet_all_soft1 ,'m--'); hold on;
        plot(taxis, absdet_all_soft2 ,'b--'); hold on;
        plot(taxis, absdet_all_soft2p ,'g--'); hold on;
        
        legend(cat(2,"|det| at "+legend_str_soft, "aha"));
    else  %one subcarrier
        legend("|det| at "+legend_str_soft);
        title('Soft bit value, in time ');
        xlabel('t (s)'); grid on;
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        
        %yyaxis right;
        plot(taxis, absDetBuffer_soft1m ,'k:'); hold on;
        plot(taxis, intBuffer_1 ,'m:'); hold on;
        plot(taxis, intBuffer_2 ,'b:'); hold on;
        plot(taxis, absDetBuffer_soft2p ,'g:'); hold on;
    end
    
    legend("|det| at "+legend_str_soft);
    title('Soft bit value, in time (+int or +per subband)');
    xlabel('t (s)'); grid on;
 
end
%% Export sound from selected index to working folder
if p.displ.sound
    audioBuffer2 = audioBuffer2./max(abs(audioBuffer2));%Normalize all channels according to max abs found among all channels
    audioBuffer1 = audioBuffer1./max(abs(audioBuffer1));%Normalize all channels according to max abs found among all channels
    FsAudio = round(1/p.tx.Tcode);%Audio sampling rate, integer value expected
    
    audiowrite(p.file.audio2,1*audioBuffer2,FsAudio); %sound(0.1*audioBuffer,FsAudio)
    audiowrite(p.file.audio1,1*audioBuffer1,FsAudio); %sound(0.1*audioBuffer,FsAudio)
end

%% Error calculation and display
% ratio |det| and expected |det|
if strcmpi(p.env,'model') && p.displ.errDet
    theorDet = ((p.fibre.Ai.* p.fibre.Ei).^2)./max(abs((p.fibre.Ai.* p.fibre.Ei).^2),[],'all');
    if p.tx.subcarriers > 1
        theorDet_bands = theorDet(1:p.tx.subcarriers:end);
        for k = 2:p.tx.subcarriers %for subband1 response in case of OFDM
            theorDet_bands = theorDet_bands + cat(2,theorDet(k:p.tx.subcarriers:end),zeros(size(theorDet_bands,2) - size(theorDet(k:p.tx.subcarriers:end),2)) );
        end
        unitarityLoss = abs(r.fullDetJonesTab)./abs(repmat(theorDet_bands./p.tx.subcarriers, 1, 1, p.rx.nbDetectedCodes));
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
        plot(dist_axis_all(1:size(theorDet,2)), mean(unitarityLoss,3)); hold on;
        xlabel('Distance(km)'); ylabel('|det|/|A_ip_i|^2'); grid on;
    else%one subcarrier
        unitarityLoss = abs(r.fullDetJonesTab(1:size(theorDet,2),:))./abs(repmat(theorDet, p.rx.nbDetectedCodes,1).');
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);%hold off;
        %plot(dist_axis_all(1:size(theorDet,2)), max(unitarityLoss,[],2)); hold on;
        plot(dist_axis_all(1:size(theorDet,2)), mean(unitarityLoss,2)); hold on;
        xlabel('Distance(km)'); ylabel('|det|/|A_ip_i|^2'); grid on;
    end
    yyaxis right;  ylabel('|A_ip_i|^2');
    semilogy(dist_axis_all(1:size(theorDet,2)),abs(theorDet));
end

% Err calc
if strcmpi(p.env,'model') && p.displ.errCalc

    
    if p.displ.lowResolFactor ~= 1
        longIdxTab = [p.displ.selectedIdxTab, zeros(1,(p.rx.nbDetectedCodes-1)*size(p.displ.selectedIdxTab,2))];
        for n = 1:p.rx.nbDetectedCodes-1
            longIdxTab(1,1+n*size(p.displ.selectedIdxTab,2):(n+1)*size(p.displ.selectedIdxTab,2)) =p.displ.selectedIdxTab +n*p.rx.nbReflectors;
        end
    else
        longIdxTab = ones(1,p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors);
    end
     p.HiTab = p.HiTab(:, 1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors); %remove extra segment if needed
    OUTIntVal = 0.5 *( p.HiTab(1,:).*conj(p.HiTab(1,:))+p.HiTab(4,:).*conj(p.HiTab(4,:))+p.HiTab(2,:).*conj(p.HiTab(2,:))+p.HiTab(3,:).*conj(p.HiTab(3,:)) );
    OUTIntVal = reshape(OUTIntVal(1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors),p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes);
    OUTdetVal = p.HiTab(1,:).*p.HiTab(4,:) - p.HiTab(3,:).*p.HiTab(2,:) ;
    OUTdetVal = abs(reshape(OUTdetVal(1:p.rx.nbDetectedCodes*p.rx.nbOvsSelectedReflectors),p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes));
   
    INIntVal =  0.5*(p.HiGen(1,:).*conj(p.HiGen(1,:))+p.HiGen(4,:).*conj(p.HiGen(4,:))+p.HiGen(2,:).*conj(p.HiGen(2,:))+p.HiGen(3,:).*conj(p.HiGen(3,:)) );
 
    INdetVal = abs( p.HiGen(1,:).*p.HiGen(4,:) - p.HiGen(3,:).*p.HiGen(2,:) );
    INbackscatter = abs(p.fibre.Ei .* p.fibre.Ai );     TheorVal = INbackscatter.^2;     TheorVal_time = repmat(INbackscatter.^2,size(OUTdetVal,2),1).';
    
    quantiles_aipi2 = quantile(INbackscatter.^2,[0.25,0.5,0.75,0.9, 1]);
    low_backscatter_index = (INbackscatter.^2 > quantiles_aipi2(1));
    
    norm_factor = max(TheorVal)./max(mean(OUTdetVal,2));
    OUTIntVal_NORM = OUTIntVal(1:size(INbackscatter,2),:)*norm_factor;
    OUTdetVal_NORM = OUTdetVal(1:size(INbackscatter,2),:)*norm_factor;
    
    %display()
    if p.displ.errCalcDispl
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold on;
        plot(TheorVal); plot(INdetVal);  plot(INIntVal); plot(mean(OUTdetVal_NORM,2)); plot(mean(OUTIntVal_NORM,2));
        legend('Aipi theor', 'IN det', 'IN int', 'OUT det', 'OUT int');
        axis([0 p.rx.nbOvsSelectedReflectors 0 max(TheorVal,[],'all')]);
    end
    %calcul des erreurs (relatives)
    errINdet = ( abs(INdetVal) - TheorVal )./TheorVal;
    errINnorm = ( abs(INIntVal) - TheorVal )./TheorVal;
    
    errTRXdet = abs( repmat(abs(INdetVal),size(OUTdetVal,2),1).' - abs(OUTdetVal(1:size(INbackscatter,2),:)) )./repmat(abs(INdetVal),size(OUTdetVal,2),1).';
    errTRXnorm =  abs( repmat(abs(INIntVal),size(OUTdetVal,2),1).' - abs(OUTIntVal(1:size(INbackscatter,2),:)) )./repmat(abs(INIntVal),size(OUTdetVal,2),1).';
    
    errTRXdet_NORM = abs( (repmat(abs(INIntVal),size(OUTdetVal,2),1).') - OUTdetVal_NORM )./TheorVal_time;
    errTRXnorm_NORM =  abs( (repmat(abs(INIntVal),size(OUTdetVal,2),1).') - OUTIntVal_NORM )./TheorVal_time;
    errOUTdet = -( TheorVal_time - OUTdetVal_NORM )./TheorVal_time;
    errOUTnorm = -( TheorVal_time - OUTIntVal_NORM )./TheorVal_time;
    errOUTdet_mean = mean(errOUTdet,2);
    errOUTnorm_mean = mean(errOUTnorm,2);
    
    errOUTdet_lowbs = errOUTdet(low_backscatter_index,:); %only for chose low backscatter index
    errOUTnorm_lowbs = errOUTnorm(low_backscatter_index,:);
    
    AerrINdet = ( abs(INdetVal) - TheorVal );
    AerrINnorm = ( abs(INIntVal) - TheorVal );
    AerrTRXdet_NORM = abs( (repmat(abs(INIntVal),size(OUTdetVal,2),1).') - OUTdetVal_NORM );
    AerrTRXnorm_NORM =  abs( (repmat(abs(INIntVal),size(OUTdetVal,2),1).') - OUTIntVal_NORM );
    AerrOUTdet = abs( TheorVal_time - OUTdetVal_NORM );
    AerrOUTnorm = abs( TheorVal_time - OUTIntVal_NORM );
    AerrOUTdet_mean = mean(AerrOUTdet,2);
    AerrOUTnorm_mean = mean(AerrOUTnorm,2);
    
    AerrOUTdet_lowbs = AerrOUTdet(low_backscatter_index,:); %only for chose low backscatter index
    AerrOUTnorm_lowbs = AerrOUTnorm(low_backscatter_index,:);

    
    %affichage des erreurs
    if p.displ.errCalcDispl
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        plot(dist_axis_all(1:size(INbackscatter,2)), errINnorm,'o-'); hold on;
        plot(dist_axis_all(1:size(INbackscatter,2)), errINdet,'x-'); hold on;
        xlabel('Distance (km)'); ylabel('relative error'); legend('Error on NORM','Error on DET');
        title('Error IN');
       %
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        plot(dist_axis_all(1:size(INbackscatter,2)),mean(errTRXnorm_NORM,2),'o-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),mean(errTRXdet_NORM,2),'x-');
        xlabel('Distance (km)'); ylabel('relative error'); legend('Error on NORM','Error on DET'); title('Error Matrix v Rx');
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        %   plot(dist_axis_all(1:size(INbackscatter,2)),errOUTdet); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),errOUTnorm);
        plot(dist_axis_all(1:size(INbackscatter,2)),errOUTnorm_mean,'k-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),errOUTdet_mean,'p--');
        xlabel('Distance (km)'); ylabel('relative error');
        %     yyaxis right; ylabel('StDv error in time');
        %     semilogy(dist_axis_all(1:size(INbackscatter,2)),-1./TheorVal);
        yyaxis right; ylabel('StDv error in time');
        plot(dist_axis_all(1:size(INbackscatter,2)),std(errOUTnorm.'),'c-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),std(errOUTdet.'),'g--');
        legend('Error on NORM','Error on DET','STD error NORM','STD error DET'); title('Error OUT');
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        %   plot(dist_axis_all(1:size(INbackscatter,2)),errOUTdet); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),errOUTnorm);
        plot(dist_axis_all(1:size(INbackscatter,2)),errOUTnorm_mean,'o-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),errOUTdet_mean,'x-');
        xlabel('Distance (km)'); ylabel('relative error');  legend('Error on NORM','Error on DET'); title('Error Tx v Rx');
    end
    
    mseDetIN = sqrt(mean(AerrINdet.^2));
    mseDetOUT = sqrt(mean(AerrOUTdet.^2, 'all'));
    mseNormIN = sqrt(mean(AerrINnorm.^2));
    mseNormOUT = sqrt(mean(AerrOUTnorm.^2, 'all'));
    
    mseDetOUTbs = sqrt(mean(AerrOUTdet_lowbs.^2, 'all'));
    mseNormOUTbs = sqrt(mean(AerrOUTnorm_lowbs.^2, 'all'));
    
    if p.displ.errCalcDispl
        M = mean(TheorVal,'all');
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        plot(dist_axis_all(1:size(INbackscatter,2)),mean(AerrTRXnorm_NORM,2),'o-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),mean(AerrTRXdet_NORM,2),'x-');
        xlabel('Distance (km)'); ylabel('Absolute error'); legend('Error on NORM','Error on DET'); title(sprintf('Error Matrix v Rx, mean %.2e',M));
        
        p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);
        %   plot(dist_axis_all(1:size(INbackscatter,2)),errOUTdet); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),errOUTnorm);
        plot(dist_axis_all(1:size(INbackscatter,2)),AerrOUTnorm_mean,'o-'); hold on; plot(dist_axis_all(1:size(INbackscatter,2)),AerrOUTdet_mean,'x-');
        xlabel('Distance (km)'); ylabel('Absolute error');  legend('Error on NORM','Error on DET'); title(sprintf('Error Tx v Rx, mean %.2e',M));
    end
    
    fprintf('Absolut MSE values : Det IN %.2e Norm IN %.2e -  Det OUT %.2e Norm OUT %.2e \n', mseDetIN, mseNormIN, mseDetOUT, mseNormOUT);
    p.mseDetOut = mseDetOUT; p.mseNormOUT = mseNormOUT;
    p.mseDetOUTbs = mseDetOUTbs; p.mseNormOUTbs = mseNormOUTbs;
   % return; 
    
    mseDetIN = sqrt(mean(errINdet.^2));
    mseDetOUT = sqrt(mean(errOUTdet.^2, 'all'));
    mseNormIN = sqrt(mean(errINnorm.^2));
    mseNormOUT = sqrt(mean(errOUTnorm.^2, 'all'));
    
    mseDetOUTbs = sqrt(mean(errOUTdet_lowbs.^2, 'all'));
    mseNormOUTbs = sqrt(mean(errOUTnorm_lowbs.^2, 'all'));
    
    fprintf('Relative MSE values : Det OUT %.2e Norm OUT %.2e -- Low BS Det OUT %.2e Norm OUT %.2e \n', mseDetOUT, mseNormOUT, mseDetOUTbs, mseNormOUTbs);
    
    p.mseDetOut = mseDetOUT; p.mseNormOUT = mseNormOUT;
    p.mseDetOUTbs = mseDetOUTbs; p.mseNormOUTbs = mseNormOUTbs;
    
    
    return;
    OUTvsINGain = max(p.rx.AvgIntensPerReflectorTab)/max(INIntVal); %
    %OUTvsINGain = (norm(p.HiTab(1,:))+norm(p.HiTab(2,:))+norm(p.HiTab(3,:))+norm(p.HiTab(4,:)))/(norm(p.HiGen(1,:))+norm(p.HiGen(2,:))+norm(p.HiGen(3,:))+norm(p.HiGen(4,:)) );
    OUTxx = (1/OUTvsINGain) * reshape(p.HiTab(1,:), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes);
    OUTxy = (1/OUTvsINGain) * reshape(p.HiTab(2,:), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes);
    OUTyx = (1/OUTvsINGain) * reshape(p.HiTab(3,:), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes);
    OUTyy = (1/OUTvsINGain) * reshape(p.HiTab(4,:), p.rx.nbOvsSelectedReflectors, p.rx.nbDetectedCodes);
    
    INdet = INxx.*INyy-INxy.*INyx;
    OUTDet = OUTxx.*OUTyy - OUTyx.*OUTxy;
    
    r.errAbsDet = abs(OUTDet*OUTvsINGain-repmat(INdet,size(OUTDet,2),1).');
    r.errAbsRelDet = abs(r.errAbsDet./repmat(INdet,size(OUTDet,2),1).');
    r.meanErrAbsRelDet = mean(r.errAbsRelDet(:,round(0.1*p.rx.nbDetectedCodes):round(0.9*p.rx.nbDetectedCodes)).');%Differential phase error of the selected backscatters after HP filtering (Edges in time ignored)
    
    INPhi = 0.5*angle(INdet);
    OUTPhi = 0.5*angle(OUTDet);
    
    r.errAbsPhi = mod(abs(OUTPhi-repmat(INPhi,size(OUTPhi,2),1).'), 2*pi);
    r.meanErrAbsPhi = mean(r.errAbsPhi(:,round(0.1*p.rx.nbDetectedCodes):round(0.9*p.rx.nbDetectedCodes)).');%Differential phase error of the selected backscatters after HP filtering (Edges in time ignored)
    r.realDet = real(repmat(INdet,size(OUTDet,2),1).'); r.imagDet = imag(repmat(INdet,size(OUTDet,2),1).');
    INdet = repmat(INxx.*INyy-INxy.*INyx,size(OUTPhi,2),1).';
    r.realErrDetTab = abs(real(OUTDet-INdet)./(real(INdet)));%CD 11/19 error relative du determinant
    
    if p.displ.ploterrors
        figure(100);
        plot(INIntVal*OUTvsINGain);hold on;  plot(p.rx.AvgIntensPerReflectorTab);
        title('Verify INint OUTint alignment');
        figure(101);
        semilogy(r.errRelJones); hold on; grid on;
        title('Relative error Intensity'); xlabel('Distance (segt index)');
        figure(1102);
        semilogy(r.errAbsPhi); hold on; grid on;
        %semilogy(mean(r.errAbsPhi(:,2:end-1),2)); hold on; grid on;
        %semilogy(mean(r.realErrDetTab,2)); hold on;
        yyaxis right;
        semilogy(p.stdDiffPhiTabSelect);
        title('Absolute error Phase (mean) copared to std'); legend('errPhi','real err', 'std');xlabel('Distance (segt index)');
        figure(103);
        semilogy(mean(r.realErrDetTab(:,2:end-1),2)); hold on; grid on;
        semilogy(mean(r.imagErrDetTab(:,2:end-1),2));
        title('Error on real and im part of det (mean)'); legend('Real', 'Im');xlabel('Distance (segt index)');
        figure(104);
        semilogy(mean(r.absDiagErrDetTab(:,2:end-1),2)); hold on; grid on;
        semilogy(mean(r.absAntiDiagErrDetTab(:,2:end-1),2));
        title('Error relative to diag and antidiag (mean)'); legend('Diag', 'Antidiag');xlabel('Distance (segt index)');
        figure(105); hold off
        plot(abs(real(INxx.*INyy - INyx.*INxy))); hold on; grid on;
        
        figure(1); hold off;
        yyaxis left;
        semilogy(r.meanErrAbsPhi); hold on; grid on; ylabel('Phi absolute error (rad)');
        yyaxis right;
        semilogy(r.meanErrAbsRelDet./abs(INdet)); hold on; grid on; ylabel('Det relative error');
        title(sprintf('Err.rate: Phi(left,mu:%.2e) & Det(right,mu:%.2e), meth.%d',mean(r.meanErrAbsPhi),mean(r.meanErrAbsRelDet./abs(INdet)),p.tx.ProbingMethod)); legend('errPhi','errDet');xlabel('Distance (segment index)');
        figure(2); hold off;
        yyaxis left;
        semilogy(r.meanErrAbsPhi); hold on; grid on; ylabel('Phi absolute error (rad)');
        
        yyaxis right;
        semilogy(p.stdDiffPhiTabSelect); hold on; grid on; ylabel('Phi stdev over time (rad)');
        title(sprintf('Err. rate: Phi(left,mu:%.2e) & std(right,mu:%.2e), meth.%d',mean(r.meanErrAbsPhi),mean(p.stdDiffPhiTabSelect),p.tx.ProbingMethod)); legend('errPhi','stdDiffPhi');xlabel('Distance (segment index)');
        
        figure(3); hold off;
        semilogy(r.meanErrAbsPhi); hold on; grid on; ylabel('Phi absolute error (rad)');
        figure(4); hold off;
        semilogy(r.meanErrAbsRelDet./abs(INdet)); hold on; grid on; ylabel('Det relative error');
        %         figure(3); hold off;
        %         plot(angle(INdet)); hold on; grid on;
    end
end %p.env.MODEL && errCalc


%% Display polarization features
%% Optical polarization standard deviation for S0, S1 S2 S3
if p.displ.polar && p.displ.polarparam
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),p.pola.stddiffPolaS1(2:end),'r.-');hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),p.pola.stddiffPolaS2(2:end),'g.-');hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),p.pola.stddiffPolaS3(2:end),'b.-');hold on;
    xlabel('Fiber distance (km)'); ylabel('Diff Pola Std (rad)');
    %axis([0 0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end-1)) 0 1]); grid on;
    title(sprintf('Optical polarization standard deviation as fct of fiber length'));
    legend(sprintf('StDv of S1 Stokes parameter'),'StDv of S2 Stokes parameter','StDv of S3 Stokes parameter');
    
    p.displ.fIdx = p.displ.fIdx+1; figure(p.displ.fIdx);
    plot(0.001*0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end-1),std(p.pola.S0_t(2:end-1,:).'),'r.-');hold on; %/max(std(p.pola.S0_t.'))
    title(sprintf('Optical polarization intensity (S0) standard deviation as fct of fiber length'));
    legend('StDv S0 parameter (Intensity)');%, 'Average intensity per reflector');
end

%% Display Poincaré sphere
if p.displ.sphere && p.displ.polar
    time_line = p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes;% size(JonesMatTab3D,3); %time_line/sizeTab = nbDetectedCodes
    clear s c;
    window = (1:ceil(2*p.rx.nbDetectedCodes/16));%(ceil(3*p.rx.nbDetectedCodes/8):ceil(1*p.rx.nbDetectedCodes/2));
    win_size = size(window, 2);
    
    S = linspace(50,45,win_size);
    C1 = linspace(1,3,win_size);
    C2 = linspace(2,4,win_size);
    p.pola.S1_t = reshape(p.pola.S1_t,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);
    p.pola.S2_t = reshape(p.pola.S2_t,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);
    p.pola.S3_t = reshape(p.pola.S3_t,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);
    intem = ceil(0.5*(selectedReflectorIdx1+selectedReflectorIdx2));
    p.displ.fIdx = p.displ.fIdx+1;     figure(p.displ.fIdx);
    scatter3(p.pola.S1_t(selectedReflectorIdx1,window),p.pola.S2_t(selectedReflectorIdx1,window),p.pola.S3_t(selectedReflectorIdx1,window),S,C1);hold on;
    legend("Segt "+legend_str_soft(2));
    
    p.displ.fIdx = p.displ.fIdx+1;     figure(p.displ.fIdx);
    scatter3(p.pola.S1_t(intem,window),p.pola.S2_t(intem,window),p.pola.S3_t(intem,window),S,C1);
    legend("Segt "+string( p.fibre.spatialRes* p.displ.selectedIdxTab( intem )  ));
    
    p.displ.fIdx = p.displ.fIdx+1;     figure(p.displ.fIdx);
    scatter3(p.pola.S1_t(selectedReflectorIdx2,window),p.pola.S2_t(selectedReflectorIdx2,window),p.pola.S3_t(selectedReflectorIdx2,window),S,C2);hold on;
    legend("Segt "+legend_str_soft(3));
    
    title('Selected segments');
    %legend(sprintf(num2str((selectedReflectorIdx1))), sprintf(num2str(intem)),sprintf(num2str((selectedReflectorIdx2))));
    clear C S time_line;
end


%% Display Polar PSD
if p.displ.polarPSD && p.displ.polar
    diffS0 = reshape(p.pola.S0_t, 1,[]); diffS0 = [diffS0(1) diff(diffS0)];
    diffS0 = reshape(diffS0, [], p.rx.nbDetectedCodes);
    POLA0audioBuffer1 = raisedCosFct(diffS0(selectedReflectorIdx1,:),p.displ.edgeRatio);
    POLA0audioBuffer2 = raisedCosFct(diffS0(selectedReflectorIdx2,:),p.displ.edgeRatio);grid on;
    
    df = 1/(p.tx.Tcode*r.max_time_idx); faxis = -1/(2*p.tx.Tcode):df:1/(2*p.tx.Tcode)-df;
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %plot(faxis,fftshift(20*log10(abs(fft(audioBuffer1)))),'b-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA0audioBuffer1)))),'m-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA0audioBuffer2)))),'b-'); hold on; grid on;
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB)');
    legend(sprintf(num2str((selectedReflectorIdx1))), sprintf(num2str((selectedReflectorIdx2))));
    xlim([0 0.5/p.tx.Tcode]); ylim([-30 50]);grid on;
    title('PSD selected buffers on S0 parameter');
    
    POLA1audioBuffer1 = raisedCosFct(p.pola.diffPolaS1(selectedReflectorIdx1,:),p.displ.edgeRatio);
    POLA1audioBuffer2 = raisedCosFct(p.pola.diffPolaS1(selectedReflectorIdx2,:),p.displ.edgeRatio);grid on;
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %plot(faxis,fftshift(20*log10(abs(fft(audioBuffer1)))),'b-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA1audioBuffer1)))),'m-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA1audioBuffer2)))),'b-'); hold on; grid on;
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB)');
    legend(sprintf(num2str((selectedReflectorIdx1))), sprintf(num2str((selectedReflectorIdx2))));
    xlim([0 0.5/p.tx.Tcode]); ylim([-30 50]);grid on;
    title('PSD selected buffers on S1 parameter');
    
    POLA2audioBuffer1 = raisedCosFct(p.pola.diffPolaS2(selectedReflectorIdx1,:),p.displ.edgeRatio);
    POLA2audioBuffer2 = raisedCosFct(p.pola.diffPolaS2(selectedReflectorIdx2,:),p.displ.edgeRatio);grid on;
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %plot(faxis,fftshift(20*log10(abs(fft(audioBuffer1)))),'b-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA2audioBuffer1)))),'m-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA2audioBuffer2)))),'b-'); hold on; grid on;
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB)');
    legend(sprintf(num2str((selectedReflectorIdx1))), sprintf(num2str((selectedReflectorIdx2))));
    xlim([0 0.5/p.tx.Tcode]); ylim([-30 50]);grid on;
    title('PSD selected buffers on S2 parameter');
    
    POLA3audioBuffer1 = raisedCosFct(p.pola.diffPolaS3(selectedReflectorIdx1,:),p.displ.edgeRatio);
    POLA3audioBuffer2 = raisedCosFct(p.pola.diffPolaS3(selectedReflectorIdx2,:),p.displ.edgeRatio);grid on;
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %plot(faxis,fftshift(20*log10(abs(fft(audioBuffer1)))),'b-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA3audioBuffer1)))),'m-'); hold on;
    plot(faxis,fftshift(20*log10(abs(fft(POLA3audioBuffer2)))),'b-'); hold on; grid on;
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB)');
    legend(sprintf(num2str((selectedReflectorIdx1))), sprintf(num2str((selectedReflectorIdx2))));
    xlim([0 0.5/p.tx.Tcode]); ylim([-20 60]);grid on;
    title('PSD selected buffers on S3 parameter');
end %polarPSD

%% Display linear retardance, DOP, (Y-X) polar Std
if p.displ.polar && p.displ.polarparam
    stdLinRetardance =  std([p.pola.delta(:,1), diff(p.pola.delta,1,2)].');
    stdCircRetardance =  std([p.pola.rho(:,1), diff(p.pola.rho,1,2)].');
    stdDOP = std(p.pola.DOP.');
    %diffDOP =  std([p.pola.DOP(:,1), diff(p.pola.DOP,1,2)].');
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),stdLinRetardance(2:end),'r.-');hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),stdCircRetardance(2:end),'g.-');hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),stdDOP(2:end),'b.-');hold on;
    xlabel('Fiber distance (km)'); ylabel('Diff Pola Std (rad)');
    %axis([0 0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end-1)) 0 1]); grid on;
    title(sprintf('Optical polarization standard deviation as fct of fiber length'));
    legend('StDv of Linear retardance \delta','StDv of circular retardance \rho','StDv of DOP');
    
    p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx);hold off;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),std(p.pola.deltaPhi(2:end,:).'));hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),std(p.pola.ellipticity(2:end,:).'));hold on;
    plot(0.001*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end),std(p.pola.psy(2:end,:).'));hold on;
    xlabel('Fiber distance (km)'); ylabel('Diff Pola Std (rad)');
    %axis([0 0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp*p.displ.selectedIdxTab(2:end-1)) 0 1]); grid on;
    title(sprintf('Optical polarization standard deviation as fct of fiber length'));
    legend('StDv of phase difference bw pola  x and y \Delta\Phi','StDv of ellipticity e','StDv of \Psi orientation of fast axis');
end

%% Display alarms (?)
if p.displ.alarms
    kTab=find(p.stdDiffPhiTabSelect>p.displ.stdPhiAlarmThres);%Get indices of alarm
    if ~isempty(kTab) &&  kTab(1)==1
        kTab = kTab(2:end); %Get rid of ref phase if selected
    end
    
    nbAlarms = length(kTab);
    if nbAlarms>0
        nbReflectorsPerAlarmMax = 20;
        p.nbReflectorsPerAlarm = NaN(nbAlarms, 1);
        p.displ.selectedIdxTabPerAlarm = NaN(nbAlarms, nbReflectorsPerAlarmMax);
        p.stdDiffPhiTabSelectPerAlarm = NaN(nbAlarms, nbReflectorsPerAlarmMax);
        
        for m=1:nbAlarms %Loop over all the reflectors with a phase std above threshold
            idxStartSearch = p.displ.selectedIdxTab(kTab(m)-1); idxStopSearch = p.displ.selectedIdxTab(kTab(m));
            [diffPhiTabSelectAlarm, p]  = getDiffPhiTabSelectAlarmRegular_cat(detJonesMatTab2D, idxStartSearch, idxStopSearch, m, p);
            if p.nbReflectorsPerAlarm(m)>1
                %% Display averaged intensity & differential phase vs time of selected backscatters (high resolution) for each detected alarm
                p.displ.fIdx=p.displ.fIdx+1;
                displayDiffPhiAlarm2D_cat(diffPhiTabSelectAlarm, idxStartSearch, idxStopSearch, m, compLevelDist, max_time_idx, p)
            end
            if m ==1
                P1 = r.diffPhiTabSelect(kTab(m),:);%diffPhiTabSelectAlarm(1,:);
                P1 = P1-mean(P1); %P1 = P1/sqrt(mean(P1.^2));
                %T = 0:p.tx.Tcode:p.tx.Tcode*(max_time_idx-1);
                figure(100); plot(P1); grid on;
                df = 1/(p.tx.Tcode*max_time_idx); faxis = -1/(2*p.tx.Tcode):df:1/(2*p.tx.Tcode)-df;
                figure(101); plot(faxis,fftshift(20*log10(abs(fft(P1))))); grid on;
                xlim([0 0.5/p.tx.Tcode]); ylim([0 70])
            else
                P1 = r.diffPhiTabSelect(kTab(m),:);%diffPhiTabSelectAlarm(1,:);
                P1 = P1-mean(P1); %P1 = P1/sqrt(mean(P1.^2));
                figure(100); hold on; plot(P1);
                %df = 1/(p.tx.Tcode*max_time_idx); faxis = -1/(2*p.tx.Tcode):df:1/(2*p.tx.Tcode)-df;
                figure(101); hold on; plot(faxis,fftshift(20*log10(abs(fft(P1)))));
            end
        end
        figure(100); xlabel(sprintf('Time index (%.2es res., %.2fs overall)',p.tx.Tcode,p.tx.Tcode*max_time_idx)); ylabel('Phase (rad)'); legend(num2str((p.displ.selectedIdxTab(kTab))'));
        figure(101); xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB)'); legend(num2str((p.displ.selectedIdxTab(kTab))'));
    end
end %alarms
end