function [StartX,StopX,EstimatedFiberLength,p] = getRayleighSpread(absSumRxCorr, ovsCodeLength, nbExpectedCodes, p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect position of Rayleigh start and termination from intensity vector at correlation output
% The values returned defines the start and stop of the first reliable backscattered code with StartX<StopX
%
% Authors:
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(p.env,'exp')
    
    resDet = sum((reshape(absSumRxCorr(1+ovsCodeLength:(nbExpectedCodes-1)*ovsCodeLength), [ovsCodeLength,(nbExpectedCodes-2)])).');%The average of absSumRxCorr over all the transmitted codes (except 1st & last)
    [PeakY,PeakX] = max(resDet); %Find the max intensity (assumed to belong to the fiber Rayleigh backscatter)
    resDetShift = circshift(resDet, ovsCodeLength/2- PeakX);%A copy of resDet but shifted to get the peak in the middle of the buffer (X=ovsCodeLength/2)
    
    %Check if peak is fiber start, end, or a random location along the fiber
    intensWinLen = 10;
    sumIntensPrior = sum(resDetShift(ovsCodeLength/2-intensWinLen-1:ovsCodeLength/2-1));
    sumIntensAfter = sum(resDetShift(ovsCodeLength/2+1:ovsCodeLength/2+intensWinLen+1));
    ratioIntens = sumIntensPrior/sumIntensAfter;%>>1: peak is fiber end; <<1: peak is fiber start; close to 1: peak along the fiber, but location undetermined (probably  a high reflexion due to a bad connection)...
    
    %Detect minimum prior and after the peak (minima are assumed to NOT belong to the fiber Rayleigh backscatter)  => Backscatter is framed
    [~,xminPrior] = min(resDetShift(ovsCodeLength/4+1:ovsCodeLength/2-1)); xminPrior = xminPrior + ovsCodeLength/4;
    [~,xminAfter] = min(resDetShift(ovsCodeLength/2+1:round(3/4*ovsCodeLength))); xminAfter = xminAfter + ovsCodeLength/2;
    
    %Get a low-resolution (averaged) version of resDetShift
    lowResFact = 8;
    tmpTabLowRes = resDetShift(xminPrior:xminAfter);
    tmpTabLowRes = mean(reshape(tmpTabLowRes(1:lowResFact*floor(length(tmpTabLowRes)/lowResFact)),[lowResFact,floor(length(tmpTabLowRes)/lowResFact)]));
    
    %Search max & min ratio between 2 consecutive low resolution samples to localize Rayleigh end & start
    ratioTmpTabLowRes = ([1 tmpTabLowRes(2:end)./tmpTabLowRes(1:end-1)]);
    [StartYLowRes,StartXLowRes] = max(ratioTmpTabLowRes);%Fiber start
    [StopYLowRes,StopXLowRes] = min(ratioTmpTabLowRes);%Fiber end
    
    [~,tmp1X] = max(resDetShift((StartXLowRes-1)*lowResFact + xminPrior - lowResFact : (StartXLowRes-1)*lowResFact + xminPrior + lowResFact));%Look for max intensity position in resDetShift around low-res intitial position
    StartX = (StartXLowRes-1)*lowResFact + xminPrior - lowResFact -1 + tmp1X - (ovsCodeLength/2-PeakX);
    [~,tmp2X] = max(resDetShift((StopXLowRes-1)*lowResFact + xminPrior - 2*lowResFact : (StopXLowRes-1)*lowResFact + xminPrior + lowResFact));%Look for max intensity position in resDetShift around low-res intitial position
    StopX  = (StopXLowRes-1)*lowResFact + xminPrior - 2*lowResFact - 1 + tmp2X - (ovsCodeLength/2-PeakX); %-2*lowResFact instead of lowResFact because of intensity spread after final peak (response of the overall setup)
    
    % if p.fibre.b2b % in case of back to back measurement, we get impulse response, centered around max/peak value
    %     StartX = PeakX-10*p.rx.ovsFactor;
    %     StopX = PeakX+10*p.rx.ovsFactor;
    % end
    
    if StartX<1  % FIXME check for experimental part
        StartX = StartX + ovsCodeLength;%TODO = ovsCodeLength
    elseif StartX> ovsCodeLength
        StartX = StartX - ovsCodeLength;
    end
    
    if StopX<1
        StopX = StopX + ovsCodeLength;
    elseif StopX> ovsCodeLength
        StopX = StopX - ovsCodeLength;
    end
    %  StartX= 464490;
    %  StopX=60987;
    
    %Find the number of Rayleigh backscattering segments per code
    nbSegments = StopX-StartX +1;
    if nbSegments < 0
        nbSegments = StopX-StartX+ovsCodeLength+1;
    end
    
    EstimatedFiberLength = nbSegments*p.fibre.spatialRes/p.rx.ovsFactor; % (StopX-StartX)*p.fibre.spatialRes*p.tx.fSymb/p.rx.fSamp;
    fprintf('Rayleigh Start & Stop detection: PeakX:%d (ratioAfterOverBefore:%.1f) StartX:%d  StopX:%d  PeakY:%.2e  StartY:%.2e  StopY:%.2e  Estimated fiber length:%.1fm\n', PeakX, ratioIntens, StartX, StopX, PeakY, StartYLowRes, StopYLowRes, EstimatedFiberLength);
    
    if p.displ.detection
        if length(tmpTabLowRes) > 1
            %Display (low resolution, shifted)
            %         p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold off;
            %         plot(tmpTabLowRes./max(tmpTabLowRes),'g.-'); hold on;
            %         plot(ratioTmpTabLowRes./StartYLowRes,'b.-'); hold on;
            %         plot(1./ratioTmpTabLowRes.*StopYLowRes,'r.-'); hold on;
            %         xlabel(sprintf('Sample index (lowResFact=%d)',lowResFact)); ylabel('Norm. V^2');
            %         axis([1 length(tmpTabLowRes) 0 1.1]); grid on;
            %         title('Low resolution detection');
            
            %         %Display (full resolution, shifted)
            %         bufPeakStart = zeros(size(resDetShift)); bufPeakStart((StartXLowRes-1)*lowResFact + xminPrior - lowResFact -1 + tmp1X) = 1;
            %         bufPeakEnd = zeros(size(resDetShift)); bufPeakEnd((StopXLowRes-1)*lowResFact + xminPrior - 2*lowResFact - 1 + tmp2X) = 1;
            %
            %         p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold off;
            %         plot(resDetShift./max(resDetShift),'g.-'); hold on;
            %         plot(bufPeakStart,'b.-'); hold on;
            %         plot(bufPeakEnd,'r.-'); hold on;
            %         xlabel('Sample index'); ylabel('Norm. V^2');
            %         axis([1 length(resDetShift) 0 1.1]); grid on;
            %         title('Full resolution, shifted');
            
            %Display (full resolution, unshifted)
            bufPeakStart = zeros(size(resDetShift)); bufPeakStart(StartX) = 1;
            bufPeakEnd = zeros(size(resDetShift)); bufPeakEnd(StopX) = 1;
            
            p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold off;
            plot(resDet./max(resDet),'g.-'); hold on;
            plot(bufPeakStart,'b.-'); hold on;
            plot(bufPeakEnd,'r.-'); hold on;
            xlabel('Sample index'); ylabel('Norm. V^2');
            axis([1 ovsCodeLength 0 1]); grid on;
            title('Full resolution, unshifted');
            
            %         if p.fibre.b2b
            %             %Full block codes for comparison
            %             fullBlockResponses= reshape(absSumRxCorr(1+ovsCodeLength:(nbExpectedCodes-1)*ovsCodeLength), [ovsCodeLength,(nbExpectedCodes-2)]);
            %             ResponsesCentered = circshift(fullBlockResponses, ovsCodeLength/2- PeakX, 1);
            %             p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold off;
            %             plot(ResponsesCentered,'.-');
            %             xlabel('Sample index'); ylabel('Norm. V^2');
            %             axis([ovsCodeLength/2-2*nbSegments  ovsCodeLength/2+2*nbSegments 0 max(ResponsesCentered,[],'all')+0.5]);
            %             title('B2B: All impulse responses');
            %
            %             impulseresp = sum(ResponsesCentered,2);
            %             p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); hold off;
            %             plot(impulseresp./max(impulseresp),'.-');
            %             xlabel('Sample index'); ylabel('Norm. V^2');
            %             axis([ovsCodeLength/2-2*nbSegments  ovsCodeLength/2+2*nbSegments 0 1]);
            %             title('B2B: Mean impulse response');
            %         end
            
            
            %         p.displ.fIdx=p.displ.fIdx+1; figure(p.displ.fIdx); %hold off;
            %         plot(0.001*0.5*p.fibre.cFiber/p.rx.fSamp*resDet./max(resDet),'g.-'); hold on;
            %         plot(0.001*0.5*p.fibre.cFiber/p.rx.fSamp*ovsCodeLength, bufPeakStart,'b.-'); hold on;
            %         plot(0.001*0.5*p.fibre.cFiber/p.rx.fSamp*ovsCodeLength, bufPeakEnd,'r.-'); hold on;
            %         xlabel('Length (km)'); ylabel('Norm. V^2');
            %         axis([0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp)*1 0.001*max(0.5*p.fibre.cFiber/p.rx.fSamp)*ovsCodeLength 0 1]); grid on;
            %         title('Full resolution, unshifted');
            
            %      p.tx.ProbingMethod = 5;% 0 CODES ou 4 TRANSLATED SWEEPou 5 CAZAC
            %      r.on=1;
            %     [p,r]=SensingModelRayleigh_NoNoise(p,r);
            
        end %env
    end %p.displ.detection
    
elseif strcmpi(p.env,'model') %We know where the fiber response will be in simulation
    if (p.rx.offset_ratio-floor(p.rx.offset_ratio))
        StartX = 1+round(p.rx.offset_ratio-floor(p.rx.offset_ratio))*ovsCodeLength;
        StopX =  round(p.rx.offset_ratio-floor(p.rx.offset_ratio))*ovsCodeLength;
    else
        StartX = 1;
        StopX =  p.fibre.nbSegments;
    end
        
    nbSegments = StopX-StartX+1;
    EstimatedFiberLength = nbSegments*p.fibre.spatialRes/p.rx.ovsFactor; % (StopX-StartX)*p.fibre.spatialRes*p.tx.fSymb/p.rx.fSamp;
end

p.rx.nbOvsReflectors = nbSegments;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





