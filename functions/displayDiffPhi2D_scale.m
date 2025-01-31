function p=displayDiffPhi2D_scale(diffPhiTabSelect, linear_scale, exp_scale, max_dist_idx, max_time_idx, p)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays 2D Phase map (scaled)
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AdjReflectorsDist = (0.5*p.fibre.cFiber/p.rx.fSamp);
taxis = (0:max_time_idx-1).*p.tx.Tcode; %default max_time_idx is p.rx.nbDetectedCodes

%% Display intensity of backscatters & differential phase vs time of a selection of highest intensity backscatters
figure(p.displ.fIdx); hold off;
set(gcf,'units','points','position',[400,150,600,500])

pos1 = [0.07 0.1 0.14 0.85];
intensityplot = subplot('Position',pos1); hold off;
plot(p.rx.AvgIntensPerReflectorTab./max(p.rx.AvgIntensPerReflectorTab),(1:p.rx.nbReflectors)*AdjReflectorsDist,'LineWidth',2); hold on;
plot(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)./max(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)),p.displ.selectedIdxTab*AdjReflectorsDist,'ro');
axis([0 1 1 ceil(max_dist_idx*AdjReflectorsDist)]);
xlabel('RBS intensity'); ylabel('distance (m)')

pos2 = [0.24 0.1 0.6 0.85];
diffphaseplot = subplot('Position',pos2); hold off;
shifted_index = repmat(AdjReflectorsDist*p.displ.selectedIdxTab.', 1,p.rx.nbDetectedCodes); 
%plot((diffPhiTabSelect+repmat(AdjReflectorsDist*p.displ.selectedIdxTab.', 1,p.rx.nbDetectedCodes)).'); hold on;
scaled_diffPhi = linear_scale*abs(diffPhiTabSelect).^exp_scale; 
plot(taxis,(scaled_diffPhi+shifted_index).'); hold on;
plot(taxis,(-scaled_diffPhi+shifted_index).'); hold on;
axis([0 taxis(end) 1 ceil(max_dist_idx*AdjReflectorsDist)]);
title(sprintf('BP filtered (%d:%dHz) phase of selected reflectors vs time', p.rx.f_cutoff, p.rx.f_cutoff_end));
title(sprintf('Optical phase standard deviation as fct of fiber length \n(BP filter:%d-%dHz) MinPhi:%.2e MaxPhi:%.2e AvgPhi:%.2e', p.rx.f_cutoff,p.rx.f_cutoff_end, min(p.stdDiffPhiTabSelect(2:end-1)), max(p.stdDiffPhiTabSelect(2:end-1)), mean(p.stdDiffPhiTabSelect(2:end-1))));
xlabel('Time (s)'); ylabel('Phase (rad)');

pos3 = [0.85 0.1 0.14 0.85];
stdplot = subplot('Position',pos3); hold off;
plot(p.stdDiffPhiTabSelect,p.displ.selectedIdxTab*AdjReflectorsDist,'Color','red','LineWidth',2);
%axis([0 max(p.stdDiffPhiTabSelect) 1 ceil(max_dist_idx*AdjReflectorsDist)]);
axis([0 max(max(p.stdDiffPhiTabSelect(2:end-2)),1e-8) 1 ceil(max_dist_idx*AdjReflectorsDist)]);%CD modif 06/20 to avoid end of array artefact in MODEL mode
xlabel('Diff phase StDv'); %ylabel('distance (m)')

linkaxes([intensityplot, diffphaseplot, stdplot],'y'); %synchronize the individual axis limits across several figures

%% POLAR Display polarization intensity & SOP vs time of a selection of highest intensity backscatters
if p.displ.polar
    p.displ.fIdx=p.displ.fIdx+1;figure(p.displ.fIdx); hold off;
    % Variables :
    % p.pola.diffPolaS1,2,3
    % p.pola.S0_t
    
    set(gcf,'units','points','position',[400,150,600,500])
    %pos1 = [0.1 0.1 0.14 0.85];
    pos1 = [0.07 0.1 0.14 0.85];
    intensityplot = subplot('Position',pos1); hold off;
    plot(p.rx.AvgIntensPerReflectorTab./max(p.rx.AvgIntensPerReflectorTab),(1:p.rx.nbReflectors)*AdjReflectorsDist,'LineWidth',2); hold on;
    %plot(p.highestIntensThresh*ones(1,p.rx.nbReflectors),(1:p.rx.nbReflectors)*compLevelDist,'r-.');
    plot(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)./max(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)),p.displ.selectedIdxTab*AdjReflectorsDist,'ro');
    axis([0 1 1 ceil(max_dist_idx*AdjReflectorsDist)]);
    xlabel('RBS intensity'); ylabel(sprintf('Distance (m)'));
    %axis([0 1 1 ceil(max_dist_idx*compLevelDist*AdjReflectorsDist)]);
    %xlabel('RBS intensity'); ylabel('distance (m)')
    
    %compLevelDist = 0.1;%Compression
    %pos2 = [0.32 0.1 0.65 0.85];
    pos2 = [0.24 0.1 0.6 0.85];
    diffphaseplot = subplot('Position',pos2); hold off;
    
    plot(taxis, (p.pola.diffPolaS1+repmat(AdjReflectorsDist*p.displ.selectedIdxTab.', 1,p.rx.nbDetectedCodes)).'); hold on;
    plot(taxis, (p.pola.diffPolaS2+repmat(AdjReflectorsDist*p.displ.selectedIdxTab.', 1,p.rx.nbDetectedCodes)).'); hold on;
    plot(taxis, (p.pola.diffPolaS3+repmat(AdjReflectorsDist*p.displ.selectedIdxTab.', 1,p.rx.nbDetectedCodes)).'); hold on;
    axis([0 taxis(end) 1 ceil(max_dist_idx*AdjReflectorsDist)]);
    
    title(sprintf('HP filtered (%dHz) POLA of selected reflectors as fct of time', p.rx.f_cutoff));
    xlabel(sprintf('Time (s)')); ylabel('Phase (rad)');
    
    
    pos3 = [0.85 0.1 0.14 0.85];
    stdplot = subplot('Position',pos3); hold off;
%     plot(p.pola.stddiffPolaS1,p.displ.selectedIdxTab*AdjReflectorsDist,'LineWidth',2); hold on; 
%     plot(p.pola.stddiffPolaS2,p.displ.selectedIdxTab*AdjReflectorsDist,'LineWidth',2); hold on; 
%     plot(p.pola.stddiffPolaS3,p.displ.selectedIdxTab*AdjReflectorsDist,'LineWidth',2);
     plot(p.pola.stddiffPolaS1+p.pola.stddiffPolaS2+p.pola.stddiffPolaS3,p.displ.selectedIdxTab*AdjReflectorsDist,'LineWidth',2);

    %axis([0 max(p.stdDiffPhiTabSelect) 1 ceil(max_dist_idx*AdjReflectorsDist)]);
    %legend('S1','S2','S3'); 
    axis([0 max(max(p.stdDiffPhiTabSelect(2:end-2)),1e-8) 1 ceil(max_dist_idx*AdjReflectorsDist)]);%CD modif 06/20 to avoid end of array artefact in MODEL mode
    xlabel('S1+S2+S3 StDv'); %ylabel('distance (m)')
    
    
    
    linkaxes([intensityplot, diffphaseplot, stdplot],'y'); %synchronize the individual axis limits across several figures
    
end %polar

end