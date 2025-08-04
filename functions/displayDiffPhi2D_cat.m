function p=displayDiffPhi2D_cat(diffPhiTabSelect, compLevelDist, max_dist_idx, max_time_idx, p)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays 2D Phase map
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display intensity of backscatters & differential phase vs time of a selection of highest intensity backscatters
figure(p.displ.fIdx); hold off;

set(gcf,'units','points','position',[400,150,600,500])
pos1 = [0.1 0.1 0.14 0.85];
subplot('Position',pos1); hold off;
plot(p.rx.AvgIntensPerReflectorTab./max(p.rx.AvgIntensPerReflectorTab),(1:p.rx.nbReflectors)*compLevelDist,'LineWidth',2); hold on;
%plot(p.highestIntensThresh*ones(1,p.rx.nbOvsRayleighReflectors),(1:p.rx.nbOvsRayleighReflectors)*compLevelDist,'r-.');
plot(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)./max(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)),p.displ.selectedIdxTab*compLevelDist,'ro');
axis([0 1 1 ceil(max_dist_idx*compLevelDist)]);
xlabel('RBS intensity'); ylabel(sprintf('Reflector index (%.2fm spatial res.)',0.5*p.fibre.cFiber/p.rx.fSamp));

%compLevelDist = 0.1;%Compression
pos2 = [0.32 0.1 0.65 0.85];
subplot('Position',pos2); hold off;
plot(diffPhiTabSelect(1,:)+p.displ.selectedIdxTab(1)*compLevelDist, 'k--'); hold on;
for n=2:length(p.displ.selectedIdxTab)
    if p.stdDiffPhiTabSelect(n) < p.displ.stdPhiAlarmThres %Plot the phase std value on top of each curve
        plot(diffPhiTabSelect(n,:)+p.displ.selectedIdxTab(n)*compLevelDist); hold on;
        text(1+round(size(diffPhiTabSelect,2)*rand(1,1)),p.displ.selectedIdxTab(n),sprintf('%.2e',p.stdDiffPhiTabSelect(n)),'Color','black','FontSize',8);
    else
        plot(diffPhiTabSelect(n,:)+p.displ.selectedIdxTab(n)*compLevelDist, 'r','LineWidth',2); hold on;
        text(1+round(size(diffPhiTabSelect,2)*rand(1,1)),p.displ.selectedIdxTab(n),sprintf('%.2e',p.stdDiffPhiTabSelect(n)),'Color','red','FontSize',12);
    end
end
axis([0 max_time_idx 1 ceil(max_dist_idx*compLevelDist)]);
title(sprintf('HP filtered (%dHz) phase of selected reflectors as fct of time', p.rx.f_cutoff));
xlabel(sprintf('Time index (%.2es res., %.2fs overall)',p.tx.Tcode,p.tx.Tcode*max_time_idx)); ylabel('Phase (rad)');

%% POLAR Display polarization intensity & SOP vs time of a selection of highest intensity backscatters
if p.displ.polar
    p.displ.fIdx=p.displ.fIdx+1;figure(p.displ.fIdx); hold off;
    % Variables : 
    % p.pola.diffPolaS1,2,3
    % p.pola.S0_t
    
    set(gcf,'units','points','position',[400,150,600,500])
    pos1 = [0.1 0.1 0.14 0.85];
    subplot('Position',pos1); hold off;
    plot(p.rx.AvgIntensPerReflectorTab./max(p.rx.AvgIntensPerReflectorTab),(1:p.rx.nbOvsRayleighReflectors)*compLevelDist,'LineWidth',2); hold on;
    %plot(p.highestIntensThresh*ones(1,p.rx.nbOvsRayleighReflectors),(1:p.rx.nbOvsRayleighReflectors)*compLevelDist,'r-.');
    plot(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)./max(p.rx.AvgIntensPerReflectorTab(p.displ.selectedIdxTab)),p.displ.selectedIdxTab*compLevelDist,'ro');
    axis([0 1 1 ceil(max_dist_idx*compLevelDist)]);
    xlabel('RBS intensity'); ylabel(sprintf('Reflector index (%.2fm spatial res.)',0.5*p.fibre.cFiber/p.rx.fSamp));
    
    %compLevelDist = 0.1;%Compression
    pos2 = [0.32 0.1 0.65 0.85];
    subplot('Position',pos2); hold off;
    plot(p.pola.diffPolaS1(1,:)+p.displ.selectedIdxTab(1)*compLevelDist, 'k--'); hold on;
    for n=2:length(p.displ.selectedIdxTab)
        plot(p.pola.diffPolaS1(n,:)+p.displ.selectedIdxTab(n)*compLevelDist, '-'); hold on;
        text(1+round(size(p.pola.diffPolaS1,2)*rand(1,1)),p.displ.selectedIdxTab(n),sprintf('%.2e',p.pola.stddiffPolaS1(n)),'Color','black','FontSize',8);
        plot(p.pola.diffPolaS2(n,:)+p.displ.selectedIdxTab(n)*compLevelDist, '-.'); hold on;
        text(1+round(size(p.pola.diffPolaS2,2)*rand(1,1)),p.displ.selectedIdxTab(n),sprintf('%.2e',p.pola.stddiffPolaS2(n)),'Color','black','FontSize',8);
        plot(p.pola.diffPolaS3(n,:)+p.displ.selectedIdxTab(n)*compLevelDist,'--'); hold on;
        text(1+round(size(p.pola.diffPolaS3,2)*rand(1,1)),p.displ.selectedIdxTab(n),sprintf('%.2e',p.pola.stddiffPolaS3(n)),'Color','black','FontSize',8);
        
    end
    axis([0 max_time_idx 1 ceil(max_dist_idx*compLevelDist)]);
    title(sprintf('HP filtered (%dHz) POLA of selected reflectors as fct of time', p.rx.f_cutoff));
    xlabel(sprintf('Time index (%.2es res., %.2fs overall)',p.tx.Tcode,p.tx.Tcode*max_time_idx)); ylabel('Phase (rad)');
    
end