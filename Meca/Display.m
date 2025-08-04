%% Display.m - Functions to show the curves

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display.m - Classe contenant les fonctions d'affichage pour les essais mécaniques
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Display
    
    methods(Static)

        %% Fonction pour afficher plusieurs courbes d'un coup : LoadCell, AccMonoZ, Selected row, Corrélation glissante entre les deux 
        function quickdisplay(TestInProgress)
            % Traitement des données
            rowStudied = TestInProgress.rowStudied;
            selectedRow = TestInProgress.DAS.diffPhiTabSelect(rowStudied, :);
            LenDAS = length(selectedRow);
            LenMeca = length(TestInProgress.IMC.LoadCell50KN);

            % Interpolation
            x_original = linspace(1, LenDAS, LenDAS);
            x_target = linspace(1, LenDAS, LenMeca);
            selectedRow_resized = interp1(x_original, selectedRow, x_target, 'linear');

            % Affichage
            figure(1);
            tiledlayout('vertical')
            ax1 = nexttile; plot(TestInProgress.IMC.LoadCell50KN_time, TestInProgress.IMC.LoadCell50KN); title('LoadCell50KN');
            ax2 = nexttile; plot(TestInProgress.IMC.Acc_Mono_Z_time, TestInProgress.IMC.Acc_Mono_Z); title('AccMono_Z');
            ax3 = nexttile; plot(TestInProgress.IMC.LoadCell50KN_time, selectedRow_resized); title('SelectedRow interpolated');

            % Corrélation glissante
            t=0.2;
            windowSize = t*length(TestInProgress.IMC.Acc_Mono_Z)/TestInProgress.timeRecord;
            step = 10;
            nPoints = length(TestInProgress.IMC.LoadCell50KN);
            nSteps = floor((nPoints - windowSize)/step) + 1;

            correlationCoeff = zeros(1, nSteps);
            time_corr = zeros(1, nSteps);

            for i = 1:nSteps
                idx = (i-1)*step + 1 : (i-1)*step + windowSize;
                segment1 = TestInProgress.IMC.LoadCell50KN(idx);
                segment2 = selectedRow_resized(idx);
                r = corr(segment1(:), segment2(:));
                correlationCoeff(i) = r;
                time_corr(i) = mean(TestInProgress.IMC.LoadCell50KN_time(idx));
            end

            ax4 = nexttile;
            plot(time_corr, correlationCoeff)
            title('Corrélation glissante Meca/DAS')
            xlabel('Temps (s)')
            ylabel('Coeff. de corrélation')
            ylim([-1 1])

            linkaxes([ax1, ax2, ax3, ax4], 'x');

            % Figure supplémentaire
            figure(2)
            plot(TestInProgress.IMC.LoadCell50KN_time, TestInProgress.IMC.LoadCell50KN)
            title('Impact mesuré en N')

            % Appel à la fonction externe computeCoherency
            % computeCoherency(TestInProgress);
        end
        
        %% Fonction pour afficher uen courbe issue de l'IMC et une autre issue du DAS en les remettant à la même taille
        function DisplayComp(time_signal1, signal1, signal2, signal1Titre, signal2Titre)
            % Vérifier que signal1 et time_signal1 ont la même taille
            if length(time_signal1) ~= length(signal1)
                error('time_signal1 et signal1 doivent avoir la même longueur.');
            end
        
            % Redimensionner signal2 pour matcher signal1
            Lensignal1 = length(signal1);
            Lensignal2 = length(signal2);
        
            % Interpolation de signal2 → taille de signal1
            x_original = linspace(1, Lensignal2, Lensignal2);
            x_target = linspace(1, Lensignal2, Lensignal1);
            signal2_resized = interp1(x_original, signal2, x_target, 'linear');
        
            % Affichage
            figure;
            tiledlayout('vertical');
                ax1 = nexttile; plot(time_signal1, signal1); title(signal1Titre);
                ax2 = nexttile; plot(time_signal1, signal2_resized); title(signal2Titre);
        
            linkaxes([ax1, ax2], 'x');

            disp("Taille time_signal1 : " + length(time_signal1));
            disp("Taille signal1 : " + length(signal1));
            disp("Taille signal2 : " + length(signal2));
            disp("Taille signal2_resized : " + length(signal2_resized));

        end

         %% Fonction pour afficher une courbe avec la distance affichée (DAS)
         function DisplayDistAjust(signal, Lchannel)
             
         end

         function PLOT_Lin(signal, titre, gridOnOff, Xaxis, Yaxis)
              
             figure; 
             plot(signal);
             title(titre);
             grid(gridOnOff);
             xlabel(Xaxis);
             ylabel(Yaxis);

         end


    end
end
