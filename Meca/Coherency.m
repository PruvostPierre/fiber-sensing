%% Coherency.m - Functions to show the curves

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display.m - Classe contenant les fonctions d'affichage pour les essais mécaniques
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Coherency
    methods(Static)

        %% Calcule le coefficient de cohérence entre segments de la fibre optique (doi:10.26443/seismica.v4i1.1488)
        function C = CohSeg(strain, twin, TestInProgress)

            [channelNUM, echNUM] = size(strain);
            %Découpage de la matrice de contrainte pour garder une taille 
            % (channelNUM x twin) :
            
            echCUT=echNUM*(twin/TestInProgress.timeRecord)
            strain_window=strain(:, 1:echCUT);

            %Boucle for pour les segments (x_channel)
            for i = 1:channelNUM - 1
                SW1 = strain_window(i, :);
                SW2 = strain_window(i + 1, :);
                
                num = dot(SW1, SW2);  % Produit scalaire
                denom = norm(SW1) * norm(SW2);  % Produit des normes
                
                if denom ~= 0
                    C(i) = abs(dot(SW1, SW2) / (norm(SW1) * norm(SW2))); %(Pour avoir C entre 0 et 1)
                else
                    C(i) = 0;  % Évite la division par zéro
                end
            end
            
            C(channelNUM) = C(channelNUM-1); % To keep the same vector size
        end


        %% Permet de trouver la meilleure cohérence entre 2 signaux par décallage circulaire
        function [bestShiftSigned, xcorrMAX] = CircCorr(signal1, signal2)

            if length(signal1) ~= length(signal2)
                error('Impossible : signal1 et signal2 doivent avoir la même longueur.');
            end
        
            L = length(signal1);
            maxCorr = -inf;
            bestShift = 0;
        
            for shift = 0:L-1
                shifted = circshift(signal2, shift);
                corrVal = max(xcorr(signal1, shifted, 'normalized'));
        
                if corrVal > maxCorr
                    maxCorr = corrVal;
                    bestShift = shift;
                end
            end
        
            % Convertir le décalage en valeur signée
            % Si bestShift > L/2, cela correspond à un décalage négatif
            if bestShift <= L/2
                bestShiftSigned = bestShift;
            else
                bestShiftSigned = bestShift - L;
            end
        
            xcorrMAX = maxCorr;
        
            % Affichage
            fprintf('Décalage optimal : %d échantillons\n', bestShiftSigned);
            fprintf('Corrélation maximale : %.4f\n', xcorrMAX);
        
        end




        function [y,varargout]=xcorAlign(x)
            % Aligns data according to cross correlation maximum. If only 1 output is 
            % specified, then dataout (a matrix the same size as data) is returned.
            %
            % USAGE:
            % [dataout,lagout] = xcorAlign(data);
            %
            % INPUT:
            % data: a matrix of data vectors stored column-wise.
            %
            % OUTPUT: 
            % y:        A matrix size(data) back, circulary shifted according to
            %           maximum correlation with x(:,1)
            % lagout:   Optional.  A vector of lag indices.  These give the amount the 
            %           column vectors of dataout are shifted with respect to the first 
            %           column of data.
            %-----------------------------------------------------------------------
            % Latest Edit: 02.March.2007
            % Joshua D Carmichael
            % josh.carmichael@gmail.com
            %
            % Edit Log
            %-----------------------------------------------------------------------
            
            [M,N] = size(x);
            y     = x;
            lagout= zeros(N,1);
            
            %make infs zeros for computing purpose
            x(isinf(x)) = 0;
            
            %align against the first column vector
            ftx1 = fft(x(:,2:end));
            ftx2 = fft(x(:,1));
            
            for k=2:N
            
            xc          = ifft(conj(ftx1(:,k-1)).*ftx2);
            [m,i]       = max(abs(xc));
            y(:,k)      = circshift(x(:,k),i);
            lagout(k)   = i;
            
            end
            
            %asign output values
            if(nargout==2)
            varargout{1}=lagout;
            end
        end
    end
end

