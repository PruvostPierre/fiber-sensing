%% DSP_Freq.m - Functions to study frequency

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MecaFunctions.m - Classe contenant les fonctions de méca
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef DSP_Freq

   methods(Static)
       
        function FFToutput = FFT_Wiwdow(signal, t1, t2, fc, DurationOfRecord,LogLog)
            
            
            s1=t1*length(signal)/DurationOfRecord
            s2=t2*length(signal)/DurationOfRecord
            
            % Assuming FFTDAS is your input signal
            segment = signal(s1:s2); % Extract the segment from sample 20000 to 42000
            
            % Compute the FFT of the segment
            Y = fft(segment);
            
            % Optional: Plot the magnitude of the FFT            
            L = length(segment);
            f = fc*(0:(L/2))/L; % Frequency axis
            
            % Compute the two-sided spectrum
            P2 = abs(Y/L);
            
            % Compute the single-sided spectrum
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            
            % Plot the Signal FFT in between t1 and t2 with cut frequency fc
            if strcmp(LogLog,'LogLog')

                figure;
                loglog(f, P1);
                title('Signal FFT in between t1 and t2 with cut frequency fc');
                xlabel('Frequency (f)');
                ylabel('|Signal(f)|');

            elseif strcmp(LogLog,'LogLin') 

                figure;
                semilogx(f, P1);
                title('Signal FFT in between t1 and t2 with cut frequency fc');
                xlabel('Frequency (f)');
                ylabel('|Signal(f)|');

            elseif strcmp(LogLog,'LinLog') 

                figure;
                semilogy(f, P1);
                title('Signal FFT in between t1 and t2 with cut frequency fc');
                xlabel('Frequency (f)');
                ylabel('|Signal(f)|');

            else                

                error('Invalid LogLog option. Choose "LogLog", "LogLin", or "LinLog".');

            end

            FFToutput = P1; % Store the output for further use
            
        end 

   end   
end

 