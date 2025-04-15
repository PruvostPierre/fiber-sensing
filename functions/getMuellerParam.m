function p = getMuellerParam(p,JonesMatTab3D)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GETMUELLERPARAM get polarization information from Jones matrix
% Stokes parameters evolution, differential polar are stored il p.pola
%
% Authors: 
% Original code by S. Guerrier, C. Dorize & E. Awwad - 2022
% Modified version by E. Awwad - 2024 elie.awwad@telecom-paris.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RxCorrXX = p.HiTab(1,1:p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes);
RxCorrXY = p.HiTab(2,1:p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes);

%Get evolution of beta polarization parameter
% Number of Jones matrices
% numJonesMatrices = size(p.HiTab, 2);
% p.fibre.lengthSegment = floor(numJonesMatrices/p.fibre.nbSegments);
% % Initialize a 3D array to store the Jones matrices
% p.pola.JonesMatrices = zeros(2, 2, floor(numJonesMatrices/p.fibre.lengthSegment));
% 
% % Loop through each column and extract the Jones matrix
% for n = 1:numJonesMatrices
%     p.pola.JonesMatrices(:, :, n) = [p.HiTab(1, n), p.HiTab(2, n); p.HiTab(3, n), p.HiTab(4, n)];
%     %[U, S, V] = svd(p.pola.JonesMatrices(:, :, n));
%     %p.pola.JonesMatrices(:, :, n) = U*V';
%     %singvalues = svd(p.pola.JonesMatrices(:, :, n));
%     %det_H = det(p.pola.JonesMatrices(:, :, n));
%     %exp_2jargp_m = det_H/(singvalues(1)^2);
%     %argp_m = angle(exp_2jargp_m)/2;
%     %argp_m = angle(det_H)/2;
%     %p_m = singvalues(1)*exp(1j*argp_m);
%     %disp(p_m);
%     p.pola.JonesMatrices(:, :, n) = p.pola.JonesMatrices(:, :, n)*(1/sqrt(det(p.pola.JonesMatrices(:, :, n))));%normalisation
% end
% 
% %list of polarization angle variations
% p.pola.polaAngle = zeros(1,floor(numJonesMatrices/p.fibre.lengthSegment));
% p.pola.diffPolaAngle = zeros(1,floor(numJonesMatrices/p.fibre.lengthSegment)-1);
% p.pola.ProdJonesMatrices = zeros(2,2,floor(numJonesMatrices/p.fibre.lengthSegment));
% p.pola.eigenvalues = zeros(2, floor(numJonesMatrices/p.fibre.lengthSegment));
% 
% %calculate H_roundtrip_z_dagger*H_roundtrip_z' for for each pair of z, z'separated by one fiber segment
% %calculate the eigenvalues of the product of the Jones matrices for each pair of z, z'separated by one fiber segment
% k=0;
% for n = 1:p.fibre.lengthSegment:numJonesMatrices-p.fibre.lengthSegment
%     k=k+1;
%     p.pola.ProdJonesMatrices(:,:,k) = p.pola.JonesMatrices(:,:,n)'*p.pola.JonesMatrices(:,:,n+p.fibre.lengthSegment);
% end
% 
% for i=1:size(p.pola.ProdJonesMatrices,3)
%     p.pola.eigenvalues(:,i)= eig(p.pola.ProdJonesMatrices(:,:,i));
%     arg_eigenvalues = angle(p.pola.eigenvalues(:,i));
%     pos_angle = zeros(1,2);
%     for j=1:2
%         pos_angle(j) = getAnglePlusMinusPiOver2(arg_eigenvalues(j));
%     end
%     p.pola.polaAngle(i) = (1/2)*pos_angle(pos_angle>0); %pola angle
%     %if p.pola.polaAngle(i) > pi/2-0.01
%     %    p.pola.polaAngle(i) = p.pola.polaAngle(i) - pi/2;
%     %end
% end
% 
% for i=1:size(p.pola.ProdJonesMatrices,3)-1
%     p.pola.diffPolaAngle(i) = p.pola.polaAngle(i+1)-p.pola.polaAngle(i);
% end

A = [1 0 0 1;1 0 0 -1;0 1 1 0;0 -1j 1j 0];
    time_line = p.rx.nbOvsSelectedReflectors*p.rx.nbDetectedCodes;% size(JonesMatTab3D,3); %time_line/sizeTab = nbDetectedCodes
    Mi = zeros(4,4,time_line);
    p.pola.S0_t = zeros(1,time_line);
    for n=1:time_line
        Mi(:,:,n) = A*kron(JonesMatTab3D(:,:,n), conj(JonesMatTab3D(:,:,n)))/A; %Jones to Mueller transformation
        p.pola.S0_t(n) = Mi(1,1,n); %see evolution of total intensity
        Mi(:,:,n) = (1/Mi(1,1,n))*Mi(:,:,n); %normalisation (s0 = 1)
    end
    % Evolution of Stokes parameters
    Stokes_t = zeros(3, time_line); %S0 to 1, here S=(s1,s2,s3)
    p.pola.S1_t = zeros(1,time_line); p.pola.S2_t = zeros(1,time_line); p.pola.S3_t = zeros(1,time_line);
    for n=1:time_line
        St = Mi(:,:,n)*[1;0;0;-1];%circular left hand; [1;0;1;0] ;%linearly 45deg %[1 1 0 0]linearly polarized input (1/sqrt(2))*
        Stokes_t(:,n) = St(2:end);%Mi(2:4,2:4,n)*[1;0;0];
        p.pola.S1_t(n)=St(2); p.pola.S2_t(n)=St(3); p.pola.S3_t(n)=St(4);
    end
    
    p.pola.S0_t = single(reshape(p.pola.S0_t, [p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]));
    
    % Differential pola / Standard deviation of SOP
    p.pola.diffPolaS1 = p.pola.S1_t; p.pola.diffPolaS1 = [Stokes_t(1,1), diff(p.pola.diffPolaS1)];
    p.pola.diffPolaS1 = single(reshape(p.pola.diffPolaS1,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]));%Reshape vector to get reflectors per row and codes per column
    p.pola.diffPolaS2 = p.pola.S2_t; p.pola.diffPolaS2 = [Stokes_t(2,1), diff(p.pola.diffPolaS2)];
    p.pola.diffPolaS2 = single(reshape(p.pola.diffPolaS2,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]));%Reshape vector to get reflectors per row and codes per column
    p.pola.diffPolaS3 = p.pola.S3_t; p.pola.diffPolaS3 = [Stokes_t(3,1), diff(p.pola.diffPolaS3)];
    p.pola.diffPolaS3 = single(reshape(p.pola.diffPolaS3,[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]));%Reshape vector to get reflectors per row and codes per column
    for n=1:p.rx.nbOvsSelectedReflectors
        p.pola.diffPolaS1(n,:) = real(ifft(p.rx.spWeightingTab.*fft(p.pola.diffPolaS1(n,:)))); %with HP filtering
        p.pola.diffPolaS2(n,:) = real(ifft(p.rx.spWeightingTab.*fft(p.pola.diffPolaS2(n,:)))); %with HP filtering
        p.pola.diffPolaS3(n,:) = real(ifft(p.rx.spWeightingTab.*fft(p.pola.diffPolaS3(n,:)))); %with HP filtering
    end
    %p.pola.globalDiff;
    p.pola.stddiffPolaS1 = single(std(p.pola.diffPolaS1.'));
    p.pola.stddiffPolaS2 = single(std(p.pola.diffPolaS2.'));
    p.pola.stddiffPolaS3 = single(std(p.pola.diffPolaS3.'));
    
    % ellipticity / fast axis / retardance
    p.pola.DOP = 0.5*reshape(sqrt(sum(Stokes_t.^2, 1)),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]); %sqrt(Stokes_t(1,:).^2+Stokes_t(2,:).^2+Stokes_t(3,:).^2)
    p.pola.eta = reshape(0.5*asin(Stokes_t(3,:)),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]); %sin(2eta) = s3/s0
    p.pola.psy = reshape(0.5*atan(Stokes_t(2,:)./Stokes_t(1,:)),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]); % tan(2psy)=s2/s1
    p.pola.ellipticity = reshape(tan(p.pola.eta),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]); %ellipticity e=minoraxis/majoraxis=tan(eta)
    p.pola.deltaPhi = reshape(angle(RxCorrXX)-angle(RxCorrXY),[p.rx.nbOvsSelectedReflectors,p.rx.nbDetectedCodes]);
    p.pola.deltaPhi = getAnglePlusMinusPiOver2((p.pola.deltaPhi).').';%Get differential phase between the selected reflectors
    %p.pola.deltaPhi = (unwrap(2*p.pola.deltaPhi.')/2).';%Unwrap to get rid of PI phase jumps between reflectors
    %p.pola.deltaPhi = std(p.pola.deltaPhi.');%Differential phase std of the selected backscatters after HP filtering (Edges in time ignored) %(:,round(0.1*p.rx.nbDetectedCodes):round(0.9*p.rx.nbDetectedCodes))
    
    p.pola.delta = p.pola.deltaPhi./sqrt(tan(2*p.pola.eta).^2+1);%see Rogers 1981 AppendixA
    p.pola.rho = 0.5*p.pola.delta.*tan(2*p.pola.eta); %Rogers 1981 Appendix A
end

